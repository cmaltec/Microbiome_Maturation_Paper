#!/usr/bin/env Rscript
#________________________________________
# - File Name : 03_growth_curve.R
# - Purpose : Cohort-mean Gompertz-with-Gaussian-dip growth-curve fit
# - Created By : Christian Maltecca
# - e-mail: cmaltec@ncsu.edu
# - Note : Organized for this repository from the author's original analysis
#          scripts with the help of Claude (Anthropic).
#      / \       _-'
#     /|  \-''- _ /
# ' { |          \
#   /             \
#  /      "o.  |o }
#     |           \
#                 ',
#      \_         __\
#        ''-_    \.//
#          / '-____'
#         /
#        '
#     -'
# ╔╗╔╔═╗╔═╗╦ ╦
# ║║║║  ╚═╗║ ║
# ╝╚╝╚═╝╚═╝╚═╝
# Mon May 25 10:09:00 2026 ------------------------------

setwd(here::here())
suppressPackageStartupMessages({
  library(arrow)
  library(minpack.lm)
})

# --- Load data ---------------------------------------------------------------

wide <- read_parquet("data/processed/mag_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
feat_cols <- setdiff(names(wide), "animal_id")
mat <- as.matrix(wide[, feat_cols])
rownames(mat) <- animal_ids

feature_names <- unique(sapply(strsplit(feat_cols, "_(?=W[0-9]+$)", perl = TRUE),
                               function(x) x[1]))
tp_labels <- c("W3", "W6", "W9", "W12", "W15", "W17", "W21", "W24")
week_map <- c(W3 = 3, W6 = 6, W9 = 9, W12 = 12, W15 = 15,
              W17 = 17, W21 = 21, W24 = 24)

ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
line_lookup <- setNames(ped$line, ped$animal_id)

# --- Compute BC-based maturity per sample ------------------------------------

w24_profiles <- list()
for (a in animal_ids) {
  cols <- paste0(feature_names, "_W24")
  v <- mat[a, cols]
  if (!all(is.na(v))) w24_profiles[[a]] <- v
}
w24_centroid <- colMeans(do.call(rbind, w24_profiles))
w24_centroid_prop <- w24_centroid / sum(w24_centroid)

bc_all <- list()
for (a in animal_ids) {
  for (tp in tp_labels) {
    cols <- paste0(feature_names, "_", tp)
    v <- mat[a, cols]
    if (all(is.na(v))) next
    p <- v / sum(v)
    bc_all[[length(bc_all) + 1]] <- data.frame(
      a = a, tp = tp, w = week_map[[tp]],
      maturity = 1 - sum(abs(p - w24_centroid_prop)) / 2,
      stringsAsFactors = FALSE)
  }
}
bc_df <- do.call(rbind, bc_all)

# Population asymptotes for 2-parameter models
pop_A_w17 <- mean(bc_df$maturity[bc_df$tp == "W17"], na.rm = TRUE)

# --- Population-level Gompertz + Gaussian dip fit ----------------------------
# Models the W21 dip explicitly:
#   y(t) = A * exp(-exp(-k*(t-ti))) - C * exp(-(t-mu)^2 / (2*sigma^2))
# 5p: mu fixed at 21 (biological anchor); 6p form (mu free) buys ~0.005 R^2.
# Fit on pooled (animal, week) observations. Population estimates of
# C, sigma, mu are reused as fixed values for the per-animal constrained fit.

fit_gompertz_dip_population <- function(weeks, maturity, mu_fixed = 21) {
  out <- list(A = NA, k = NA, ti = NA, C = NA, sigma = NA, mu = mu_fixed,
              r2 = NA, rmse = NA, aic = NA, converged = FALSE,
              method = "none", n_pts = length(weeks))
  if (length(weeks) < 5) return(out)

  best_rss <- Inf
  df <- data.frame(y = maturity, t = weeks)

  A_starts <- c(0.85, 0.9, 0.95)
  k_starts <- c(0.05, 0.1, 0.2)
  ti_starts <- c(0, 3, 6)
  C_starts <- c(0.1, 0.2, 0.3)
  sigma_starts <- c(1.5, 2.5, 4.0)

  for (A0 in A_starts) for (k0 in k_starts) for (ti0 in ti_starts) {
    for (C0 in C_starts) for (s0 in sigma_starts) {
      fit <- tryCatch(suppressWarnings(
        nlsLM(y ~ A * exp(-exp(-k * (t - ti))) -
                  C * exp(-(t - mu_fixed)^2 / (2 * sigma^2)),
              data = df,
              start = list(A = A0, k = k0, ti = ti0, C = C0, sigma = s0),
              lower = c(A = 0.3, k = 0.001, ti = -10, C = 0,    sigma = 0.3),
              upper = c(A = 1.5, k = 5.0,   ti = 30,  C = 0.6, sigma = 10),
              control = nls.lm.control(maxiter = 1000, ftol = 1e-10, ptol = 1e-10))
      ), error = function(e) NULL)
      if (!is.null(fit)) {
        rss <- sum(residuals(fit)^2)
        if (rss < best_rss) {
          ss_tot <- sum((maturity - mean(maturity))^2)
          out$r2 <- if (ss_tot > 0) 1 - rss / ss_tot else NA
          out$rmse <- sqrt(mean(residuals(fit)^2))
          out$aic <- tryCatch(AIC(fit), error = function(e) NA)
          out$converged <- TRUE; out$method <- "nlsLM_5p_mu_fixed"
          cc <- coef(fit)
          out$A <- cc[["A"]]; out$k <- cc[["k"]]; out$ti <- cc[["ti"]]
          out$C <- cc[["C"]]; out$sigma <- cc[["sigma"]]
          best_rss <- rss
        }
      }
    }
  }
  out
}

cat("Fitting population-level Gompertz + Gaussian dip (5p, mu=21 fixed)...\n")
# Fit on the per-week population mean trajectory (8 points). This recovers
# the *typical* dip shape (C, sigma) cleanly. Pooled per-observation fit
# would dilute C and sigma toward 0 because animals that don't dip pull the
# fitted dip magnitude down.
pop_mean_df <- aggregate(maturity ~ w, data = bc_df, FUN = mean)
pop_dip <- fit_gompertz_dip_population(pop_mean_df$w, pop_mean_df$maturity,
                                       mu_fixed = 21)
cat(sprintf("  A=%.4f  k=%.4f  ti=%.4f  C=%.4f  sigma=%.4f  mu=%.0f\n",
            pop_dip$A, pop_dip$k, pop_dip$ti, pop_dip$C, pop_dip$sigma, pop_dip$mu))
cat(sprintf("  R^2=%.4f  RMSE=%.4f  n=%d (population mean trajectory)\n",
            pop_dip$r2, pop_dip$rmse, pop_dip$n_pts))

# Also report fit quality on pooled data with these population parameters,
# for honest comparison with vanilla Gompertz 3p full.
fit_yhat <- pop_dip$A * exp(-exp(-pop_dip$k * (bc_df$w - pop_dip$ti))) -
  pop_dip$C * exp(-(bc_df$w - pop_dip$mu)^2 / (2 * pop_dip$sigma^2))
pooled_rss <- sum((bc_df$maturity - fit_yhat)^2)
pooled_tss <- sum((bc_df$maturity - mean(bc_df$maturity))^2)
pooled_r2 <- 1 - pooled_rss / pooled_tss
pooled_rmse <- sqrt(pooled_rss / nrow(bc_df))
cat(sprintf("  Pooled (n=%d) R^2=%.4f  RMSE=%.4f using population params\n",
            nrow(bc_df), pooled_r2, pooled_rmse))

pop_dip_df <- data.frame(
  param = c("A", "k", "ti", "C", "sigma", "mu",
            "r2_mean_traj", "rmse_mean_traj", "aic_mean_traj", "n_pts_mean",
            "r2_pooled", "rmse_pooled", "n_pts_pooled"),
  value = c(pop_dip$A, pop_dip$k, pop_dip$ti, pop_dip$C, pop_dip$sigma,
            pop_dip$mu,
            pop_dip$r2, pop_dip$rmse, pop_dip$aic, pop_dip$n_pts,
            pooled_r2, pooled_rmse, nrow(bc_df)),
  stringsAsFactors = FALSE)
write.csv(pop_dip_df,
          "results/metrics/maturation_growth_curve_population_dip.csv",
          row.names = FALSE)
cat("  Saved: maturation_growth_curve_population_dip.csv\n")

# --- Per-animal constrained Gompertz + dip fit (3 free params) ---------------
# Full model is fit, but C, sigma, mu are fixed at the population estimates.
# Only A, k, ti vary per animal. Feasible with 4-6 obs per animal.

fit_gompertz_dip_constrained <- function(weeks, maturity,
                                         C_fix, sigma_fix, mu_fix) {
  out <- list(A = NA, k = NA, ti = NA, r2 = NA, rmse = NA, aic = NA,
              converged = FALSE, method = "none", n_pts = length(weeks))
  if (length(weeks) < 3 || sd(maturity) < 1e-6) return(out)

  best_rss <- Inf
  df <- data.frame(y = maturity, t = weeks,
                   C_fix = C_fix, sigma_fix = sigma_fix, mu_fix = mu_fix)

  A_starts <- c(max(maturity) + 0.05, max(maturity) + 0.2, 0.85, 1.0)
  k_starts <- c(0.05, 0.1, 0.2, 0.5, 1.0)
  ti_starts <- c(0, 3, 6, 10, 15)

  for (A0 in A_starts) for (k0 in k_starts) for (ti0 in ti_starts) {
    fit <- tryCatch(suppressWarnings(
      nlsLM(y ~ A * exp(-exp(-k * (t - ti))) -
                C_fix * exp(-(t - mu_fix)^2 / (2 * sigma_fix^2)),
            data = df,
            start = list(A = A0, k = k0, ti = ti0),
            lower = c(A = 0.3, k = 0.001, ti = -10),
            upper = c(A = 1.5, k = 5.0,   ti = 30),
            control = nls.lm.control(maxiter = 1000, ftol = 1e-10, ptol = 1e-10))
    ), error = function(e) NULL)
    if (!is.null(fit)) {
      rss <- sum(residuals(fit)^2)
      if (rss < best_rss) {
        ss_tot <- sum((maturity - mean(maturity))^2)
        out$r2 <- if (ss_tot > 0) 1 - rss / ss_tot else NA
        out$rmse <- sqrt(mean(residuals(fit)^2))
        out$aic <- tryCatch(AIC(fit), error = function(e) NA)
        out$converged <- TRUE; out$method <- "nlsLM_3p_constrained"
        cc <- coef(fit)
        out$A <- cc[["A"]]; out$k <- cc[["k"]]; out$ti <- cc[["ti"]]
        best_rss <- rss
      }
    }
  }
  out
}

# --- Robust Gompertz fitter --------------------------------------------------
# Tries nlsLM (Levenberg-Marquardt), nls port (bounded), and SSgompertz
# over a grid of starting values. Keeps the fit with lowest RSS.

fit_gompertz_robust <- function(weeks, maturity, fixed_A = NULL) {
  out <- list(k = NA, ti = NA, A = NA, r2 = NA, rmse = NA, aic = NA,
              converged = FALSE, method = "none", n_pts = length(weeks))
  if (length(weeks) < 3 || sd(maturity) < 1e-6) return(out)

  best_rss <- Inf
  df <- data.frame(y = maturity, t = weeks)

  A_starts <- if (is.null(fixed_A)) c(max(maturity) + 0.05, max(maturity) + 0.2, 0.8, 1.0)
  k_starts <- c(0.05, 0.1, 0.2, 0.5, 1.0)
  ti_starts <- c(0, 3, 6, 10, 15)

  try_fit <- function(fit, method_name) {
    if (is.null(fit)) return()
    rss <- sum(residuals(fit)^2)
    if (rss < best_rss) {
      ss_tot <- sum((maturity - mean(maturity))^2)
      out$r2 <<- if (ss_tot > 0) 1 - rss / ss_tot else NA
      out$rmse <<- sqrt(mean(residuals(fit)^2))
      out$aic <<- tryCatch(AIC(fit), error = function(e) NA)
      out$converged <<- TRUE
      out$method <<- method_name
      best_rss <<- rss
      cc <- coef(fit)
      out$k <<- cc[["k"]]
      out$ti <<- cc[["ti"]]
      out$A <<- if (is.null(fixed_A)) cc[["A"]] else fixed_A
    }
  }

  for (k0 in k_starts) {
    for (ti0 in ti_starts) {
      if (is.null(fixed_A)) {
        for (A0 in A_starts) {
          fit <- tryCatch(suppressWarnings(
            nlsLM(y ~ A * exp(-exp(-k * (t - ti))), data = df,
                  start = list(A = A0, k = k0, ti = ti0),
                  lower = c(A = 0.3, k = 0.001, ti = -10),
                  upper = c(A = 1.5, k = 5.0, ti = 30),
                  control = nls.lm.control(maxiter = 1000, ftol = 1e-10, ptol = 1e-10))
          ), error = function(e) NULL)
          try_fit(fit, "nlsLM_3p")

          fit2 <- tryCatch(suppressWarnings(
            nls(y ~ A * exp(-exp(-k * (t - ti))), data = df,
                start = list(A = A0, k = k0, ti = ti0), algorithm = "port",
                lower = c(A = 0.3, k = 0.001, ti = -10),
                upper = c(A = 1.5, k = 5.0, ti = 30),
                control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE))
          ), error = function(e) NULL)
          try_fit(fit2, "nls_port_3p")
        }
      } else {
        df$A_fix <- fixed_A
        fit <- tryCatch(suppressWarnings(
          nlsLM(y ~ A_fix * exp(-exp(-k * (t - ti))), data = df,
                start = list(k = k0, ti = ti0),
                lower = c(k = 0.001, ti = -10), upper = c(k = 5.0, ti = 30),
                control = nls.lm.control(maxiter = 1000, ftol = 1e-10, ptol = 1e-10))
        ), error = function(e) NULL)
        try_fit(fit, "nlsLM_2p")

        fit2 <- tryCatch(suppressWarnings(
          nls(y ~ A_fix * exp(-exp(-k * (t - ti))), data = df,
              start = list(k = k0, ti = ti0), algorithm = "port",
              lower = c(k = 0.001, ti = -10), upper = c(k = 5.0, ti = 30),
              control = nls.control(maxiter = 1000, tol = 1e-6, warnOnly = TRUE))
        ), error = function(e) NULL)
        try_fit(fit2, "nls_port_2p")
      }
    }
  }

  # SSgompertz self-starting fallback (3p only, different parameterization)
  if (is.null(fixed_A)) {
    fit3 <- tryCatch(suppressWarnings(
      nls(y ~ SSgompertz(t, Asym, b2, b3), data = df,
          control = nls.control(maxiter = 500, warnOnly = TRUE))
    ), error = function(e) NULL)
    if (!is.null(fit3)) {
      cc3 <- coef(fit3)
      Asym <- cc3[["Asym"]]; b2 <- cc3[["b2"]]; b3 <- cc3[["b3"]]
      # Convert: SSgompertz y = Asym * exp(-b2 * b3^t) to standard parameterization
      if (b3 > 0 && b3 < 1 && b2 > 0) {
        rss3 <- sum(residuals(fit3)^2)
        if (rss3 < best_rss) {
          ss_tot <- sum((maturity - mean(maturity))^2)
          out$r2 <- if (ss_tot > 0) 1 - rss3 / ss_tot else NA
          out$rmse <- sqrt(mean(residuals(fit3)^2))
          out$aic <- tryCatch(AIC(fit3), error = function(e) NA)
          out$converged <- TRUE
          out$method <- "SSgompertz"
          out$A <- Asym
          out$k <- -log(b3)
          out$ti <- log(b2) / (-log(b3))
        }
      }
    }
  }

  out
}

# --- Robust logistic fitter --------------------------------------------------

fit_logistic_robust <- function(weeks, maturity, fixed_L = NULL) {
  out <- list(k = NA, t0 = NA, L = NA, r2 = NA, rmse = NA, aic = NA,
              converged = FALSE, method = "none", n_pts = length(weeks))
  if (length(weeks) < 3 || sd(maturity) < 1e-6) return(out)

  best_rss <- Inf
  df <- data.frame(y = maturity, t = weeks)

  L_starts <- if (is.null(fixed_L)) c(max(maturity) + 0.05, 0.8, 1.0)
  k_starts <- c(0.05, 0.1, 0.3, 0.5, 1.0)
  t0_starts <- c(0, 5, 10)

  for (k0 in k_starts) {
    for (t00 in t0_starts) {
      if (is.null(fixed_L)) {
        for (L0 in L_starts) {
          fit <- tryCatch(suppressWarnings(
            nlsLM(y ~ L / (1 + exp(-k * (t - t0))), data = df,
                  start = list(L = L0, k = k0, t0 = t00),
                  lower = c(L = 0.3, k = 0.001, t0 = -10),
                  upper = c(L = 1.5, k = 5.0, t0 = 30),
                  control = nls.lm.control(maxiter = 1000))
          ), error = function(e) NULL)
          if (!is.null(fit)) {
            rss <- sum(residuals(fit)^2)
            if (rss < best_rss) {
              ss_tot <- sum((maturity - mean(maturity))^2)
              out$r2 <- if (ss_tot > 0) 1 - rss / ss_tot else NA
              out$rmse <- sqrt(mean(residuals(fit)^2))
              out$aic <- tryCatch(AIC(fit), error = function(e) NA)
              out$converged <- TRUE; out$method <- "nlsLM_3p"
              cc <- coef(fit)
              out$k <- cc[["k"]]; out$t0 <- cc[["t0"]]; out$L <- cc[["L"]]
              best_rss <- rss
            }
          }
        }
      } else {
        df$L_fix <- fixed_L
        fit <- tryCatch(suppressWarnings(
          nlsLM(y ~ L_fix / (1 + exp(-k * (t - t0))), data = df,
                start = list(k = k0, t0 = t00),
                lower = c(k = 0.001, t0 = -10), upper = c(k = 5.0, t0 = 30),
                control = nls.lm.control(maxiter = 1000))
        ), error = function(e) NULL)
        if (!is.null(fit)) {
          rss <- sum(residuals(fit)^2)
          if (rss < best_rss) {
            ss_tot <- sum((maturity - mean(maturity))^2)
            out$r2 <- if (ss_tot > 0) 1 - rss / ss_tot else NA
            out$rmse <- sqrt(mean(residuals(fit)^2))
            out$aic <- tryCatch(AIC(fit), error = function(e) NA)
            out$converged <- TRUE; out$method <- "nlsLM_2p"
            cc <- coef(fit)
            out$k <- cc[["k"]]; out$t0 <- cc[["t0"]]; out$L <- fixed_L
            best_rss <- rss
          }
        }
      }
    }
  }
  out
}

# --- Fit all animals ---------------------------------------------------------

cat("Fitting robust growth curves...\n")
results <- list()

for (ai in seq_along(animal_ids)) {
  a <- animal_ids[ai]
  adf <- bc_df[bc_df$a == a, ]
  if (nrow(adf) < 3) next

  full <- adf[order(adf$w), ]
  mono <- full[full$w <= 17, ]

  s <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)
  sb <- function(x) if (is.null(x) || length(x) == 0) FALSE else as.logical(x)
  sc <- function(x) if (is.null(x) || length(x) == 0) "none" else as.character(x)

  g3m <- fit_gompertz_robust(mono$w, mono$maturity)
  g2m <- fit_gompertz_robust(mono$w, mono$maturity, fixed_A = pop_A_w17)
  g3f <- fit_gompertz_robust(full$w, full$maturity)
  l3m <- fit_logistic_robust(mono$w, mono$maturity)
  l2m <- fit_logistic_robust(mono$w, mono$maturity, fixed_L = pop_A_w17)
  gdc <- fit_gompertz_dip_constrained(full$w, full$maturity,
                                      C_fix = pop_dip$C,
                                      sigma_fix = pop_dip$sigma,
                                      mu_fix = pop_dip$mu)

  lin <- lm(maturity ~ w, data = mono)

  results[[a]] <- data.frame(
    animal_id = a, line = line_lookup[[a]],
    n_obs_full = nrow(full), n_obs_mono = nrow(mono),
    g3m_k = s(g3m$k), g3m_ti = s(g3m$ti), g3m_A = s(g3m$A),
    g3m_r2 = s(g3m$r2), g3m_rmse = s(g3m$rmse), g3m_aic = s(g3m$aic),
    g3m_conv = sb(g3m$converged), g3m_method = sc(g3m$method),
    g2m_k = s(g2m$k), g2m_ti = s(g2m$ti), g2m_A = s(g2m$A),
    g2m_r2 = s(g2m$r2), g2m_rmse = s(g2m$rmse), g2m_aic = s(g2m$aic),
    g2m_conv = sb(g2m$converged), g2m_method = sc(g2m$method),
    g3f_k = s(g3f$k), g3f_ti = s(g3f$ti), g3f_A = s(g3f$A),
    g3f_r2 = s(g3f$r2), g3f_rmse = s(g3f$rmse), g3f_aic = s(g3f$aic),
    g3f_conv = sb(g3f$converged), g3f_method = sc(g3f$method),
    l3m_k = s(l3m$k), l3m_t0 = s(l3m$t0), l3m_L = s(l3m$L),
    l3m_r2 = s(l3m$r2), l3m_rmse = s(l3m$rmse), l3m_aic = s(l3m$aic),
    l3m_conv = sb(l3m$converged),
    l2m_k = s(l2m$k), l2m_t0 = s(l2m$t0), l2m_L = s(l2m$L),
    l2m_r2 = s(l2m$r2), l2m_rmse = s(l2m$rmse), l2m_aic = s(l2m$aic),
    l2m_conv = sb(l2m$converged),
    gd_constr_A = s(gdc$A), gd_constr_k = s(gdc$k), gd_constr_ti = s(gdc$ti),
    gd_constr_r2 = s(gdc$r2), gd_constr_rmse = s(gdc$rmse),
    gd_constr_aic = s(gdc$aic), gd_constr_conv = sb(gdc$converged),
    lin_slope = coef(lin)[2], lin_r2 = summary(lin)$r.squared,
    lin_rmse = sqrt(mean(residuals(lin)^2)),
    stringsAsFactors = FALSE
  )

  if (ai %% 20 == 0) cat(sprintf("  %d/%d\n", ai, length(animal_ids)))
}

df <- do.call(rbind, results)
cat(sprintf("\nFitted %d animals\n", nrow(df)))

# --- Model comparison summary ------------------------------------------------

cat("\n=== Model comparison ===\n")
cat(sprintf("  %-25s %8s %8s %10s %10s %10s\n",
            "Model", "Conv%", "n_conv", "med_R2", "med_RMSE", "med_AIC"))

models <- list(
  "Gompertz 3p W3-17"       = list(c = "g3m_conv", r = "g3m_r2", m = "g3m_rmse", a = "g3m_aic"),
  "Gompertz 2p W3-17"       = list(c = "g2m_conv", r = "g2m_r2", m = "g2m_rmse", a = "g2m_aic"),
  "Gompertz 3p full"        = list(c = "g3f_conv", r = "g3f_r2", m = "g3f_rmse", a = "g3f_aic"),
  "Gompertz+dip constr full"= list(c = "gd_constr_conv", r = "gd_constr_r2",
                                   m = "gd_constr_rmse", a = "gd_constr_aic"),
  "Logistic 3p W3-17"       = list(c = "l3m_conv", r = "l3m_r2", m = "l3m_rmse", a = "l3m_aic"),
  "Logistic 2p W3-17"       = list(c = "l2m_conv", r = "l2m_r2", m = "l2m_rmse", a = "l2m_aic"),
  "Linear W3-17"            = list(c = NULL,       r = "lin_r2",  m = "lin_rmse", a = NULL))

for (nm in names(models)) {
  mo <- models[[nm]]
  if (!is.null(mo$c)) {
    ok <- df[[mo$c]]
    cat(sprintf("  %-25s %7.0f%% %8d %10.4f %10.4f %10.1f\n",
                nm, 100 * mean(ok, na.rm = TRUE), sum(ok, na.rm = TRUE),
                median(df[[mo$r]][ok], na.rm = TRUE),
                median(df[[mo$m]][ok], na.rm = TRUE),
                median(df[[mo$a]][ok], na.rm = TRUE)))
  } else {
    cat(sprintf("  %-25s %7.0f%% %8d %10.4f %10.4f %10s\n",
                nm, 100, nrow(df),
                median(df[[mo$r]], na.rm = TRUE),
                median(df[[mo$m]], na.rm = TRUE), "---"))
  }
}

# --- Fitting method breakdown ------------------------------------------------

cat("\n=== Gompertz 3p W3-17: fitting method used ===\n")
print(table(df$g3m_method))

# --- SEN vs TOL comparison ---------------------------------------------------

cat("\n=== SEN vs TOL (converged animals) ===\n")
model_specs <- list(
  list(name = "Gompertz 3p W3-17", conv = "g3m_conv",
       params = c("g3m_k", "g3m_ti", "g3m_A", "g3m_r2")),
  list(name = "Gompertz 2p W3-17", conv = "g2m_conv",
       params = c("g2m_k", "g2m_ti", "g2m_r2")),
  list(name = "Gompertz+dip constr full", conv = "gd_constr_conv",
       params = c("gd_constr_A", "gd_constr_k", "gd_constr_ti", "gd_constr_r2")),
  list(name = "Logistic 3p W3-17", conv = "l3m_conv",
       params = c("l3m_k", "l3m_t0", "l3m_L", "l3m_r2")),
  list(name = "Linear W3-17", conv = NULL,
       params = c("lin_slope", "lin_r2")))

for (spec in model_specs) {
  cat(sprintf("\n  --- %s ---\n", spec$name))
  ok <- if (!is.null(spec$conv)) df[[spec$conv]] else rep(TRUE, nrow(df))
  for (p in spec$params) {
    vs <- df[[p]][df$line == "SEN" & ok]
    vt <- df[[p]][df$line == "TOL" & ok]
    if (sum(!is.na(vs)) < 3 || sum(!is.na(vt)) < 3) next
    tt <- t.test(vs, vt)
    cat(sprintf("    %-12s SEN=%7.4f(n=%d) TOL=%7.4f(n=%d) p=%.4f%s\n",
                p, mean(vs, na.rm = TRUE), sum(!is.na(vs)),
                mean(vt, na.rm = TRUE), sum(!is.na(vt)),
                tt$p.value,
                ifelse(tt$p.value < 0.05, " *", ifelse(tt$p.value < 0.1, " .", ""))))
  }
}

# --- Save --------------------------------------------------------------------

write.csv(df, "results/metrics/maturation_growth_curve_fits.csv", row.names = FALSE)
cat(sprintf("\nSaved: maturation_growth_curve_fits.csv (%d animals, %d cols)\n",
            nrow(df), ncol(df)))
