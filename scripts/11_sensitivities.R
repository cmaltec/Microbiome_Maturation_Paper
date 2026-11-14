#!/usr/bin/env Rscript
#________________________________________
# - File Name : 11_sensitivities.R
# - Purpose : Sensitivity refits and cross-estimator concordance
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
# Mon May 25 10:11:00 2026 ------------------------------

setwd(here::here())

# --- Load h² panels ---------------------------------------------------------

obs <- read.csv("results/metrics/maturation_h2_animal_model.csv",
                stringsAsFactors = FALSE)
imp <- read.csv("results/metrics/maturation_h2_animal_model_imputed.csv",
                stringsAsFactors = FALSE)
lin <- read.csv("results/metrics/maturation_h2_line_sensitivity.csv",
                stringsAsFactors = FALSE)
w6  <- read.csv("results/metrics/maturation_h2_w6_scalars.csv",
                stringsAsFactors = FALSE)
# The W6+ file names the Aitchison slope aitchison_to_adult_slope; the other
# three estimator files use aitchison_slope.
w6$metric[w6$metric == "aitchison_to_adult_slope"] <- "aitchison_slope"

# Map W6 metric names back to canonical names used in the obs/imp panels.
# (13_w3_artifact_audit.R already maps internally; the saved file uses the
#  full-window canonical names, so no further mapping is needed.)

# --- Align on a common metric set ------------------------------------------

m_obs <- obs[, c("metric", "h2")];   colnames(m_obs)[2] <- "h2_obs"
m_imp <- imp[, c("metric", "h2")];   colnames(m_imp)[2] <- "h2_imp"
m_lin <- lin[, c("metric", "h2_with_line")]; colnames(m_lin)[2] <- "h2_within_line"
m_w6  <- w6[,  c("metric", "h2_w6")]

mat <- Reduce(function(x, y) merge(x, y, by = "metric", all = TRUE),
              list(m_obs, m_imp, m_lin, m_w6))

# Drop rows where every estimator is NA
mat <- mat[rowSums(!is.na(mat[, -1])) > 0, ]
cat(sprintf("Concordance matrix: %d metrics × %d estimators\n",
            nrow(mat), ncol(mat) - 1))
print(mat)

write.csv(mat, "results/metrics/maturation_h2_concordance.csv", row.names = FALSE)

# --- Metric set in all four estimators (full pairwise) ---------------------

mat_full <- mat[complete.cases(mat[, -1]), ]
cat(sprintf("\nMetrics with all 4 estimators non-NA: %d\n", nrow(mat_full)))
print(mat_full)

# --- Pairwise Pearson r and Lin's CCC ---------------------------------------

lin_ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(c(r = NA, ccc = NA, bias = NA, n = length(x)))
  r  <- cor(x, y)
  mx <- mean(x); my <- mean(y); sx <- var(x); sy <- var(y)
  ccc <- 2 * r * sqrt(sx) * sqrt(sy) /
         (sx + sy + (mx - my)^2)
  c(r = r, ccc = ccc, bias = mx - my, n = length(x))
}

ests <- c("h2_obs", "h2_imp", "h2_within_line", "h2_w6")
n_e <- length(ests)
pair_out <- data.frame()
cat("\n=== Pairwise concordance (h² across estimators) ===\n")
for (i in seq_len(n_e - 1)) {
  for (j in (i + 1):n_e) {
    a <- ests[i]; b <- ests[j]
    res <- lin_ccc(mat[[a]], mat[[b]])
    pair_out <- rbind(pair_out, data.frame(
      est_A = a, est_B = b,
      n = res["n"], pearson_r = res["r"],
      lin_CCC = res["ccc"], mean_bias = res["bias"],
      stringsAsFactors = FALSE))
    cat(sprintf("  %-18s vs %-18s  n=%d  r=%+.3f  CCC=%+.3f  bias=%+.3f\n",
                a, b, res["n"], res["r"], res["ccc"], res["bias"]))
  }
}
write.csv(pair_out, "results/metrics/maturation_h2_concordance_pairwise.csv",
          row.names = FALSE)

# --- Kendall's W (coefficient of concordance) ------------------------------
# Uses metrics common to ALL 4 estimators. Each estimator ranks the metrics
# by h². W = 1 → perfect agreement; W = 0 → no agreement above chance.

kendall_W <- function(M) {
  # M: n_metrics × n_estimators
  R <- apply(M, 2, rank, ties.method = "average")
  Ri <- rowSums(R)
  Rbar <- mean(Ri)
  S <- sum((Ri - Rbar)^2)
  m <- ncol(M); n <- nrow(M)
  W <- 12 * S / (m^2 * (n^3 - n))
  chi2 <- m * (n - 1) * W
  p_chi <- pchisq(chi2, df = n - 1, lower.tail = FALSE)
  list(W = W, chi2 = chi2, df = n - 1, p_chi2 = p_chi, m = m, n = n)
}

if (nrow(mat_full) >= 3) {
  M <- as.matrix(mat_full[, ests])
  kw <- kendall_W(M)
  cat(sprintf("\nKendall's W (n=%d metrics × m=%d estimators):  W = %.3f  ",
              kw$n, kw$m, kw$W))
  cat(sprintf("χ²(%d) = %.2f  p_χ² = %.4g\n", kw$df, kw$chi2, kw$p_chi2))

  # Permutation p-value (within-estimator shuffle preserves estimator's
  # rank distribution but breaks metric identity)
  set.seed(42); B <- 9999
  W_perm <- numeric(B)
  for (b in seq_len(B)) {
    Mp <- apply(M, 2, sample)
    W_perm[b] <- kendall_W(Mp)$W
  }
  p_perm <- mean(W_perm >= kw$W)
  cat(sprintf("Permutation p (B=%d): p = %.4f\n", B, p_perm))

  # --- Subset W on the W3-robust panel only ---
  # Metrics with all four estimators that are heritable (BH q < 0.10) on W6+.
  w3_robust <- c("maz_mean", "intra_bc_var", "volatility_mean",
                 "bc_to_adult_slope", "aitchison_slope",
                 "shannon_slope", "t70", "t75")
  rob_idx <- mat_full$metric %in% w3_robust
  if (sum(rob_idx) >= 3) {
    Mr <- as.matrix(mat_full[rob_idx, ests])
    kwr <- kendall_W(Mr)
    cat(sprintf("Kendall's W on W3-robust panel only (n=%d, m=%d): W = %.3f  ",
                kwr$n, kwr$m, kwr$W))
    cat(sprintf("χ²(%d) = %.2f  p = %.4g\n",
                kwr$df, kwr$chi2, kwr$p_chi2))
    # permutation p for the robust panel (same scheme as the full panel)
    set.seed(42)
    Wr_perm <- replicate(B, kendall_W(apply(Mr, 2, sample))$W)
    p_perm_r <- mean(Wr_perm >= kwr$W)
    cat(sprintf("Permutation p (B=%d): p = %.4f\n", B, p_perm_r))
  } else {
    p_perm_r <- NA
    kwr <- list(W = NA, chi2 = NA, df = NA, p_chi2 = NA, n = sum(rob_idx))
  }

  summary_out <- data.frame(
    panel = c("all_metrics", "w3_robust"),
    n_metrics = c(kw$n, kwr$n),
    n_estimators = c(kw$m, kw$m),
    Kendall_W = c(kw$W, kwr$W),
    chi2 = c(kw$chi2, kwr$chi2),
    df = c(kw$df, kwr$df),
    p_chi2 = c(kw$p_chi2, kwr$p_chi2),
    p_perm_9999 = c(p_perm, p_perm_r),
    stringsAsFactors = FALSE
  )
  write.csv(summary_out,
            "results/metrics/maturation_h2_concordance_summary.csv",
            row.names = FALSE)
} else {
  cat("\nNot enough metrics for Kendall's W (need ≥ 3 with all 4 estimators)\n")
}

cat("\nOutputs:\n")
cat("  results/metrics/maturation_h2_concordance.csv\n")
cat("  results/metrics/maturation_h2_concordance_pairwise.csv\n")
cat("  results/metrics/maturation_h2_concordance_summary.csv\n")
cat("Done.\n")

# PART 2: Sampling-cohort (age in days at W3) covariate sensitivity
setwd(here::here())
suppressPackageStartupMessages(library(asreml))

ped3      <- read.csv("data/metadata/pedigree_3gen_3col.csv",   stringsAsFactors = FALSE)
Ainv      <- readRDS("data/processed/A_inverse_3gen.rds")
all_ids   <- ped3$id

metrics   <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                      stringsAsFactors = FALSE)
gc_fits   <- read.csv("results/metrics/maturation_growth_curve_fits.csv",
                      stringsAsFactors = FALSE)
ped_orig  <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

gc_cols <- c("animal_id", "g3m_k", "g3m_ti", "g3m_A", "l3m_k", "l3m_t0", "l3m_L")
metrics <- merge(metrics, gc_fits[, gc_cols], by = "animal_id", all.x = TRUE)
metrics <- merge(metrics,
                 ped_orig[, c("animal_id", "litter_id", "age_days_w3")],
                 by = "animal_id", all.x = TRUE)

metrics$animal_id <- factor(metrics$animal_id, levels = all_ids)
metrics$litter_id <- factor(metrics$litter_id)

cat(sprintf("Data: %d animals, age_days_w3 range %d-%d (mean %.1f, SD %.2f)\n",
            nrow(metrics),
            min(metrics$age_days_w3, na.rm = TRUE),
            max(metrics$age_days_w3, na.rm = TRUE),
            mean(metrics$age_days_w3, na.rm = TRUE),
            sd(metrics$age_days_w3, na.rm = TRUE)))

metric_defs <- list(
  list(name = "maz_mean",          class = "1-Age prediction",    col = "maz_mean"),
  list(name = "maz_slope",         class = "1-Age prediction",    col = "maz_slope"),
  list(name = "bc_to_adult_mean",  class = "2-Distance to adult", col = "bc_to_adult_mean"),
  list(name = "bc_to_adult_slope", class = "2-Distance to adult", col = "bc_to_adult_slope"),
  list(name = "aitchison_slope",   class = "2-Distance to adult", col = "aitchison_to_adult_slope"),
  list(name = "rel_maturity_w24",  class = "2-Distance to adult", col = "rel_maturity_w24"),
  list(name = "gompertz_k",        class = "3-Growth curve",      col = "g3m_k"),
  list(name = "gompertz_ti",       class = "3-Growth curve",      col = "g3m_ti"),
  list(name = "gompertz_A",        class = "3-Growth curve",      col = "g3m_A"),
  list(name = "logistic_k",        class = "3-Growth curve",      col = "l3m_k"),
  list(name = "logistic_t0",       class = "3-Growth curve",      col = "l3m_t0"),
  list(name = "volatility_mean",   class = "4-Stability",         col = "volatility_mean"),
  list(name = "volatility_slope",  class = "4-Stability",         col = "volatility_slope"),
  list(name = "intra_bc_var",      class = "4-Stability",         col = "intra_bc_var"),
  list(name = "shannon_cv",        class = "4-Stability",         col = "shannon_cv"),
  list(name = "shannon_mean",      class = "5-Diversity",         col = "shannon_mean"),
  list(name = "shannon_w24",       class = "5-Diversity",         col = "shannon_w24"),
  list(name = "shannon_slope",     class = "5-Diversity",         col = "shannon_slope"),
  list(name = "BC_FPC1",           class = "6-FPCA",              col = "BC_FPC1"),
  list(name = "BC_FPC2",           class = "6-FPCA",              col = "BC_FPC2"),
  list(name = "BC_FPC3",           class = "6-FPCA",              col = "BC_FPC3"),
  list(name = "shan_FPC1",         class = "6-FPCA",              col = "shan_FPC1"),
  list(name = "shan_FPC2",         class = "6-FPCA",              col = "shan_FPC2"),
  list(name = "shan_FPC3",         class = "6-FPCA",              col = "shan_FPC3"))

fit_h2 <- function(df_fit, fixed_formula) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed    = fixed_formula,
           random   = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units),
           data     = df_fit, na.action = na.method(y = "include"),
           trace = FALSE, maxit = 100)
  )), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  h2_vp <- tryCatch(vpredict(fit, h2 ~ V2 / (V1 + V2 + V3)),
                    error = function(e) NULL)
  list(h2    = if (!is.null(h2_vp)) h2_vp$Estimate[1] else NA_real_,
       h2_se = if (!is.null(h2_vp)) h2_vp$SE[1]       else NA_real_,
       fit = fit)
}

results <- data.frame()
for (md in metric_defs) {
  y_raw <- metrics[[md$col]]
  ok    <- !is.na(y_raw) & !is.na(metrics$age_days_w3)
  if (sum(ok) < 30) {
    cat(sprintf("  %-22s SKIP (n=%d)\n", md$name, sum(ok))); next
  }

  df_fit <- data.frame(
    y           = as.numeric(scale(y_raw)),
    animal_id   = metrics$animal_id,
    litter_id   = metrics$litter_id,
    age_days_w3 = metrics$age_days_w3)

  base <- fit_h2(df_fit, y ~ 1)
  date <- fit_h2(df_fit, y ~ 1 + age_days_w3)
  if (is.null(base) || is.null(date)) {
    cat(sprintf("  %-22s ERROR\n", md$name)); next
  }

  # Wald-style p on age_days_w3
  age_p <- NA_real_
  wald <- tryCatch(wald(date$fit), error = function(e) NULL)
  if (!is.null(wald)) {
    rn <- rownames(wald)
    idx <- grep("age_days_w3", rn)
    if (length(idx) == 1 && "Pr(Chisq)" %in% colnames(wald)) {
      age_p <- wald[idx, "Pr(Chisq)"]
    } else if (length(idx) == 1 && "Pr(>F)" %in% colnames(wald)) {
      age_p <- wald[idx, "Pr(>F)"]
    }
  }

  delta <- date$h2 - base$h2
  flag  <- ifelse(abs(delta) >= 0.05, "*", "")
  cat(sprintf("  %-22s h2_base=%.3f  h2_date=%.3f  Δ=%+.3f%s  age_p=%.3f\n",
              md$name, base$h2, date$h2, delta, flag,
              ifelse(is.na(age_p), NA_real_, age_p)))

  results <- rbind(results, data.frame(
    metric        = md$name,
    class         = md$class,
    n             = sum(ok),
    h2_baseline   = base$h2,
    h2_se_base    = base$h2_se,
    h2_date       = date$h2,
    h2_se_date    = date$h2_se,
    delta_h2      = delta,
    age_p         = age_p,
    stringsAsFactors = FALSE))
}

write.csv(results,
          "results/metrics/maturation_h2_sampling_date_sensitivity.csv",
          row.names = FALSE)
cat("\nSaved: results/metrics/maturation_h2_sampling_date_sensitivity.csv\n")
cat(sprintf("\nMetrics with |Δh²| ≥ 0.05: %d / %d\n",
            sum(abs(results$delta_h2) >= 0.05, na.rm = TRUE),
            nrow(results)))
cat(sprintf("Metrics with age_days_w3 fixed-effect p < 0.05: %d / %d\n",
            sum(results$age_p < 0.05, na.rm = TRUE),
            sum(!is.na(results$age_p))))
print(results, row.names = FALSE, digits = 3)

# PART 3: W3-W17 pre-perturbation sensitivity refit
setwd(here::here())
suppressPackageStartupMessages({
  library(arrow); library(asreml); library(dplyr)
})

# --- Load data (mirrors 01) --------------------------------------------------
wide <- read_parquet("data/processed/mag_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
feat_cols <- setdiff(names(wide), "animal_id")
mat <- as.matrix(wide[, feat_cols])
rownames(mat) <- animal_ids
col_parts <- strsplit(feat_cols, "_(?=W[0-9]+$)", perl = TRUE)
feature_names <- unique(sapply(col_parts, function(x) x[1]))
tp_labels <- unique(sapply(col_parts, function(x) x[2]))
week_map <- c(W3=3, W6=6, W9=9, W12=12, W15=15, W17=17, W21=21, W24=24)

ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
line_lookup <- setNames(ped$line, ped$animal_id)
curves <- read.csv("results/metrics/maturation_curves_per_sample.csv",
                   stringsAsFactors = FALSE)
maz_sample <- read.csv("results/metrics/maz_per_sample.csv",
                       stringsAsFactors = FALSE)

cat(sprintf("Loaded: %d animals, %d features\n", length(animal_ids), length(feature_names)))

# --- Helper (from 01) --------------------------------------------------------
get_profile <- function(animal, tp) {
  cols <- paste0(feature_names, "_", tp)
  v <- mat[animal, cols]
  if (all(is.na(v))) return(NULL)
  v
}
clr <- function(x) {
  x[x == 0] <- 0.5
  lx <- log(x); lx - mean(lx)
}

# --- W24 adult centroid (uses full data; this is a population reference) ----
w24_profiles <- list()
for (a in animal_ids) {
  p <- get_profile(a, "W24")
  if (!is.null(p)) w24_profiles[[a]] <- p
}
w24_mat <- do.call(rbind, w24_profiles)
w24_centroid <- colMeans(w24_mat)
w24_centroid_prop <- w24_centroid / sum(w24_centroid)

# --- Compute per-animal scalars on a given window ---------------------------
compute_metrics <- function(window_max_week) {
  cat(sprintf("\nComputing per-animal metrics for W3-W%d window...\n", window_max_week))
  results <- list()
  for (ai in seq_along(animal_ids)) {
    a <- animal_ids[ai]
    line <- line_lookup[[a]]
    # Observed timepoints within the window
    obs_tps <- character(0)
    for (tp in tp_labels) {
      if (week_map[[tp]] <= window_max_week &&
          !is.na(mat[a, paste0(feature_names[1], "_", tp)])) {
        obs_tps <- c(obs_tps, tp)
      }
    }
    obs_weeks <- week_map[obs_tps]
    n_obs <- length(obs_tps)
    if (n_obs < 3) next

    # -- Class 1: MAZ aggregated from per-sample MAZ on this window --
    maz_rows <- maz_sample[maz_sample$animal_id == a & maz_sample$week <= window_max_week, ]
    maz_mean  <- if (nrow(maz_rows) >= 2) mean(maz_rows$maz, na.rm = TRUE) else NA
    maz_slope <- if (nrow(maz_rows) >= 2) coef(lm(maz ~ week, data = maz_rows))[2] else NA

    # -- Class 2: Distance to adult --
    bc_vals <- numeric(n_obs); aitchison_vals <- numeric(n_obs)
    for (i in seq_along(obs_tps)) {
      p <- get_profile(a, obs_tps[i])
      p_prop <- p / sum(p)
      bc_vals[i] <- sum(abs(p_prop - w24_centroid_prop)) / 2
      aitchison_vals[i] <- sqrt(sum((clr(p) - clr(w24_centroid))^2))
    }
    bc_to_adult_mean <- mean(bc_vals)
    bc_to_adult_slope <- coef(lm(bc_vals ~ obs_weeks))[2]
    aitchison_to_adult_mean <- mean(aitchison_vals)
    aitchison_to_adult_slope <- coef(lm(aitchison_vals ~ obs_weeks))[2]

    # -- Class 4: Stability / volatility (canonical: BC between consecutive profiles) --
    consec_bc <- numeric(n_obs - 1)
    for (i in seq_len(n_obs - 1)) {
      p1 <- get_profile(a, obs_tps[i])
      p2 <- get_profile(a, obs_tps[i + 1])
      consec_bc[i] <- sum(abs(p1 / sum(p1) - p2 / sum(p2))) / 2
    }
    volatility_mean <- mean(consec_bc)
    volatility_slope <- if (length(consec_bc) >= 3) {
      mid_weeks <- (obs_weeks[-n_obs] + obs_weeks[-1]) / 2
      coef(lm(consec_bc ~ mid_weeks))[2]
    } else NA

    # Intra-individual BC variance
    all_profiles <- list()
    for (tp in obs_tps) {
      p <- get_profile(a, tp); all_profiles[[tp]] <- p / sum(p)
    }
    pairwise_bc <- c()
    for (i in seq_len(n_obs - 1)) for (j in (i + 1):n_obs) {
      pairwise_bc <- c(pairwise_bc,
        sum(abs(all_profiles[[obs_tps[i]]] - all_profiles[[obs_tps[j]]])) / 2)
    }
    intra_bc_var <- var(pairwise_bc)

    # -- Class 5: Shannon --
    animal_curves <- curves[curves$animal_id == a & curves$tp %in% obs_tps, ]
    shannon_vals <- animal_curves$shannon[match(obs_tps, animal_curves$tp)]
    shannon_cv    <- sd(shannon_vals, na.rm = TRUE) / mean(shannon_vals, na.rm = TRUE)
    shannon_slope <- coef(lm(shannon_vals ~ obs_weeks))[2]
    shannon_mean  <- mean(shannon_vals, na.rm = TRUE)

    results[[ai]] <- data.frame(
      animal_id = a, line = line, n_obs = n_obs,
      maz_mean = maz_mean, maz_slope = maz_slope,
      bc_to_adult_mean = bc_to_adult_mean,
      bc_to_adult_slope = bc_to_adult_slope,
      aitchison_to_adult_mean = aitchison_to_adult_mean,
      aitchison_to_adult_slope = aitchison_to_adult_slope,
      volatility_mean = volatility_mean,
      volatility_slope = volatility_slope,
      intra_bc_var = intra_bc_var,
      shannon_cv = shannon_cv,
      shannon_mean = shannon_mean,
      shannon_slope = shannon_slope,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, results)
}

full_scalars <- compute_metrics(24)
w17_scalars  <- compute_metrics(17)

write.csv(w17_scalars, "results/metrics/maturation_all_metrics_per_animal_w17.csv",
          row.names = FALSE)
cat(sprintf("\nSaved: results/metrics/maturation_all_metrics_per_animal_w17.csv (%d animals)\n",
            nrow(w17_scalars)))

# --- Fit univariate animal model (mirrors 03) --------------------------------
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
ped_ids <- attr(Ainv, "rowNames")
ped_orig <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

fit_h2 <- function(scalars, trait) {
  y_vec <- scalars[[trait]]
  if (sum(!is.na(y_vec)) < 30 || sd(y_vec, na.rm = TRUE) < 1e-9) return(NULL)
  df <- data.frame(
    y = as.numeric(scale(y_vec)),
    animal_id = factor(as.character(scalars$animal_id), levels = ped_ids),
    litter_id = factor(ped_orig$litter_id[match(scalars$animal_id, ped_orig$animal_id)]))
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed = y ~ 1,
           random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = df,
           na.action = na.method(y = "include"),
           trace = FALSE, maxit = 200))),
    error = function(e) NULL)
  if (is.null(fit) || !fit$converge) return(NULL)
  vc <- summary(fit)$varcomp; rn <- rownames(vc)
  s2_a <- max(vc[grep("vm.animal_id", rn)[1], "component"], 0)
  s2_c <- max(vc[which(rn == "litter_id"), "component"], 0)
  s2_e <- vc[grep("units", rn)[1], "component"]
  tot <- s2_a + s2_c + s2_e
  if (tot <= 0) return(NULL)
  h2_vp <- tryCatch(vpredict(fit, h2 ~ V2/(V1+V2+V3)), error = function(e) NULL)
  list(n = sum(!is.na(y_vec)),
       h2 = s2_a/tot,
       h2_se = if (!is.null(h2_vp)) h2_vp$SE[1] else NA,
       c2 = s2_c/tot)
}

traits <- c("maz_mean", "maz_slope",
            "bc_to_adult_mean", "bc_to_adult_slope",
            "aitchison_to_adult_mean", "aitchison_to_adult_slope",
            "volatility_mean", "volatility_slope", "intra_bc_var",
            "shannon_cv", "shannon_mean", "shannon_slope")

cat("\n=== W3-W24 vs W3-W17 h² comparison (3-gen Ainv, canonical scalars) ===\n")
cat(sprintf("%-25s | %-22s | %-22s | %s\n",
            "trait", "full W3-W24", "W3-W17 (no W21)", "Δh²"))
cat(strrep("-", 95), "\n", sep = "")
results <- data.frame()
for (tr in traits) {
  rf <- fit_h2(full_scalars, tr)
  rw <- fit_h2(w17_scalars,  tr)
  if (is.null(rf) || is.null(rw)) {
    cat(sprintf("%-25s | FIT FAILED\n", tr)); next
  }
  d <- rw$h2 - rf$h2
  flag <- if (abs(d) >= 0.05) " *" else ""
  cat(sprintf("%-25s | h²=%.3f (SE %.3f) n=%d | h²=%.3f (SE %.3f) n=%d | %+.3f%s\n",
              tr, rf$h2, rf$h2_se, rf$n, rw$h2, rw$h2_se, rw$n, d, flag))
  results <- rbind(results, data.frame(
    trait = tr,
    h2_full = rf$h2, h2_se_full = rf$h2_se, n_full = rf$n,
    h2_w17  = rw$h2, h2_se_w17  = rw$h2_se, n_w17  = rw$n,
    delta_h2 = d, delta_abs = abs(d)))
}
write.csv(results, "results/metrics/maturation_h2_full_vs_w17.csv", row.names = FALSE)
cat(sprintf("\nSaved: results/metrics/maturation_h2_full_vs_w17.csv\n"))

cat("\n=== Summary ===\n")
cat(sprintf("Mean |Δh²|:                %.4f\n", mean(results$delta_abs)))
cat(sprintf("Max |Δh²|:                 %.4f (%s)\n",
            max(results$delta_abs),
            results$trait[which.max(results$delta_abs)]))
cat(sprintf("Traits with |Δh²| ≥ 0.05:  %d / %d\n",
            sum(results$delta_abs >= 0.05), nrow(results)))
cat(sprintf("Rank-order Spearman ρ:     %.3f\n",
            cor(results$h2_full, results$h2_w17, method = "spearman")))

# Cross-check: full-window h² from this script should match canonical Table S2
cat("\n=== Cross-check: full-window h² (this script) vs canonical Table S2 ===\n")
canonical <- read.csv("results/metrics/maturation_h2_animal_model.csv")
for (tr in traits) {
  canon_h2 <- canonical$h2[canonical$metric == tr]
  if (length(canon_h2) == 1 && tr %in% results$trait) {
    my_h2 <- results$h2_full[results$trait == tr]
    cat(sprintf("  %-25s canonical=%.3f  this=%.3f  Δ=%+.3f\n",
                tr, canon_h2, my_h2, my_h2 - canon_h2))
  }
}

# PART 4: Pre-weaning growth as a fixed covariate (Methods S.10)
setwd(here::here())
suppressPackageStartupMessages({library(asreml); library(dplyr)})

Ainv    <- readRDS("data/processed/A_inverse_3gen.rds")
ped_ids <- attr(Ainv, "rowNames")
gmet <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv", stringsAsFactors = FALSE)
gthr <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv", stringsAsFactors = FALSE)
gped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
grow <- read.csv("data/metadata/growth_phenotypes.csv", stringsAsFactors = FALSE)

gdat <- gmet %>%
  inner_join(gped[, c("animal_id", "litter_id")], by = "animal_id") %>%
  left_join(gthr[, c("animal_id", "t70", "t75")], by = "animal_id") %>%
  inner_join(grow[, c("animal_id", "birthweight", "weaning_weight")], by = "animal_id") %>%
  filter(!is.na(birthweight) & !is.na(weaning_weight))
gdat$bw_c   <- as.numeric(scale(gdat$birthweight,    scale = FALSE))
gdat$wnwt_c <- as.numeric(scale(gdat$weaning_weight, scale = FALSE))
cat(sprintf("\n=== Growth-adjusted refit (n = %d animals with growth data) ===\n", nrow(gdat)))

gfit <- function(d, fx) {
  tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed = fx, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = d,
           na.action = na.method(y = "include"), trace = FALSE, maxit = 200)
  )), error = function(e) NULL)
}
gh2 <- function(fit) {
  vc <- summary(fit)$varcomp; rn <- rownames(vc)
  a <- max(vc[grep("vm.animal_id", rn)[1], "component"], 0)
  c <- max(vc[which(rn == "litter_id"), "component"], 0)
  e <- vc[grep("units", rn)[1], "component"]
  a / (a + c + e)
}
growth_rows <- list()
for (tr in c("maz_mean", "aitchison_to_adult_slope", "volatility_mean", "maz_slope",
             "bc_to_adult_slope", "intra_bc_var", "shannon_slope", "BC_FPC2", "t70", "t75")) {
  if (!tr %in% names(gdat)) next
  d <- gdat[, c("animal_id", "litter_id", "bw_c", "wnwt_c")]; d$y <- gdat[[tr]]
  d$animal_id <- factor(as.character(d$animal_id), levels = ped_ids)
  d$litter_id <- factor(d$litter_id)
  m0 <- gfit(d, y ~ 1); m1 <- gfit(d, y ~ 1 + bw_c + wnwt_c)
  if (is.null(m0) || is.null(m1)) next
  w  <- tryCatch(wald(m1), error = function(e) NULL)
  pc <- if (!is.null(w)) {
    cn <- intersect(c("Pr(Chisq)", "Pr(>F)"), colnames(w))[1]
    c(w["bw_c", cn], w["wnwt_c", cn])
  } else c(NA, NA)
  growth_rows[[tr]] <- data.frame(trait = tr, n = nrow(d),
    h2_unadj = gh2(m0), h2_adj = gh2(m1), delta_h2 = gh2(m1) - gh2(m0),
    p_birthweight = pc[1], p_weaning_weight = pc[2], stringsAsFactors = FALSE)
  cat(sprintf("  %-26s h²=%.3f -> %.3f  p(bw)=%.3f  p(wn)=%.3f\n",
              tr, gh2(m0), gh2(m1), pc[1], pc[2]))
}
growth_adj <- do.call(rbind, growth_rows)
write.csv(growth_adj, "results/metrics/maturation_h2_growth_adjusted.csv", row.names = FALSE)
cat(sprintf("Wrote %s\n", "results/metrics/maturation_h2_growth_adjusted.csv"))
sidak3 <- growth_adj[growth_adj$trait %in%
  c("maz_mean", "aitchison_to_adult_slope", "volatility_mean"), ]
cat(sprintf("Smallest covariate p among the three Sidak descriptors: %.4f\n",
            min(c(sidak3$p_birthweight, sidak3$p_weaning_weight), na.rm = TRUE)))
