#!/usr/bin/env Rscript
#________________________________________
# - File Name : 06b_rrm_permutation.R
# - Purpose : Permutation test of pedigree structure in the random regression; Table S5A and S5B
# - Created By : Christian Maltecca
# - e-mail: cmaltec@ncsu.edu
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
# Sun Nov 01 15:30:00 2026 ------------------------------

setwd(here::here())

suppressPackageStartupMessages({ library(asreml); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 999L
# Residual structure per metric, chosen by AIC (06_rrm_h2t.R): week-specific for
# BC-to-adult and Shannon, homogeneous for MAZ (standardized within week).
RES <- c(bc_to_adult = "het", shannon = "het", maz = "hom")
res_f <- function(tr) if (RES[[tr]] == "het") ~ dsum(~ units | week_f) else ~ idv(units)
set.seed(20261101)

IN  <- "results/revision/inputs"
OUT <- "results/revision/rrm_permutation"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

Ainv <- readRDS(file.path(IN, "Ainv_3gen.rds"))
ids  <- attr(Ainv, "rowNames")
long <- read.csv(file.path(IN, "per_sample_long.csv"), stringsAsFactors = FALSE)

WEEKS <- c(3, 6, 9, 12, 15, 17, 21, 24)
WMEAN <- mean(WEEKS)
long$week_f    <- factor(long$week)
long$line      <- factor(long$line)
long$litter_id <- factor(long$litter_id)
long$wc        <- long$week - WMEAN
animals <- sort(unique(long$animal_id))

fit_A <- function(df, tr) {
  df$animal_id <- factor(df$ped_id, levels = ids)
  df <- df[order(df$week_f, df$animal_id), ]
  tryCatch(suppressWarnings(suppressMessages(asreml(
    fixed = y ~ 1 + week_f + line,
    random = ~ str(~ vm(animal_id, Ainv) + wc:vm(animal_id, Ainv),
                   ~ us(2):vm(animal_id, Ainv)) + litter_id,
    residual = res_f(tr), data = df, trace = FALSE, maxit = 300,
    na.action = na.method(y = "include")))), error = function(e) NULL)
}
fit_P <- function(df, tr) {
  df$pe <- factor(df$animal_id)
  df <- df[order(df$week_f, df$pe), ]
  suppressWarnings(suppressMessages(asreml(
    fixed = y ~ 1 + week_f + line,
    random = ~ str(~ pe + wc:pe, ~ us(2):id(78)) + litter_id,
    residual = res_f(tr), data = df, trace = FALSE, maxit = 300,
    na.action = na.method(y = "include"))))
}
share_A <- function(fit, tr) {
  vc <- summary(fit)$varcomp; rn <- rownames(vc); cp <- vc$component
  k11 <- cp[grep("us\\(2\\)_1:1$", rn)]; k21 <- cp[grep("us\\(2\\)_2:1$", rn)]
  k22 <- cp[grep("us\\(2\\)_2:2$", rn)]; sc <- cp[rn == "litter_id"]
  sapply(WEEKS, function(w) {
    t <- w - WMEAN; sa <- k11 + 2 * t * k21 + t^2 * k22
    se <- if (RES[[tr]] == "het") cp[rn == sprintf("week_f_%d!R", w)] else cp[rn == "units!units"]
    sa / (sa + sc + se)
  })
}

traits <- c(bc_to_adult = "BC-to-adult", shannon = "Shannon", maz = "MAZ")
summ <- list(); nulls <- list()
for (tr in names(traits)) {
  if (tr == "maz") set.seed(20261101)   # MAZ: own random stream, as published
  df <- long[!is.na(long[[tr]]), ]
  df$y <- as.numeric(scale(df[[tr]]))
  df$ped_id <- df$animal_id
  fP <- fit_P(df, tr); stopifnot(fP$converge)
  fA <- fit_A(df, tr); stopifnot(!is.null(fA), fA$converge)
  obs_dLL <- fA$loglik - fP$loglik
  obs_sh  <- share_A(fA, tr)
  cat(sprintf("\n=== %s: observed dLL = %.4f (A %.4f, PE %.4f) ===\n",
              traits[tr], obs_dLL, fA$loglik, fP$loglik))
  t0 <- Sys.time()
  res <- matrix(NA_real_, B, 3 + length(WEEKS),
                dimnames = list(NULL, c("dLL", "converged", "dLL_final", paste0("W", WEEKS))))
  for (b in seq_len(B)) {
    perm <- setNames(sample(animals), animals)
    dfb <- df; dfb$ped_id <- perm[dfb$animal_id]
    fb <- fit_A(dfb, tr)
    if (!is.null(fb) && !isTRUE(fb$converge))
      fb <- tryCatch(suppressWarnings(suppressMessages(update(fb, maxit = 500))),
                     error = function(e) fb)
    ok <- !is.null(fb) && isTRUE(fb$converge)
    res[b, "converged"] <- ok
    # last-iteration log-likelihood, kept for every fit that returned
    if (!is.null(fb)) res[b, "dLL_final"] <- fb$loglik - fP$loglik
    if (ok) { res[b, "dLL"] <- fb$loglik - fP$loglik; res[b, -(1:3)] <- share_A(fb, tr) }
    if (b %% 100 == 0) cat(sprintf("  %d/%d  %.1f min\n", b, B,
                                   as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  ok <- res[, "converged"] == 1
  nd <- res[ok, "dLL"]
  p  <- (1 + sum(nd >= obs_dLL)) / (sum(ok) + 1)
  p_cons <- (1 + sum(nd >= obs_dLL) + sum(!ok)) / (B + 1)
  nf <- res[!is.na(res[, "dLL_final"]), "dLL_final"]
  p_final <- (1 + sum(nf >= obs_dLL)) / (length(nf) + 1)
  n_nc_above <- sum(res[!ok & !is.na(res[, "dLL_final"]), "dLL_final"] >= obs_dLL)
  q95 <- apply(res[ok, -(1:3), drop = FALSE], 2, quantile, 0.95)
  pw  <- sapply(seq_along(WEEKS), function(i)
    (1 + sum(res[ok, 3 + i] >= obs_sh[i])) / (sum(ok) + 1))
  summ[[tr]] <- data.frame(row.names = NULL, trait = traits[tr], residual = RES[[tr]], obs_dLL = obs_dLL,
                           B_converged = sum(ok), B_total = B, p_dLL = p,
                           p_dLL_conservative = p_cons, p_dLL_all_final = p_final,
                           B_returned = length(nf), nonconverged_above_obs = n_nc_above,
                           null_dLL_q95 = quantile(nd, 0.95))
  nulls[[tr]] <- data.frame(row.names = NULL, trait = traits[tr], week = WEEKS, obs_share_A = obs_sh,
                            null_q95 = q95, p_week = pw)
  write.csv(data.frame(trait = traits[tr], res, row.names = NULL), file.path(OUT, sprintf("null_%s.csv", tr)),
            row.names = FALSE)
  print(summ[[tr]], row.names = FALSE); print(nulls[[tr]], row.names = FALSE, digits = 3)
}
write.csv(do.call(rbind, summ),  file.path(OUT, "perm_dLL_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, nulls), file.path(OUT, "perm_week_share.csv"),  row.names = FALSE)

# ---- Table S5 ------------------------------------------------------------------
fits <- read.csv("results/revision/rrm/rrm_model_fit.csv", stringsAsFactors = FALSE)
byw  <- read.csv("results/revision/rrm/rrm_by_week.csv", stringsAsFactors = FALSE)
pw   <- read.csv("results/revision/rrm_permutation/perm_week_share.csv", stringsAsFactors = FALSE)
pd   <- read.csv("results/revision/rrm_permutation/perm_dLL_summary.csv", stringsAsFactors = FALSE)

# ---- A. model comparison ------------------------------------------------------
lab <- c(M1_Arr_hom = "Additive only, homogeneous residual",
         M2_Arr_het = "Additive only, week-specific residual",
         M2L2_Aleg2_het = "Additive only, week-specific residual, second order",
         M3_Prr_het = "Permanent environment only, week-specific residual",
         M3h_Prr_hom = "Permanent environment only, homogeneous residual",
         M4_ArrPrr_het = "Additive and permanent environment, week-specific residual",
         M4h_ArrPrr_hom = "Additive and permanent environment, homogeneous residual")
A <- fits[fits$model %in% names(lab), c("trait", "model", "converge", "logLik", "npar", "AIC")]
A$Model <- lab[A$model]
A$logLik <- ifelse(is.na(A$logLik), "—", sprintf("%.2f", A$logLik))
A$AIC    <- ifelse(is.na(A$AIC), "—", sprintf("%.1f", A$AIC))
A$npar   <- ifelse(is.na(A$npar), "—", as.character(A$npar))
A$Converged <- ifelse(is.na(A$converge), "not fitted", ifelse(A$converge, "yes", "no"))
A <- A[order(match(A$trait, c("BC-to-adult", "Shannon", "MAZ")), match(A$model, names(lab))),
       c("trait", "Model", "Converged", "logLik", "npar", "AIC")]
write.csv(A, "results/revision/rrm/table_S5A.csv", row.names = FALSE)

# ---- B. per-week proportions --------------------------------------------------
MOD <- list("BC-to-adult" = c(A = "M2_Arr_het", P = "M3_Prr_het"),
            "Shannon"     = c(A = "M2_Arr_het", P = "M3_Prr_het"),
            "MAZ"         = c(A = "M1_Arr_hom", P = "M3h_Prr_hom"))
B <- do.call(rbind, lapply(names(MOD), function(tr) {
  a <- byw[byw$trait == tr & byw$model == MOD[[tr]]["A"], ]
  p <- byw[byw$trait == tr & byw$model == MOD[[tr]]["P"], ]
  q <- pw[pw$trait == tr, ]
  data.frame(trait = tr, week = a$week,
             additive_only = sprintf("%.3f (%.3f)", a$h2, a$h2_se),
             null_q95 = sprintf("%.3f", q$null_q95[match(a$week, q$week)]),
             p_perm = sprintf("%.3f", q$p_week[match(a$week, q$week)]),
             pe_only = sprintf("%.3f (%.3f)", p$pe2[match(a$week, p$week)], p$pe2_se[match(a$week, p$week)]))
}))
write.csv(B, "results/revision/rrm/table_S5B.csv", row.names = FALSE)
write.csv(pd, "results/revision/rrm/table_S5_perm_summary.csv", row.names = FALSE)
print(A, row.names = FALSE); print(B, row.names = FALSE); print(pd, row.names = FALSE)
