#!/usr/bin/env Rscript
#________________________________________
# - File Name : 08_genetic_correlations.R
# - Purpose : Sire BLUPs, Calo sire-product correlations, Bayesian bivariate fits
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
# Mon May 25 10:11:00 2026 ------------------------------

setwd(here::here())
suppressPackageStartupMessages(library(asreml))

# --- Load pedigree ----------------------------------------------------------

ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
ped_orig <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
all_ped_ids <- ped3$id

# "Cohort sires" = the 7 sires whose progeny are the 78 phenotyped animals.
# (Pulling sires from the 3-gen pedigree would include grandsires; we want only
# the direct sire population so BLUP downstream matches the 2-gen baseline.)
sire_ids <- unique(ped_orig$sire_id[!is.na(ped_orig$sire_id)])
cat(sprintf("Pedigree: %d individuals, %d cohort sires\n", length(all_ped_ids), length(sire_ids)))

# --- Per-metric source assignment -------------------------------------------

src <- read.csv("paper/paper2_metric_source.csv", stringsAsFactors = FALSE)
src_lookup <- setNames(src$source, src$metric)
get_source <- function(metric) {
  s <- src_lookup[[metric]]
  if (is.null(s) || is.na(s) || s == "either") "obs" else s
}

# --- Metric definitions: name -> column in CSV ------------------------------

metric_defs <- list(
  list(name = "maz_mean",          col = "maz_mean",                src_csv = "metrics"),
  list(name = "maz_slope",         col = "maz_slope",               src_csv = "metrics"),
  list(name = "bc_to_adult_mean",  col = "bc_to_adult_mean",        src_csv = "metrics"),
  list(name = "bc_to_adult_slope", col = "bc_to_adult_slope",       src_csv = "metrics"),
  list(name = "aitchison_slope",   col = "aitchison_to_adult_slope",src_csv = "metrics"),
  list(name = "rel_maturity_w24",  col = "rel_maturity_w24",        src_csv = "metrics"),
  list(name = "gompertz_k",        col = "g3m_k",                   src_csv = "gc"),
  list(name = "gompertz_ti",       col = "g3m_ti",                  src_csv = "gc"),
  list(name = "gompertz_A",        col = "g3m_A",                   src_csv = "gc"),
  list(name = "logistic_k",        col = "l3m_k",                   src_csv = "gc"),
  list(name = "logistic_t0",       col = "l3m_t0",                  src_csv = "gc"),
  list(name = "gd_constr_A",       col = "gd_constr_A",             src_csv = "gc"),
  list(name = "gd_constr_k",       col = "gd_constr_k",             src_csv = "gc"),
  list(name = "gd_constr_ti",      col = "gd_constr_ti",            src_csv = "gc"),
  list(name = "volatility_mean",   col = "volatility_mean",         src_csv = "metrics"),
  list(name = "volatility_slope",  col = "volatility_slope",        src_csv = "metrics"),
  list(name = "intra_bc_var",      col = "intra_bc_var",            src_csv = "metrics"),
  list(name = "shannon_cv",        col = "shannon_cv",              src_csv = "metrics"),
  list(name = "shannon_mean",      col = "shannon_mean",            src_csv = "metrics"),
  list(name = "shannon_slope",     col = "shannon_slope",           src_csv = "metrics"),
  list(name = "BC_FPC1",           col = "BC_FPC1",                 src_csv = "metrics"),
  list(name = "BC_FPC2",           col = "BC_FPC2",                 src_csv = "metrics"))

# --- Load both observed and imputed CSVs ------------------------------------

m_obs <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",        stringsAsFactors = FALSE)
m_imp <- read.csv("results/metrics/maturation_all_metrics_per_animal_imputed.csv",stringsAsFactors = FALSE)
g_obs <- read.csv("results/metrics/maturation_growth_curve_fits.csv",             stringsAsFactors = FALSE)
g_imp <- read.csv("results/metrics/maturation_growth_curve_fits_imputed.csv",     stringsAsFactors = FALSE)

# --- Build consolidated per-animal data frame -------------------------------

animal_ids <- m_obs$animal_id
df <- data.frame(animal_id = animal_ids, stringsAsFactors = FALSE)
df$litter_id <- ped_orig$litter_id[match(df$animal_id, ped_orig$animal_id)]
df$line      <- ped_orig$line     [match(df$animal_id, ped_orig$animal_id)]

source_used <- character(0)
for (md in metric_defs) {
  s <- get_source(md$name)
  source_used[md$name] <- s
  if (md$src_csv == "metrics") {
    src_df <- if (s == "imp") m_imp else m_obs
  } else {
    src_df <- if (s == "imp") g_imp else g_obs
  }
  if (!md$col %in% names(src_df)) {
    cat(sprintf("  WARN: %s not in %s (%s)\n", md$col, md$src_csv, s))
    df[[md$name]] <- NA_real_
  } else {
    df[[md$name]] <- src_df[[md$col]][match(df$animal_id, src_df$animal_id)]
  }
}

cat(sprintf("\nSource per metric:\n"))
for (m in names(source_used)) cat(sprintf("  %-20s %s\n", m, source_used[m]))

# --- Set up factors + ASReml-friendly structure -----------------------------

df$animal_id <- factor(df$animal_id, levels = all_ped_ids)
df$litter_id <- factor(df$litter_id)
df$line      <- factor(df$line)

# --- Univariate animal model + sire BLUP extraction -------------------------

cat("\n=== Univariate models + sire BLUP / PEV extraction ===\n")
sire_blup_list <- list()
sire_acc_rows  <- data.frame()   # long: sire_id, metric, blup, se, pev, h2, sigma2_A, accuracy
trait_h2 <- list()

for (md in metric_defs) {
  y <- df[[md$name]]
  if (sum(!is.na(y)) < 30 || sd(y, na.rm = TRUE) < 1e-6) {
    cat(sprintf("  %-20s skip (insufficient data)\n", md$name)); next
  }
  d_fit <- df; d_fit$y <- y

  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = d_fit,
           na.action = na.method(y = "include"), trace = FALSE, maxit = 100)
  )), error = function(e) NULL)
  if (is.null(fit)) {
    cat(sprintf("  %-20s asreml failed\n", md$name)); next
  }

  # h^2
  vc <- summary(fit)$varcomp
  s2_a <- max(vc[grep("vm\\(animal_id", rownames(vc))[1], "component"], 0)
  s2_c <- max(vc[which(rownames(vc) == "litter_id"), "component"], 0)
  s2_e <- vc[grep("units", rownames(vc))[1], "component"]
  h2   <- s2_a / (s2_a + s2_c + s2_e)

  # Animal BLUPs + SE (from the MME coefficient-matrix inverse via ASReml).
  cs <- summary(fit, coef = TRUE)$coef.random
  rn <- rownames(cs)
  blup_rows <- grep("^vm\\(animal_id, Ainv\\)_", rn)
  blup_vals <- as.numeric(cs[blup_rows, "solution"])
  blup_se   <- as.numeric(cs[blup_rows, "std.error"])
  blup_ids  <- sub("^vm\\(animal_id, Ainv\\)_", "", rn[blup_rows])
  names(blup_vals) <- blup_ids
  names(blup_se)   <- blup_ids

  sire_blups <- blup_vals[as.character(sire_ids)]
  sire_se    <- blup_se[as.character(sire_ids)]
  sire_pev   <- sire_se^2
  # Accuracy r_i = sqrt(1 - PEV_i / sigma^2_A); guard sigma^2_A == 0.
  if (s2_a > 1e-8) {
    sire_acc <- sqrt(pmax(0, 1 - sire_pev / s2_a))
  } else {
    sire_acc <- rep(NA_real_, length(sire_blups))
  }
  sire_blup_list[[md$name]] <- sire_blups
  trait_h2[[md$name]] <- h2

  sire_acc_rows <- rbind(sire_acc_rows, data.frame(
    sire_id   = as.character(sire_ids),
    metric    = md$name,
    blup      = sire_blups,
    se_blup   = sire_se,
    pev       = sire_pev,
    h2        = h2,
    sigma2_A  = s2_a,
    accuracy  = sire_acc,
    stringsAsFactors = FALSE))

  cat(sprintf("  %-20s h2=%.3f  sire BLUP range [%.3f, %.3f]  acc mean=%.3f range [%.3f, %.3f]\n",
              md$name, h2,
              min(sire_blups, na.rm=TRUE), max(sire_blups, na.rm=TRUE),
              mean(sire_acc, na.rm = TRUE),
              suppressWarnings(min(sire_acc, na.rm = TRUE)),
              suppressWarnings(max(sire_acc, na.rm = TRUE))))
}

# --- Threshold descriptors (t55, t60, t70, t75) -----------------------------
# Fit the same univariate animal model on the PACE-imputed time-to-threshold
# scalars so Table S17 can use a PEV-based accuracy for them too.

thr_csv <- "results/metrics/maturation_threshold_metrics_imputed.csv"
if (file.exists(thr_csv)) {
  thr_df <- read.csv(thr_csv, stringsAsFactors = FALSE)
  thr_metrics <- c("t55", "t60", "t70", "t75")
  for (mname in thr_metrics) {
    if (!mname %in% names(thr_df)) next
    d_fit <- data.frame(animal_id = factor(thr_df$animal_id, levels = all_ped_ids),
                        litter_id = factor(thr_df$litter_id),
                        y         = thr_df[[mname]],
                        stringsAsFactors = FALSE)
    if (sum(!is.na(d_fit$y)) < 30 || sd(d_fit$y, na.rm = TRUE) < 1e-6) next
    fit <- tryCatch(suppressWarnings(suppressMessages(
      asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
             residual = ~ idv(units), data = d_fit,
             na.action = na.method(y = "include"), trace = FALSE, maxit = 100)
    )), error = function(e) NULL)
    if (is.null(fit)) { cat(sprintf("  %-20s asreml failed\n", mname)); next }

    vc <- summary(fit)$varcomp
    s2_a <- max(vc[grep("vm\\(animal_id", rownames(vc))[1], "component"], 0)
    s2_c <- max(vc[which(rownames(vc) == "litter_id"), "component"], 0)
    s2_e <- vc[grep("units", rownames(vc))[1], "component"]
    h2   <- s2_a / (s2_a + s2_c + s2_e)

    cs <- summary(fit, coef = TRUE)$coef.random
    rn <- rownames(cs)
    blup_rows <- grep("^vm\\(animal_id, Ainv\\)_", rn)
    blup_vals <- as.numeric(cs[blup_rows, "solution"])
    blup_se   <- as.numeric(cs[blup_rows, "std.error"])
    blup_ids  <- sub("^vm\\(animal_id, Ainv\\)_", "", rn[blup_rows])
    names(blup_vals) <- blup_ids
    names(blup_se)   <- blup_ids

    sire_blups <- blup_vals[as.character(sire_ids)]
    sire_se    <- blup_se[as.character(sire_ids)]
    sire_pev   <- sire_se^2
    sire_acc   <- if (s2_a > 1e-8) sqrt(pmax(0, 1 - sire_pev / s2_a)) else rep(NA_real_, length(sire_blups))

    sire_blup_list[[mname]] <- sire_blups
    trait_h2[[mname]]        <- h2
    sire_acc_rows <- rbind(sire_acc_rows, data.frame(
      sire_id   = as.character(sire_ids),
      metric    = mname,
      blup      = sire_blups,
      se_blup   = sire_se,
      pev       = sire_pev,
      h2        = h2,
      sigma2_A  = s2_a,
      accuracy  = sire_acc,
      stringsAsFactors = FALSE))
    cat(sprintf("  %-20s h2=%.3f  sire BLUP range [%.3f, %.3f]  acc mean=%.3f range [%.3f, %.3f]\n",
                mname, h2,
                min(sire_blups, na.rm=TRUE), max(sire_blups, na.rm=TRUE),
                mean(sire_acc, na.rm = TRUE),
                suppressWarnings(min(sire_acc, na.rm = TRUE)),
                suppressWarnings(max(sire_acc, na.rm = TRUE))))
  }
}

# Save sire BLUPs (wide, one column per trait, as before)
sire_blup_df <- data.frame(sire_id = as.character(sire_ids), stringsAsFactors = FALSE)
for (m in names(sire_blup_list)) sire_blup_df[[m]] <- sire_blup_list[[m]]
write.csv(sire_blup_df, "results/metrics/maturation_sire_blups.csv", row.names = FALSE)
cat(sprintf("\nSaved sire BLUPs (%d sires x %d traits)\n",
            nrow(sire_blup_df), length(sire_blup_list)))

# Save per-trait per-sire PEV-based accuracy (long format)
write.csv(sire_acc_rows,
          "results/metrics/maturation_sire_blup_accuracy.csv",
          row.names = FALSE)
cat(sprintf("Saved sire BLUP accuracy (%d rows, %d traits x %d sires)\n",
            nrow(sire_acc_rows),
            length(unique(sire_acc_rows$metric)),
            length(unique(sire_acc_rows$sire_id))))

# --- Bivariate model + rg extraction ----------------------------------------

fit_rg <- function(y1, y2, animal_ids_vec, litter_vec, Ainv,
                   trait_levels = c("y1","y2"), ped_levels) {
  # Long format with trait factor
  n <- length(animal_ids_vec)
  d_long <- rbind(
    data.frame(animal_id = animal_ids_vec, litter_id = litter_vec,
               trait = trait_levels[1], value = y1, stringsAsFactors = FALSE),
    data.frame(animal_id = animal_ids_vec, litter_id = litter_vec,
               trait = trait_levels[2], value = y2, stringsAsFactors = FALSE))
  d_long$animal_id <- factor(d_long$animal_id, levels = ped_levels)
  d_long$litter_id <- factor(d_long$litter_id)
  d_long$trait     <- factor(d_long$trait, levels = trait_levels)

  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(value ~ trait,
           random = ~ us(trait):vm(animal_id, Ainv) + diag(trait):litter_id,
           residual = ~ id(animal_id):us(trait),
           data = d_long,
           na.action = na.method(y = "include"),
           trace = FALSE, maxit = 200)
  )), error = function(e) NULL)
  if (is.null(fit) || !fit$converge) return(list(rg = NA, rg_se = NA, converged = FALSE))

  vc <- summary(fit)$varcomp
  rn <- rownames(vc)
  # Genetic variance components
  va1 <- vc[rn == sprintf("trait:vm(animal_id, Ainv)!trait_%s:%s", trait_levels[1], trait_levels[1]),"component"]
  va2 <- vc[rn == sprintf("trait:vm(animal_id, Ainv)!trait_%s:%s", trait_levels[2], trait_levels[2]),"component"]
  cov_g <- vc[rn == sprintf("trait:vm(animal_id, Ainv)!trait_%s:%s", trait_levels[2], trait_levels[1]),"component"]
  if (length(va1) == 0 || length(va2) == 0 || length(cov_g) == 0 ||
      va1 <= 0 || va2 <= 0) {
    return(list(rg = NA, rg_se = NA, converged = TRUE,
                note = "boundary / nonpositive var_g"))
  }
  rg <- cov_g / sqrt(va1 * va2)
  rg <- max(min(rg, 1), -1)

  rg_se <- tryCatch(
    vpredict(fit, rg ~ V3 / sqrt(V1 * V5))$SE[1],
    error = function(e) NA)

  list(rg = rg, rg_se = rg_se, converged = TRUE,
       va1 = va1, va2 = va2, cov_g = cov_g)
}

# --- Pair list --------------------------------------------------------------

pairs <- list(
  c("maz_slope",       "aitchison_slope"),
  c("maz_mean",        "volatility_mean"),
  c("maz_mean",        "intra_bc_var"),
  c("maz_slope",       "volatility_mean"),
  c("volatility_mean", "intra_bc_var"),
  c("maz_mean",        "BC_FPC2"),
  c("maz_slope",       "BC_FPC2"),
  c("bc_to_adult_slope","aitchison_slope"),
  c("gd_constr_ti",    "maz_slope"),
  c("gd_constr_A",     "maz_mean"))

cat("\n=== Genetic correlations (r_pheno, r_EBV, rg_bivariate) ===\n")
results <- data.frame()
for (p in pairs) {
  m1 <- p[1]; m2 <- p[2]
  y1 <- df[[m1]]; y2 <- df[[m2]]
  if (is.null(y1) || is.null(y2)) {
    cat(sprintf("  %s x %s: missing column\n", m1, m2)); next
  }

  ok <- !is.na(y1) & !is.na(y2)
  r_pheno <- if (sum(ok) >= 10) cor(y1[ok], y2[ok]) else NA
  r_pheno_p <- if (sum(ok) >= 10) cor.test(y1[ok], y2[ok])$p.value else NA

  b1 <- sire_blup_list[[m1]]; b2 <- sire_blup_list[[m2]]
  if (!is.null(b1) && !is.null(b2)) {
    r_EBV <- cor(b1, b2, use = "complete.obs")
  } else r_EBV <- NA

  biv <- fit_rg(y1, y2, df$animal_id, df$litter_id, Ainv,
                trait_levels = c("y1","y2"), ped_levels = all_ped_ids)

  results <- rbind(results, data.frame(
    pair = paste(m1, "x", m2),
    metric_1 = m1, source_1 = source_used[m1],
    metric_2 = m2, source_2 = source_used[m2],
    n_animals = sum(ok),
    r_pheno = r_pheno, r_pheno_p = r_pheno_p,
    r_EBV = r_EBV,
    rg = biv$rg, rg_se = biv$rg_se,
    rg_converged = biv$converged,
    stringsAsFactors = FALSE))

  cat(sprintf("  %-22s r_pheno=%+.3f  r_EBV=%+.3f  rg=%s\n",
              paste(m1, "x", m2),
              ifelse(is.na(r_pheno), NA, r_pheno),
              ifelse(is.na(r_EBV),   NA, r_EBV),
              ifelse(is.na(biv$rg),
                     ifelse(biv$converged, "boundary", "no-converge"),
                     sprintf("%+.3f (SE %.3f)", biv$rg,
                             ifelse(is.na(biv$rg_se), NA, biv$rg_se)))))
}

# --- Save -------------------------------------------------------------------

write.csv(results, "results/metrics/maturation_genetic_correlations.csv", row.names = FALSE)
cat(sprintf("\nSaved: maturation_genetic_correlations.csv (%d pairs)\n", nrow(results)))
cat("Done.\n")

# PART 2: Sire-level bootstrap CIs on the Calo r_EBV statistic
setwd(here::here())

blups <- read.csv("results/metrics/maturation_sire_blups.csv",
                  stringsAsFactors = FALSE)
pairs <- read.csv("results/metrics/maturation_genetic_correlations.csv",
                  stringsAsFactors = FALSE)

n_sires <- nrow(blups)
cat(sprintf("Sire BLUP matrix: %d sires x %d traits\n",
            n_sires, ncol(blups) - 1))
stopifnot(n_sires == 7)

set.seed(20260518)
B <- 9999

bca_ci <- function(theta_hat, theta_boot, theta_jack, alpha = 0.05) {
  # Bias-correction
  p0 <- mean(theta_boot < theta_hat, na.rm = TRUE)
  if (p0 <= 0 || p0 >= 1) return(c(NA_real_, NA_real_))
  z0 <- qnorm(p0)
  # Acceleration from jackknife
  jbar <- mean(theta_jack, na.rm = TRUE)
  num <- sum((jbar - theta_jack)^3, na.rm = TRUE)
  den <- 6 * (sum((jbar - theta_jack)^2, na.rm = TRUE))^1.5
  a <- if (den > 0) num / den else 0
  z_lo <- qnorm(alpha / 2)
  z_hi <- qnorm(1 - alpha / 2)
  alpha_lo <- pnorm(z0 + (z0 + z_lo) / (1 - a * (z0 + z_lo)))
  alpha_hi <- pnorm(z0 + (z0 + z_hi) / (1 - a * (z0 + z_hi)))
  quantile(theta_boot, c(alpha_lo, alpha_hi), na.rm = TRUE, names = FALSE)
}

out <- data.frame()
for (i in seq_len(nrow(pairs))) {
  m1 <- pairs$metric_1[i]
  m2 <- pairs$metric_2[i]
  if (!m1 %in% names(blups) || !m2 %in% names(blups)) next
  u1 <- blups[[m1]]
  u2 <- blups[[m2]]
  ok <- is.finite(u1) & is.finite(u2)
  if (sum(ok) < 4) next
  u1 <- u1[ok]
  u2 <- u2[ok]
  n_eff <- length(u1)

  r_hat <- cor(u1, u2)

  # Bootstrap
  r_boot <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n_eff, n_eff, replace = TRUE)
    if (length(unique(idx)) < 3) { r_boot[b] <- NA_real_; next }
    s1 <- u1[idx]; s2 <- u2[idx]
    if (sd(s1) == 0 || sd(s2) == 0) { r_boot[b] <- NA_real_; next }
    r_boot[b] <- cor(s1, s2)
  }
  n_valid <- sum(!is.na(r_boot))

  # Jackknife for BCa
  r_jack <- vapply(seq_len(n_eff), function(k) {
    s1 <- u1[-k]; s2 <- u2[-k]
    if (sd(s1) == 0 || sd(s2) == 0) NA_real_ else cor(s1, s2)
  }, numeric(1))

  pct_ci <- quantile(r_boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  bca <- bca_ci(r_hat, r_boot, r_jack)

  # Two-sided percentile-bootstrap p-value vs r = 0
  p_boot <- 2 * min(mean(r_boot <= 0, na.rm = TRUE),
                    mean(r_boot >= 0, na.rm = TRUE))

  out <- rbind(out, data.frame(
    pair        = pairs$pair[i],
    metric_1    = m1,
    metric_2    = m2,
    n_sires     = n_eff,
    r_EBV       = r_hat,
    pct_lo      = pct_ci[1],
    pct_hi      = pct_ci[2],
    bca_lo      = bca[1],
    bca_hi      = bca[2],
    p_boot_2s   = p_boot,
    n_valid_boot = n_valid,
    stringsAsFactors = FALSE))

  cat(sprintf("  %-45s r=%+0.3f  pct[%+0.2f, %+0.2f]  bca[%+0.2f, %+0.2f]  p=%.3f\n",
              pairs$pair[i], r_hat, pct_ci[1], pct_ci[2],
              bca[1], bca[2], p_boot))
}

write.csv(out,
          "results/metrics/maturation_calo_bootstrap_ci.csv",
          row.names = FALSE)
cat(sprintf("\nSaved: results/metrics/maturation_calo_bootstrap_ci.csv (%d pairs)\n",
            nrow(out)))

# PART 3: Bayesian bivariate animal model (MCMCglmm) for two illustrative pairs
setwd(here::here())
suppressPackageStartupMessages({
  library(MCMCglmm)
})

obs <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                stringsAsFactors = FALSE)
imp <- read.csv("results/metrics/maturation_all_metrics_per_animal_imputed.csv",
                stringsAsFactors = FALSE)
ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

# Pedigree in MCMCglmm format: id, sire (dam), dam (litter proxy)
ped_mcmc <- data.frame(
  id   = ped3$id,
  dam  = ped3$dam,
  sire = ped3$sire,
  stringsAsFactors = FALSE
)

# Per-metric source rule (per Supp §S2)
get_metric <- function(name, source) {
  src <- if (source == "obs") obs else imp
  if (name == "aitchison_slope") name <- "aitchison_to_adult_slope"
  src[, c("animal_id", name)]
}

pairs <- list(
  list(name = "BC_slope_x_Aitchison_slope",
       m1 = "bc_to_adult_slope", s1 = "obs",
       m2 = "aitchison_slope",   s2 = "imp"),
  list(name = "MAZ_slope_x_Aitchison_slope",
       m1 = "maz_slope",         s1 = "obs",
       m2 = "aitchison_slope",   s2 = "imp")
)

run_bivariate <- function(name, m1, s1, m2, s2,
                          nitt = 60000, burnin = 10000, thin = 50) {
  d1 <- get_metric(m1, s1); colnames(d1)[2] <- "y1"
  d2 <- get_metric(m2, s2); colnames(d2)[2] <- "y2"
  df <- merge(d1, d2, by = "animal_id")
  df <- merge(df, ped[, c("animal_id", "litter_id")], by = "animal_id", all.x = TRUE)
  df <- df[complete.cases(df[, c("y1", "y2")]), ]
  df$animal <- as.character(df$animal_id)
  df$litter_id <- as.factor(df$litter_id)
  df$y1 <- as.numeric(scale(df$y1))
  df$y2 <- as.numeric(scale(df$y2))

  cat(sprintf("\n=== Pair: %s ===\n", name))
  cat(sprintf("  n_animals = %d\n", nrow(df)))

  # Weakly-informative prior centered at zero
  prior_us <- list(
    G = list(
      G1 = list(V = diag(2) * 0.5, nu = 1.002),  # animal genetic
      G2 = list(V = diag(2) * 0.5, nu = 1.002)   # litter
    ),
    R = list(V = diag(2) * 0.5, nu = 1.002)
  )

  set.seed(42)
  fit <- tryCatch(
    MCMCglmm(
      cbind(y1, y2) ~ trait - 1,
      random = ~ us(trait):animal + us(trait):litter_id,
      rcov   = ~ us(trait):units,
      family = c("gaussian", "gaussian"),
      data = df, pedigree = ped_mcmc, prior = prior_us,
      nitt = nitt, burnin = burnin, thin = thin, verbose = FALSE
    ),
    error = function(e) { message("MCMCglmm error: ", conditionMessage(e)); NULL }
  )

  if (is.null(fit)) {
    return(data.frame(
      pair = name, n = nrow(df),
      var_y1_g = NA, var_y2_g = NA, cov_y1y2_g = NA,
      var_y1_l = NA, var_y2_l = NA, cov_y1y2_l = NA,
      var_y1_e = NA, var_y2_e = NA, cov_y1y2_e = NA,
      h2_y1_post_mean = NA, h2_y1_HPD_low = NA, h2_y1_HPD_high = NA,
      h2_y2_post_mean = NA, h2_y2_HPD_low = NA, h2_y2_HPD_high = NA,
      rg_post_mean = NA, rg_HPD_low = NA, rg_HPD_high = NA,
      converged = FALSE, stringsAsFactors = FALSE
    ))
  }

  # Posterior summaries
  vc <- fit$VCV
  # Animal genetic covariance: columns named "trait...:trait...animal"
  g_y1y1 <- vc[, "traity1:traity1.animal"]
  g_y2y2 <- vc[, "traity2:traity2.animal"]
  g_y1y2 <- vc[, "traity1:traity2.animal"]
  l_y1y1 <- vc[, "traity1:traity1.litter_id"]
  l_y2y2 <- vc[, "traity2:traity2.litter_id"]
  l_y1y2 <- vc[, "traity1:traity2.litter_id"]
  e_y1y1 <- vc[, "traity1:traity1.units"]
  e_y2y2 <- vc[, "traity2:traity2.units"]
  e_y1y2 <- vc[, "traity1:traity2.units"]

  rg_post  <- g_y1y2 / sqrt(g_y1y1 * g_y2y2)
  h2_y1_post <- g_y1y1 / (g_y1y1 + l_y1y1 + e_y1y1)
  h2_y2_post <- g_y2y2 / (g_y2y2 + l_y2y2 + e_y2y2)

  cat(sprintf("  rg posterior mean = %+.3f, 95%% HPD [%+.3f, %+.3f]\n",
              mean(rg_post), HPDinterval(rg_post)[1], HPDinterval(rg_post)[2]))
  cat(sprintf("  h²(%s) post mean = %.3f, HPD [%.3f, %.3f]\n",
              m1, mean(h2_y1_post), HPDinterval(h2_y1_post)[1], HPDinterval(h2_y1_post)[2]))
  cat(sprintf("  h²(%s) post mean = %.3f, HPD [%.3f, %.3f]\n",
              m2, mean(h2_y2_post), HPDinterval(h2_y2_post)[1], HPDinterval(h2_y2_post)[2]))

  data.frame(
    pair = name, n = nrow(df),
    var_y1_g = mean(g_y1y1), var_y2_g = mean(g_y2y2), cov_y1y2_g = mean(g_y1y2),
    var_y1_l = mean(l_y1y1), var_y2_l = mean(l_y2y2), cov_y1y2_l = mean(l_y1y2),
    var_y1_e = mean(e_y1y1), var_y2_e = mean(e_y2y2), cov_y1y2_e = mean(e_y1y2),
    h2_y1_post_mean = mean(h2_y1_post),
    h2_y1_HPD_low  = HPDinterval(h2_y1_post)[1],
    h2_y1_HPD_high = HPDinterval(h2_y1_post)[2],
    h2_y2_post_mean = mean(h2_y2_post),
    h2_y2_HPD_low  = HPDinterval(h2_y2_post)[1],
    h2_y2_HPD_high = HPDinterval(h2_y2_post)[2],
    rg_post_mean = mean(rg_post),
    rg_HPD_low   = HPDinterval(rg_post)[1],
    rg_HPD_high  = HPDinterval(rg_post)[2],
    converged = TRUE, stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(pairs, function(p) {
  run_bivariate(p$name, p$m1, p$s1, p$m2, p$s2)
}))

write.csv(results, "results/metrics/maturation_mcmcglmm_bivariate.csv",
          row.names = FALSE)
cat("\n=== output ===\nresults/metrics/maturation_mcmcglmm_bivariate.csv\n")
print(results)
