#!/usr/bin/env Rscript
#________________________________________
# - File Name : 05_heritability_panel.R
# - Purpose : Univariate animal-model h2 panel, permutation tests, W6-W24 refit
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
suppressPackageStartupMessages({library(asreml); library(dplyr)})

Ainv_3g  <- readRDS("data/processed/A_inverse_3gen.rds")
ped_3g   <- attr(Ainv_3g, "rowNames")

Ainv_2g  <- readRDS("data/processed/A_inverse_proper.rds")
ped_2g   <- attr(Ainv_2g, "rowNames")

metrics  <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                     stringsAsFactors = FALSE)
thresh   <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv",
                     stringsAsFactors = FALSE)
ped_orig <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

dat <- metrics %>%
  inner_join(ped_orig[, c("animal_id", "litter_id", "line")], by = "animal_id") %>%
  left_join(thresh[, c("animal_id", "t70", "t75")], by = "animal_id")
cat(sprintf("Primary panel: %d animals\n", nrow(dat)))

fit_one <- function(y_vec, df, ped_ids, Ainv) {
  if (sum(!is.na(y_vec)) < 30 || sd(y_vec, na.rm = TRUE) < 1e-9) return(NULL)
  df$y <- y_vec
  df$animal_id <- factor(as.character(df$animal_id), levels = ped_ids)
  df$litter_id <- factor(df$litter_id)
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed = y ~ 1,
           random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = df,
           na.action = na.method(y = "include"),
           trace = FALSE, maxit = 200)
  )), error = function(e) NULL)
  if (is.null(fit) || !fit$converge) return(NULL)
  vc <- summary(fit)$varcomp; rn <- rownames(vc)
  s2_a <- max(vc[grep("vm.animal_id", rn)[1], "component"], 0)
  s2_c <- max(vc[which(rn == "litter_id"), "component"], 0)
  s2_e <- vc[grep("units", rn)[1], "component"]
  tot <- s2_a + s2_c + s2_e
  if (tot <= 0) return(NULL)
  h2_vp <- tryCatch(vpredict(fit, h2 ~ V2/(V1+V2+V3)), error = function(e) NULL)
  c2_vp <- tryCatch(vpredict(fit, c2 ~ V1/(V1+V2+V3)), error = function(e) NULL)
  list(n = sum(!is.na(y_vec)),
       h2 = s2_a/tot, h2_se = if (!is.null(h2_vp)) h2_vp$SE[1] else NA,
       c2 = s2_c/tot, c2_se = if (!is.null(c2_vp)) c2_vp$SE[1] else NA,
       s2_a = s2_a, s2_c = s2_c, s2_e = s2_e)
}

traits <- c("maz_mean", "maz_slope",
            "bc_to_adult_mean", "bc_to_adult_slope", "aitchison_to_adult_slope",
            "rel_maturity_w24",
            "volatility_mean", "volatility_slope", "intra_bc_var", "shannon_cv",
            "shannon_mean", "shannon_w24", "shannon_slope",
            "BC_FPC1", "BC_FPC2", "BC_FPC3",
            "shan_FPC1", "shan_FPC2", "shan_FPC3",
            "t70", "t75")

cat("\n=== 21-trait h² panel: 3-gen (primary) vs 2-gen pedigree ===\n")
primary <- data.frame(); compare <- data.frame()
for (tr in traits) {
  if (!tr %in% names(dat)) {
    cat(sprintf("  %-25s SKIP (not in dat)\n", tr)); next
  }
  y_vec <- dat[[tr]]
  df <- dat[, c("animal_id", "litter_id")]
  r3 <- fit_one(y_vec, df, ped_3g, Ainv_3g)
  r2 <- fit_one(y_vec, df, ped_2g, Ainv_2g)
  if (is.null(r3)) { cat(sprintf("  %-25s FIT FAILED (3g)\n", tr)); next }
  delta <- r3$h2 - (if (is.null(r2)) NA else r2$h2)
  cat(sprintf("  %-25s n=%d  h²_3g=%.3f (SE %.3f)  c²_3g=%.3f  |  Δ vs 2g = %+.3f\n",
              tr, r3$n, r3$h2, r3$h2_se, r3$c2, delta))
  primary <- rbind(primary, data.frame(
    trait = tr, n = r3$n,
    h2 = r3$h2, h2_se = r3$h2_se, c2 = r3$c2, c2_se = r3$c2_se,
    sigma2_A = r3$s2_a, sigma2_C = r3$s2_c, sigma2_E = r3$s2_e
  ))
  compare <- rbind(compare, data.frame(
    trait = tr,
    h2_3gen = r3$h2, h2_se_3gen = r3$h2_se,
    h2_2gen = if (is.null(r2)) NA else r2$h2,
    h2_se_2gen = if (is.null(r2)) NA else r2$h2_se,
    delta_h2 = delta, delta_abs = abs(delta)
  ))
}

write.csv(primary, "results/metrics/maturation_h2_primary_3gen.csv", row.names = FALSE)
write.csv(compare, "results/metrics/maturation_h2_3gen_vs_2gen.csv", row.names = FALSE)
cat("\nSaved:\n  results/metrics/maturation_h2_primary_3gen.csv  (new Table S1/S2 source)\n",
    "  results/metrics/maturation_h2_3gen_vs_2gen.csv     (3-vs-2-gen comparison)\n")

cat(sprintf("\n=== Summary (3-gen vs 2-gen, %d traits with valid fits) ===\n",
            sum(!is.na(compare$delta_h2))))
ok <- !is.na(compare$delta_abs)
cat(sprintf("Mean |Δh²|: %.4f\n", mean(compare$delta_abs[ok])))
cat(sprintf("Max |Δh²|:  %.4f (%s)\n",
            max(compare$delta_abs[ok]),
            compare$trait[ok][which.max(compare$delta_abs[ok])]))
cat(sprintf("Traits with |Δh²| ≥ 0.02: %d / %d\n",
            sum(compare$delta_abs[ok] >= 0.02), sum(ok)))
cat(sprintf("Traits with |Δh²| ≥ 0.05: %d / %d\n",
            sum(compare$delta_abs[ok] >= 0.05), sum(ok)))

# PART 2: Permutation test (full W3-W24 window)
setwd(here::here())
suppressPackageStartupMessages(library(asreml))

ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
ped_orig <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
all_ped_ids <- ped3$id

obs <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                stringsAsFactors = FALSE)
imp <- read.csv("results/metrics/maturation_all_metrics_per_animal_imputed.csv",
                stringsAsFactors = FALSE)
thr <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv",
                stringsAsFactors = FALSE)
src_map <- read.csv("paper/paper2_metric_source.csv", stringsAsFactors = FALSE)

# --- Full h^2 panel driven by metric-source map -----------------------------
# Resolve each metric's source and the column name in that source CSV.
# Convention: "either" -> "obs"; "drop" -> skipped.
# Special cases: aitchison_slope -> "aitchison_to_adult_slope"; t70_imp/t75_imp -> "t70"/"t75" in thr.

resolve_source <- function(s) if (s == "either") "obs" else s
resolve_col <- function(metric) {
  if (metric == "aitchison_slope") return("aitchison_to_adult_slope")
  if (metric == "t70_imp")          return("t70")
  if (metric == "t75_imp")          return("t75")
  metric
}

get_src <- function(name) {
  if (name == "obs") obs else if (name == "imp") imp else if (name == "thr") thr else stop("unknown source: ", name)
}

traits <- list()
for (i in seq_len(nrow(src_map))) {
  s <- src_map$source[i]
  if (s == "drop") next
  if (s %in% c("t70_imp", "t75_imp")) {
    # not used in this CSV column; the metric name itself encodes the file
  }
  src_label <- if (src_map$metric[i] %in% c("t70_imp", "t75_imp")) "thr" else resolve_source(s)
  traits[[length(traits) + 1]] <- list(
    name   = src_map$metric[i],
    source = src_label,
    col    = resolve_col(src_map$metric[i]))
}
cat(sprintf("Panel: %d non-drop metrics\n", length(traits)))

B <- 1999       # permutations per trait, per test (h^2 / c^2)
seed <- 20260518
set.seed(seed)

fit_vc <- function(df_fit) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed    = y ~ 1,
           random   = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units),
           data     = df_fit, na.action = na.method(y = "include"),
           trace = FALSE, maxit = 100))),
    error = function(e) NULL)
  if (is.null(fit)) return(c(h2 = NA_real_, c2 = NA_real_))
  vp_h2 <- tryCatch(vpredict(fit, h2 ~ V2 / (V1 + V2 + V3)), error = function(e) NULL)
  vp_c2 <- tryCatch(vpredict(fit, c2 ~ V1 / (V1 + V2 + V3)), error = function(e) NULL)
  c(h2 = if (is.null(vp_h2)) NA_real_ else vp_h2$Estimate[1],
    c2 = if (is.null(vp_c2)) NA_real_ else vp_c2$Estimate[1])
}

summary_rows <- data.frame()
null_dist_rows <- data.frame()

for (tr in traits) {
  cat(sprintf("\n=== %s (source = %s) ===\n", tr$name, tr$source))
  src <- get_src(tr$source)
  if (!(tr$col %in% names(src))) {
    cat(sprintf("  SKIP: column %s not found in %s\n", tr$col, tr$source))
    next
  }
  y_raw <- src[[tr$col]]
  animal_id <- src$animal_id
  ok <- !is.na(y_raw)
  if (sum(ok) < 30) {
    cat(sprintf("  SKIP: n_ok=%d\n", sum(ok))); next
  }
  y_raw <- y_raw[ok]
  animal_id <- animal_id[ok]

  lid <- ped_orig$litter_id[match(animal_id, ped_orig$animal_id)]

  df0 <- data.frame(
    y         = as.numeric(scale(y_raw)),
    animal_id = factor(animal_id, levels = all_ped_ids),
    litter_id = factor(lid))

  vc_obs <- fit_vc(df0)
  h2_obs <- vc_obs["h2"]; c2_obs <- vc_obs["c2"]
  cat(sprintf("  observed h^2 = %.3f, c^2 = %.3f  (n = %d)\n",
              h2_obs, c2_obs, nrow(df0)))

  if (is.na(h2_obs) && is.na(c2_obs)) next

  # --- Permutation A: animal_id shuffled (tests sigma^2_A > 0) -------------
  t0 <- proc.time()
  h2_perm_A <- numeric(B); c2_perm_A <- numeric(B)
  for (b in seq_len(B)) {
    df_p <- df0
    df_p$animal_id <- factor(sample(df0$animal_id), levels = all_ped_ids)
    vc_p <- fit_vc(df_p)
    h2_perm_A[b] <- vc_p["h2"]; c2_perm_A[b] <- vc_p["c2"]
    if (b %% 200 == 0) {
      el <- (proc.time() - t0)[3]
      cat(sprintf("    permA %d/%d (%.1f s, %.2f s/fit)\n", b, B, el, el / b))
    }
  }

  # --- Permutation C: litter_id shuffled (tests sigma^2_c > 0) -------------
  t0 <- proc.time()
  c2_perm_C <- numeric(B)
  for (b in seq_len(B)) {
    df_p <- df0
    df_p$litter_id <- factor(sample(df0$litter_id))
    vc_p <- fit_vc(df_p)
    c2_perm_C[b] <- vc_p["c2"]
    if (b %% 200 == 0) {
      el <- (proc.time() - t0)[3]
      cat(sprintf("    permC %d/%d (%.1f s, %.2f s/fit)\n", b, B, el, el / b))
    }
  }

  n_valid_A <- sum(!is.na(h2_perm_A))
  n_valid_C <- sum(!is.na(c2_perm_C))
  p_h2 <- (1 + sum(h2_perm_A >= h2_obs, na.rm = TRUE)) / (n_valid_A + 1)
  p_c2 <- (1 + sum(c2_perm_C >= c2_obs, na.rm = TRUE)) / (n_valid_C + 1)
  cat(sprintf("  h^2: null mean=%.3f sd=%.3f max=%.3f  p=%.4f  (n_valid=%d)\n",
              mean(h2_perm_A, na.rm = TRUE), sd(h2_perm_A, na.rm = TRUE),
              max(h2_perm_A, na.rm = TRUE), p_h2, n_valid_A))
  cat(sprintf("  c^2: null mean=%.3f sd=%.3f max=%.3f  p=%.4f  (n_valid=%d)\n",
              mean(c2_perm_C, na.rm = TRUE), sd(c2_perm_C, na.rm = TRUE),
              max(c2_perm_C, na.rm = TRUE), p_c2, n_valid_C))

  summary_rows <- rbind(summary_rows, data.frame(
    trait              = tr$name,
    source             = tr$source,
    n                  = nrow(df0),
    B                  = B,
    h2_obs             = h2_obs,
    h2_null_mean       = mean(h2_perm_A, na.rm = TRUE),
    h2_null_sd         = sd(h2_perm_A, na.rm = TRUE),
    h2_null_max        = max(h2_perm_A, na.rm = TRUE),
    h2_n_valid_perm    = n_valid_A,
    p_h2_one_tail      = p_h2,
    c2_obs             = c2_obs,
    c2_null_mean       = mean(c2_perm_C, na.rm = TRUE),
    c2_null_sd         = sd(c2_perm_C, na.rm = TRUE),
    c2_null_max        = max(c2_perm_C, na.rm = TRUE),
    c2_n_valid_perm    = n_valid_C,
    p_c2_one_tail      = p_c2,
    stringsAsFactors = FALSE))

  null_dist_rows <- rbind(null_dist_rows, data.frame(
    trait      = tr$name,
    perm_index = seq_len(B),
    h2_perm_A  = h2_perm_A,
    c2_perm_C  = c2_perm_C,
    stringsAsFactors = FALSE))
}

write.csv(summary_rows,
          "results/metrics/maturation_h2_permutation_test.csv",
          row.names = FALSE)
write.csv(null_dist_rows,
          "results/metrics/maturation_h2_permutation_null_dist.csv",
          row.names = FALSE)

cat("\n=== Summary ===\n")
print(summary_rows, row.names = FALSE, digits = 3)
cat("\nSaved: results/metrics/maturation_h2_permutation_test.csv\n")
cat("Saved: results/metrics/maturation_h2_permutation_null_dist.csv\n")

# PART 3: W3-balanced sensitivity audit (W6-W24 obs scalars + refits)
setwd(here::here())
suppressPackageStartupMessages({
  library(arrow)
  library(asreml)
})

ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
line_lookup <- setNames(ped$line, ped$animal_id)

# --- Load wide MAG matrix and curves ----------------------------------------

wide <- read_parquet("data/processed/mag_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
cols <- setdiff(names(wide), "animal_id")
mat <- as.matrix(wide[, cols]); rownames(mat) <- animal_ids
parts <- strsplit(cols, "_(?=W[0-9]+$)", perl = TRUE)
feat <- sapply(parts, `[`, 1)
tp   <- sapply(parts, `[`, 2)
features <- unique(feat)
# Chronological, not lexical: the consecutive-distance loop below steps through
# this order, so sort() would pair W24 with W6 and skip W9-W12.
tps_all  <- unique(tp)[order(as.integer(sub("W", "", unique(tp))))]
tps_w6   <- tps_all[as.integer(sub("W", "", tps_all)) >= 6]
week_map <- setNames(as.integer(sub("W", "", tps_all)), tps_all)

curves <- read.csv("results/metrics/maturation_curves_per_sample.csv",
                   stringsAsFactors = FALSE)
maz    <- read.csv("results/metrics/maz_per_sample.csv",
                   stringsAsFactors = FALSE)
maz$week <- as.integer(sub("^W", "", maz$tp))

# --- W24 centroid (across animals at W24) -----------------------------------

w24_cols <- paste0(features, "_W24")
w24_mat  <- mat[, w24_cols]
w24_centroid <- colMeans(w24_mat, na.rm = TRUE)
w24_centroid_prop <- w24_centroid / sum(w24_centroid)

clr <- function(x) {
  x <- pmax(x, 1e-9)
  log(x) - mean(log(x))
}
w24_clr <- clr(w24_centroid)

get_profile <- function(a, t) {
  unname(mat[a, paste0(features, "_", t)])
}

# --- Per-animal W6-W24 obs metrics ------------------------------------------

w6_results <- list()
for (a in animal_ids) {
  obs_tps <- character(0)
  for (t in tps_w6) {
    if (!is.na(mat[a, paste0(features[1], "_", t)])) obs_tps <- c(obs_tps, t)
  }
  if (length(obs_tps) < 3) next
  obs_weeks <- week_map[obs_tps]

  # bc / aitchison to W24 centroid
  bc_vals <- numeric(length(obs_tps))
  ait_vals <- numeric(length(obs_tps))
  profs <- list()
  for (i in seq_along(obs_tps)) {
    p <- get_profile(a, obs_tps[i])
    p_prop <- p / sum(p)
    profs[[obs_tps[i]]] <- p_prop
    bc_vals[i]  <- sum(abs(p_prop - w24_centroid_prop)) / 2
    ait_vals[i] <- sqrt(sum((clr(p) - w24_clr)^2))
  }

  # Stability (consecutive BC) and intra_bc_var (all-pairwise BC)
  consec_bc <- numeric(length(obs_tps) - 1)
  for (i in seq_len(length(obs_tps) - 1)) {
    consec_bc[i] <- sum(abs(profs[[obs_tps[i]]] - profs[[obs_tps[i + 1]]])) / 2
  }
  pairwise_bc <- c()
  for (i in seq_len(length(obs_tps) - 1)) {
    for (j in (i + 1):length(obs_tps)) {
      pairwise_bc <- c(pairwise_bc,
                       sum(abs(profs[[obs_tps[i]]] - profs[[obs_tps[j]]])) / 2)
    }
  }
  volatility_slope_w6 <- if (length(consec_bc) >= 3) {
    mid_weeks <- (obs_weeks[-length(obs_weeks)] + obs_weeks[-1]) / 2
    coef(lm(consec_bc ~ mid_weeks))[2]
  } else NA

  # Shannon slope and CV on W6+
  cur <- curves[curves$animal_id == a & curves$tp %in% obs_tps, ]
  shannon_vals <- cur$shannon[match(obs_tps, cur$tp)]
  shannon_slope_w6 <- if (sum(!is.na(shannon_vals)) >= 3) {
    coef(lm(shannon_vals ~ obs_weeks))[2]
  } else NA
  shannon_cv_w6 <- if (sum(!is.na(shannon_vals)) >= 3 &&
                       mean(shannon_vals, na.rm = TRUE) != 0) {
    sd(shannon_vals, na.rm = TRUE) / mean(shannon_vals, na.rm = TRUE)
  } else NA

  # MAZ slope on W6+
  m <- maz[maz$animal_id == a & maz$tp %in% obs_tps, ]
  maz_slope_w6 <- if (nrow(m) >= 3) {
    coef(lm(m$maz ~ m$week))[2]
  } else NA

  w6_results[[a]] <- data.frame(
    animal_id = a, line = line_lookup[[a]], n_obs_w6 = length(obs_tps),
    bc_to_adult_slope_w6      = coef(lm(bc_vals ~ obs_weeks))[2],
    aitchison_to_adult_slope_w6 = coef(lm(ait_vals ~ obs_weeks))[2],
    volatility_mean_w6        = mean(consec_bc),
    volatility_slope_w6       = volatility_slope_w6,
    intra_bc_var_w6           = var(pairwise_bc),
    shannon_slope_w6          = shannon_slope_w6,
    shannon_cv_w6             = shannon_cv_w6,
    maz_slope_w6              = maz_slope_w6,
    stringsAsFactors = FALSE
  )
}
w6_df <- do.call(rbind, w6_results)
write.csv(w6_df, "results/metrics/maturation_w6_obs_metrics.csv", row.names = FALSE)
cat(sprintf("Wrote %s (%d animals)\n",
            "results/metrics/maturation_w6_obs_metrics.csv", nrow(w6_df)))

# --- Line tests: full vs W6+ ------------------------------------------------

obs_full <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                     stringsAsFactors = FALSE)

pairs <- list(
  c("maz_slope",                "maz_slope_w6"),
  c("bc_to_adult_slope",        "bc_to_adult_slope_w6"),
  c("aitchison_to_adult_slope", "aitchison_to_adult_slope_w6"),
  c("volatility_mean",          "volatility_mean_w6"),
  c("volatility_slope",         "volatility_slope_w6"),
  c("intra_bc_var",             "intra_bc_var_w6"),
  c("shannon_slope",            "shannon_slope_w6"),
  c("shannon_cv",               "shannon_cv_w6")
)

merged <- merge(obs_full[, c("animal_id", "line",
                             unique(sapply(pairs, `[`, 1)))],
                w6_df[, c("animal_id",
                          unique(sapply(pairs, `[`, 2)))],
                by = "animal_id")

cat("\n=== SEN vs TOL: full window vs W6-W24 balanced window ===\n")
out_lines <- data.frame()
for (p in pairs) {
  full_col <- p[1]; w6_col <- p[2]
  tt_full <- t.test(merged[[full_col]] ~ merged$line)
  tt_w6   <- t.test(merged[[w6_col]]   ~ merged$line)
  cat(sprintf("  %-28s  full p=%.4f (SEN=%+.4f, TOL=%+.4f)  W6+ p=%.4f (SEN=%+.4f, TOL=%+.4f)\n",
              full_col,
              tt_full$p.value,
              tt_full$estimate[1], tt_full$estimate[2],
              tt_w6$p.value,
              tt_w6$estimate[1], tt_w6$estimate[2]))
  out_lines <- rbind(out_lines, data.frame(
    metric          = full_col,
    p_full          = tt_full$p.value,
    sen_mean_full   = tt_full$estimate[1],
    tol_mean_full   = tt_full$estimate[2],
    p_w6            = tt_w6$p.value,
    sen_mean_w6     = tt_w6$estimate[1],
    tol_mean_w6     = tt_w6$estimate[2],
    survives_w6     = tt_w6$p.value < 0.05,
    stringsAsFactors = FALSE
  ))
}
write.csv(out_lines,
          "results/metrics/maturation_w3_artifact_line_tests.csv",
          row.names = FALSE)
cat(sprintf("Wrote %s\n",
            "results/metrics/maturation_w3_artifact_line_tests.csv"))

# --- W6+ multivariate test: Hotelling's T^2 ---------------------------------
# Per-metric source rule: stability=obs(W6+), distance=imp (PACE-robust,
# unchanged from 11), threshold=imp, age/diversity=obs(W6+), FPCA=obs (FPCA
# fit on observed trajectories — weight on W3 sparse points is small under
# default fdapace bandwidth, treated as W3-robust here).

imp_full <- read.csv("results/metrics/maturation_all_metrics_per_animal_imputed.csv",
                     stringsAsFactors = FALSE)
thr      <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv",
                     stringsAsFactors = FALSE)

m1 <- merge(w6_df[, c("animal_id", "line", "maz_slope_w6", "intra_bc_var_w6", "shannon_slope_w6")],
            imp_full[, c("animal_id", "aitchison_to_adult_slope")],
            by = "animal_id")
m1 <- merge(m1, thr[, c("animal_id", "t75")], by = "animal_id")
fpc <- obs_full[, c("animal_id", "BC_FPC2")]
m1 <- merge(m1, fpc, by = "animal_id")
m1 <- m1[complete.cases(m1[, c("maz_slope_w6", "aitchison_to_adult_slope",
                                "t75", "intra_bc_var_w6",
                                "shannon_slope_w6", "BC_FPC2")]), ]

X <- as.matrix(scale(m1[, c("maz_slope_w6", "aitchison_to_adult_slope",
                             "t75", "intra_bc_var_w6",
                             "shannon_slope_w6", "BC_FPC2")]))
g <- factor(m1$line)

hotelling_T2 <- function(X, g) {
  lvl <- levels(g)
  X1 <- X[g == lvl[1], , drop = FALSE]
  X2 <- X[g == lvl[2], , drop = FALSE]
  n1 <- nrow(X1); n2 <- nrow(X2); p <- ncol(X)
  m1 <- colMeans(X1); m2 <- colMeans(X2)
  S1 <- cov(X1); S2 <- cov(X2)
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  delta <- m1 - m2
  T2 <- (n1 * n2 / (n1 + n2)) * t(delta) %*% MASS::ginv(Sp) %*% delta
  Fstat <- ((n1 + n2 - p - 1) / (p * (n1 + n2 - 2))) * T2
  pval <- pf(Fstat, p, n1 + n2 - p - 1, lower.tail = FALSE)
  list(T2 = as.numeric(T2), F = as.numeric(Fstat),
       df1 = p, df2 = n1 + n2 - p - 1, p = as.numeric(pval),
       n1 = n1, n2 = n2)
}
ht <- hotelling_T2(X, g)
cat(sprintf("\nHotelling T^2 on W6+ vector: T2=%.2f  F(%d,%d)=%.2f  p=%.4f  (n_SEN=%d, n_TOL=%d)\n",
            ht$T2, ht$df1, ht$df2, ht$F, ht$p, ht$n1, ht$n2))

# Mahalanobis permutation
S <- cov(X)
Sinv <- MASS::ginv(S)
mahal_obs <- {
  m_sen <- colMeans(X[g == "SEN", ]); m_tol <- colMeans(X[g == "TOL", ])
  d <- m_sen - m_tol
  as.numeric(t(d) %*% Sinv %*% d)
}
set.seed(42); B <- 9999
perm_d <- numeric(B)
for (b in seq_len(B)) {
  gp <- sample(g)
  d <- colMeans(X[gp == "SEN", ]) - colMeans(X[gp == "TOL", ])
  perm_d[b] <- as.numeric(t(d) %*% Sinv %*% d)
}
p_perm <- mean(perm_d >= mahal_obs)
cat(sprintf("Mahalanobis permutation (B=%d): observed=%.3f  p=%.4f\n",
            B, mahal_obs, p_perm))

mv_out <- data.frame(
  test = c("Hotelling_T2_W6", "Mahalanobis_perm_W6"),
  statistic = c(ht$T2, mahal_obs),
  df1 = c(ht$df1, NA), df2 = c(ht$df2, NA),
  Fstat = c(ht$F, NA),
  p_value = c(ht$p, p_perm),
  n_SEN = ht$n1, n_TOL = ht$n2,
  stringsAsFactors = FALSE
)
write.csv(mv_out, "results/metrics/maturation_w6_multivariate_tests.csv",
          row.names = FALSE)
cat(sprintf("Wrote %s\n",
            "results/metrics/maturation_w6_multivariate_tests.csv"))

cat("\n=== Summary ===\n")
cat("Metrics whose line test SURVIVES W6+ refit:\n")
print(out_lines[out_lines$survives_w6, c("metric", "p_full", "p_w6")])
cat("\nMetrics whose line test FAILS the W6+ refit (likely W3-driven):\n")
print(out_lines[!out_lines$survives_w6, c("metric", "p_full", "p_w6")])

# --- Alpha-diversity W6+ h² robustness check --------------------------------

a_full <- read.csv("results/metrics/alpha_diversity_extended_h2.csv",
                   stringsAsFactors = FALSE)
a_per_sample <- read.csv("results/metrics/alpha_diversity_extended_per_sample.csv",
                         stringsAsFactors = FALSE)
a_per_sample$week <- as.integer(sub("^W", "", a_per_sample$tp))
a_w6 <- a_per_sample[a_per_sample$week >= 6, ]

ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
litter_lookup <- setNames(as.character(ped$litter_id), ped$animal_id)

per_animal_summary <- function(df, col) {
  fit <- tryCatch(stats::lm(df[[col]] ~ df$week), error = function(e) NULL)
  list(
    mean_v = mean(df[[col]], na.rm = TRUE),
    slope  = if (!is.null(fit) && nrow(coef(summary(fit))) >= 2) coef(fit)[2] else NA,
    cv     = sd(df[[col]], na.rm = TRUE) / max(abs(mean(df[[col]], na.rm = TRUE)), 1e-9),
    w24    = if (any(df$tp == "W24")) df[[col]][df$tp == "W24"][1] else NA
  )
}

metrics_to_summarise <- c("shannon", "simpson", "invsimpson",
                          "richness", "pielou", "faith_pd")
metrics_to_summarise <- metrics_to_summarise[metrics_to_summarise %in% names(a_w6)]

per_animal_w6 <- data.frame(animal_id = unique(a_w6$animal_id))
for (m in metrics_to_summarise) {
  vals_mean <- vals_slope <- vals_cv <- vals_w24 <- numeric(nrow(per_animal_w6))
  for (i in seq_len(nrow(per_animal_w6))) {
    sub <- a_w6[a_w6$animal_id == per_animal_w6$animal_id[i], ]
    s <- per_animal_summary(sub, m)
    vals_mean[i]  <- s$mean_v
    vals_slope[i] <- s$slope
    vals_cv[i]    <- s$cv
    vals_w24[i]   <- s$w24
  }
  per_animal_w6[[paste0(m, "_mean")]]  <- vals_mean
  per_animal_w6[[paste0(m, "_slope")]] <- vals_slope
  per_animal_w6[[paste0(m, "_cv")]]    <- vals_cv
  per_animal_w6[[paste0(m, "_w24")]]   <- vals_w24
}
per_animal_w6$animal_id <- factor(per_animal_w6$animal_id, levels = ped3$id)
per_animal_w6$litter_id <- factor(as.character(litter_lookup[as.character(per_animal_w6$animal_id)]))

extract_pieces <- function(fit) {
  vc <- summary(fit)$varcomp; rn <- rownames(vc)
  s2_A <- {a <- grep("vm\\(animal_id", rn); if (length(a)) max(vc[a[1], "component"], 0) else 0}
  s2_c <- {c <- which(rn == "litter_id");   if (length(c)) max(vc[c, "component"], 0) else 0}
  s2_e <- {e <- grep("units", rn);          if (length(e)) vc[e[1], "component"] else NA}
  total <- s2_A + s2_c + s2_e
  h2vp <- tryCatch(vpredict(fit, h2 ~ V2 / (V1 + V2 + V3)), error = function(e) NULL)
  list(
    h2 = if (!is.null(h2vp)) h2vp$Estimate[1] else s2_A / total,
    h2_se = if (!is.null(h2vp)) h2vp$SE[1] else NA
  )
}

cat("\n=== Alpha-diversity W6+ h² robustness ===\n")
target_cols <- setdiff(names(per_animal_w6), c("animal_id", "litter_id"))
h2_w6 <- data.frame()
for (col in target_cols) {
  y_raw <- per_animal_w6[[col]]
  if (sum(!is.na(y_raw)) < 30) next
  if (sd(y_raw, na.rm = TRUE) < 1e-9) next
  df_fit <- data.frame(
    y = as.numeric(scale(y_raw)),
    animal_id = per_animal_w6$animal_id,
    litter_id = per_animal_w6$litter_id
  )
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = df_fit, trace = FALSE, maxit = 100,
           na.action = na.method(y = "include"))
  )), error = function(e) NULL)
  if (is.null(fit)) next
  p <- extract_pieces(fit)
  h2_full_row <- a_full[a_full$metric == col, ]
  h2_full_val <- if (nrow(h2_full_row) == 1) h2_full_row$h2 else NA
  h2_full_se  <- if (nrow(h2_full_row) == 1) h2_full_row$h2_se else NA
  h2_w6 <- rbind(h2_w6, data.frame(
    metric = col,
    h2_full = h2_full_val, h2_full_se = h2_full_se,
    h2_w6   = p$h2,        h2_w6_se   = p$h2_se,
    delta   = p$h2 - (if (is.na(h2_full_val)) NA else h2_full_val),
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  %-22s  full h²=%.3f±%.3f  W6+ h²=%.3f±%.3f  Δ=%+.3f\n",
              col, h2_full_val, h2_full_se, p$h2, p$h2_se,
              p$h2 - (if (is.na(h2_full_val)) NA else h2_full_val)))
}
write.csv(h2_w6, "results/metrics/alpha_diversity_w6_h2.csv", row.names = FALSE)
cat(sprintf("Wrote %s\n", "results/metrics/alpha_diversity_w6_h2.csv"))

# --- Per-animal scalar h² on W6+ (growth-trajectory scalars) ---------------
# Refits the headline per-animal h² estimates on the W6+ window to confirm
# whether the growth-trajectory heritabilities (maz_mean, intra_bc_var,
# volatility_mean, maz_slope) survive the W3 strip.

# Recompute maz_mean on W6+
maz_mean_w6 <- numeric(nrow(w6_df))
for (i in seq_len(nrow(w6_df))) {
  m <- maz[maz$animal_id == w6_df$animal_id[i] &
           maz$tp %in% tps_w6, ]
  maz_mean_w6[i] <- if (nrow(m) >= 3) mean(m$maz, na.rm = TRUE) else NA
}
w6_df$maz_mean_w6 <- maz_mean_w6

w6_for_h2 <- w6_df
w6_for_h2$animal_id <- factor(w6_for_h2$animal_id, levels = ped3$id)
w6_for_h2$litter_id <- factor(as.character(litter_lookup[as.character(w6_for_h2$animal_id)]))

cat("\n=== Per-animal scalar h² on W6+ window ===\n")
scalar_h2 <- data.frame()
scalars <- c("maz_mean_w6", "maz_slope_w6", "intra_bc_var_w6",
             "volatility_mean_w6", "shannon_slope_w6", "bc_to_adult_slope_w6",
             "aitchison_to_adult_slope_w6")
# Pull full-window references from obs_full
full_ref <- list(
  maz_mean_w6 = "maz_mean", maz_slope_w6 = "maz_slope",
  intra_bc_var_w6 = "intra_bc_var", volatility_mean_w6 = "volatility_mean",
  shannon_slope_w6 = "shannon_slope",
  bc_to_adult_slope_w6 = "bc_to_adult_slope",
  aitchison_to_adult_slope_w6 = "aitchison_to_adult_slope")

for (col in scalars) {
  y_raw <- w6_for_h2[[col]]
  if (sum(!is.na(y_raw)) < 30) next
  if (sd(y_raw, na.rm = TRUE) < 1e-9) next
  df_fit <- data.frame(
    y = as.numeric(scale(y_raw)),
    animal_id = w6_for_h2$animal_id,
    litter_id = w6_for_h2$litter_id)
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = df_fit, trace = FALSE, maxit = 100,
           na.action = na.method(y = "include"))
  )), error = function(e) NULL)
  if (is.null(fit)) next
  p <- extract_pieces(fit)
  scalar_h2 <- rbind(scalar_h2, data.frame(
    metric = full_ref[[col]],
    n = sum(!is.na(y_raw)),
    h2_w6 = p$h2, h2_w6_se = p$h2_se,
    stringsAsFactors = FALSE))
  cat(sprintf("  %-25s  n=%d  W6+ h²=%.3f±%.3f\n",
              full_ref[[col]], sum(!is.na(y_raw)), p$h2, p$h2_se))
}
write.csv(scalar_h2, "results/metrics/maturation_h2_w6_scalars.csv",
          row.names = FALSE)
cat(sprintf("Wrote %s\n", "results/metrics/maturation_h2_w6_scalars.csv"))

# --- Community-level per-week h² on bc_to_adult and shannon ----------------
# Refits per-week intercept animal model on (1) full data, (2) W6+ window only.

cur <- read.csv("results/metrics/maturation_curves_per_sample.csv",
                stringsAsFactors = FALSE)
cur$week <- as.integer(sub("^W", "", cur$tp))
cur <- merge(cur, ped[, c("animal_id", "litter_id")], by = "animal_id", all.x = TRUE)

fit_per_week_h2 <- function(df, trait_col, label) {
  out <- data.frame()
  weeks <- sort(unique(df$week))
  for (w in weeks) {
    sub <- df[df$week == w, ]
    sub <- sub[!is.na(sub[[trait_col]]), ]
    if (nrow(sub) < 30) next
    if (sd(sub[[trait_col]]) < 1e-9) next
    sub$animal_id <- factor(sub$animal_id, levels = ped3$id)
    sub$litter_id <- factor(as.character(sub$litter_id))
    sub$y <- as.numeric(scale(sub[[trait_col]]))
    fit <- tryCatch(suppressWarnings(suppressMessages(
      asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
             residual = ~ idv(units), data = sub, trace = FALSE, maxit = 100,
             na.action = na.method(y = "include"))
    )), error = function(e) NULL)
    if (is.null(fit)) next
    p <- extract_pieces(fit)
    out <- rbind(out, data.frame(
      label = label, trait = trait_col, week = w, n = nrow(sub),
      h2 = p$h2, h2_se = p$h2_se, stringsAsFactors = FALSE))
  }
  out
}

cat("\n=== Per-week community h²(t) ===\n")
pw_full <- rbind(
  fit_per_week_h2(cur,                      "bc_to_adult", "full"),
  fit_per_week_h2(cur,                      "shannon",     "full")
)
pw_w6 <- rbind(
  fit_per_week_h2(cur[cur$week >= 6, ],     "bc_to_adult", "w6"),
  fit_per_week_h2(cur[cur$week >= 6, ],     "shannon",     "w6")
)
pw_all <- rbind(pw_full, pw_w6)
write.csv(pw_all, "results/metrics/maturation_per_week_h2.csv",
          row.names = FALSE)

for (tr in c("bc_to_adult", "shannon")) {
  cat(sprintf("\n  %s\n", tr))
  for (lab in c("full", "w6")) {
    sub <- pw_all[pw_all$trait == tr & pw_all$label == lab, ]
    cat(sprintf("    %-5s window:  ", lab))
    cat(paste(sprintf("W%02d=%.2f", sub$week, sub$h2), collapse = "  "))
    cat("\n")
  }
}
cat(sprintf("\nWrote %s\n", "results/metrics/maturation_per_week_h2.csv"))

# --- Alpha-diversity line tests (full + W6+) --------------------------------

a_per_animal_full <- read.csv("results/metrics/alpha_diversity_extended_per_animal.csv",
                              stringsAsFactors = FALSE)
a_per_animal_full <- merge(a_per_animal_full,
                           ped[, c("animal_id", "line")],
                           by = "animal_id")
a_per_animal_w6 <- per_animal_w6
a_per_animal_w6$line <- ped$line[match(a_per_animal_w6$animal_id, ped$animal_id)]

cat("\n=== Alpha-diversity line tests (SEN vs TOL) ===\n")
alpha_lines <- data.frame()
for (col in target_cols) {
  if (!(col %in% names(a_per_animal_full))) next
  yf <- a_per_animal_full[[col]]
  yw <- a_per_animal_w6[[col]]
  if (sum(!is.na(yf)) < 30 || sum(!is.na(yw)) < 30) next
  if (sd(yf, na.rm = TRUE) < 1e-9 || sd(yw, na.rm = TRUE) < 1e-9) next
  ttf <- tryCatch(t.test(yf ~ a_per_animal_full$line), error = function(e) NULL)
  ttw <- tryCatch(t.test(yw ~ a_per_animal_w6$line),  error = function(e) NULL)
  if (is.null(ttf) || is.null(ttw)) next
  alpha_lines <- rbind(alpha_lines, data.frame(
    metric = col,
    p_full = ttf$p.value, sen_full = ttf$estimate[1], tol_full = ttf$estimate[2],
    p_w6   = ttw$p.value, sen_w6   = ttw$estimate[1], tol_w6   = ttw$estimate[2],
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  %-22s  full p=%.4f (SEN=%+.3f, TOL=%+.3f)  W6+ p=%.4f (SEN=%+.3f, TOL=%+.3f)\n",
              col, ttf$p.value, ttf$estimate[1], ttf$estimate[2],
              ttw$p.value, ttw$estimate[1], ttw$estimate[2]))
}
write.csv(alpha_lines, "results/metrics/alpha_diversity_w3_artifact_line_tests.csv",
          row.names = FALSE)
cat(sprintf("Wrote %s\n", "results/metrics/alpha_diversity_w3_artifact_line_tests.csv"))

# PART 4: Permutation test (W6-W24 balanced refit)
setwd(here::here())
suppressPackageStartupMessages(library(asreml))

ped3     <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
ped_orig <- read.csv("data/metadata/pedigree.csv",      stringsAsFactors = FALSE)
Ainv     <- readRDS("data/processed/A_inverse_3gen.rds")
all_ped_ids <- ped3$id

w6_obs <- read.csv("results/metrics/maturation_w6_obs_metrics.csv",
                   stringsAsFactors = FALSE)
maz    <- read.csv("results/metrics/maz_per_sample.csv",
                   stringsAsFactors = FALSE)

# --- Reconstruct maz_mean_w6 (per script 13_w3_artifact_audit.R) ------------
tps_w6 <- c("W6", "W9", "W12", "W15", "W17", "W21", "W24")
maz_mean_w6 <- numeric(nrow(w6_obs))
for (i in seq_len(nrow(w6_obs))) {
  m <- maz[maz$animal_id == w6_obs$animal_id[i] & maz$tp %in% tps_w6, ]
  maz_mean_w6[i] <- if (nrow(m) >= 3) mean(m$maz, na.rm = TRUE) else NA
}
w6_obs$maz_mean_w6 <- maz_mean_w6

# --- Panel definition -------------------------------------------------------
traits <- list(
  list(name = "maz_mean_w6",                 col = "maz_mean_w6"),
  list(name = "maz_slope_w6",                col = "maz_slope_w6"),
  list(name = "intra_bc_var_w6",             col = "intra_bc_var_w6"),
  list(name = "volatility_mean_w6",          col = "volatility_mean_w6"),
  list(name = "volatility_slope_w6",         col = "volatility_slope_w6"),
  list(name = "shannon_slope_w6",            col = "shannon_slope_w6"),
  list(name = "shannon_cv_w6",               col = "shannon_cv_w6"),
  list(name = "bc_to_adult_slope_w6",        col = "bc_to_adult_slope_w6"),
  list(name = "aitchison_to_adult_slope_w6", col = "aitchison_to_adult_slope_w6"))

B    <- 1999
seed <- 20260518
set.seed(seed)

fit_vc <- function(df_fit) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed    = y ~ 1,
           random   = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units),
           data     = df_fit, na.action = na.method(y = "include"),
           trace = FALSE, maxit = 100))),
    error = function(e) NULL)
  if (is.null(fit)) return(c(h2 = NA_real_, c2 = NA_real_))
  vp_h2 <- tryCatch(vpredict(fit, h2 ~ V2 / (V1 + V2 + V3)), error = function(e) NULL)
  vp_c2 <- tryCatch(vpredict(fit, c2 ~ V1 / (V1 + V2 + V3)), error = function(e) NULL)
  c(h2 = if (is.null(vp_h2)) NA_real_ else vp_h2$Estimate[1],
    c2 = if (is.null(vp_c2)) NA_real_ else vp_c2$Estimate[1])
}

summary_rows <- data.frame()

for (tr in traits) {
  cat(sprintf("\n=== %s ===\n", tr$name))
  if (!(tr$col %in% names(w6_obs))) {
    cat(sprintf("  SKIP: column %s not found\n", tr$col)); next
  }
  y_raw     <- w6_obs[[tr$col]]
  animal_id <- w6_obs$animal_id
  ok        <- !is.na(y_raw)
  if (sum(ok) < 30) { cat(sprintf("  SKIP: n_ok=%d\n", sum(ok))); next }
  y_raw     <- y_raw[ok]
  animal_id <- animal_id[ok]

  lid <- ped_orig$litter_id[match(animal_id, ped_orig$animal_id)]
  df0 <- data.frame(
    y         = as.numeric(scale(y_raw)),
    animal_id = factor(animal_id, levels = all_ped_ids),
    litter_id = factor(lid))

  vc_obs <- fit_vc(df0)
  h2_obs <- vc_obs["h2"]; c2_obs <- vc_obs["c2"]
  cat(sprintf("  observed h^2 = %.3f, c^2 = %.3f  (n = %d)\n",
              h2_obs, c2_obs, nrow(df0)))

  if (is.na(h2_obs) && is.na(c2_obs)) next

  # --- Permutation A: animal_id shuffled (tests sigma^2_A > 0) ---
  t0 <- proc.time()
  h2_perm_A <- numeric(B)
  for (b in seq_len(B)) {
    df_p <- df0
    df_p$animal_id <- factor(sample(df0$animal_id), levels = all_ped_ids)
    h2_perm_A[b] <- fit_vc(df_p)["h2"]
    if (b %% 200 == 0) {
      el <- (proc.time() - t0)[3]
      cat(sprintf("    permA %d/%d (%.1f s, %.2f s/fit)\n", b, B, el, el / b))
    }
  }

  # --- Permutation C: litter_id shuffled (tests sigma^2_c > 0) ---
  t0 <- proc.time()
  c2_perm_C <- numeric(B)
  for (b in seq_len(B)) {
    df_p <- df0
    df_p$litter_id <- factor(sample(df0$litter_id))
    c2_perm_C[b] <- fit_vc(df_p)["c2"]
    if (b %% 200 == 0) {
      el <- (proc.time() - t0)[3]
      cat(sprintf("    permC %d/%d (%.1f s, %.2f s/fit)\n", b, B, el, el / b))
    }
  }

  n_valid_A <- sum(!is.na(h2_perm_A))
  n_valid_C <- sum(!is.na(c2_perm_C))
  p_h2 <- (1 + sum(h2_perm_A >= h2_obs, na.rm = TRUE)) / (n_valid_A + 1)
  p_c2 <- (1 + sum(c2_perm_C >= c2_obs, na.rm = TRUE)) / (n_valid_C + 1)
  cat(sprintf("  h^2: null mean=%.3f max=%.3f  p=%.4f  (n_valid=%d)\n",
              mean(h2_perm_A, na.rm = TRUE), max(h2_perm_A, na.rm = TRUE),
              p_h2, n_valid_A))
  cat(sprintf("  c^2: null mean=%.3f max=%.3f  p=%.4f  (n_valid=%d)\n",
              mean(c2_perm_C, na.rm = TRUE), max(c2_perm_C, na.rm = TRUE),
              p_c2, n_valid_C))

  summary_rows <- rbind(summary_rows, data.frame(
    trait              = tr$name,
    n                  = nrow(df0),
    B                  = B,
    h2_obs             = h2_obs,
    h2_null_mean       = mean(h2_perm_A, na.rm = TRUE),
    h2_null_max        = max(h2_perm_A, na.rm = TRUE),
    h2_n_valid_perm    = n_valid_A,
    p_h2_one_tail      = p_h2,
    c2_obs             = c2_obs,
    c2_null_mean       = mean(c2_perm_C, na.rm = TRUE),
    c2_null_max        = max(c2_perm_C, na.rm = TRUE),
    c2_n_valid_perm    = n_valid_C,
    p_c2_one_tail      = p_c2,
    stringsAsFactors = FALSE))
}

write.csv(summary_rows,
          "results/metrics/maturation_h2_permutation_test_w6.csv",
          row.names = FALSE)

cat("\n=== W6+ permutation summary ===\n")
print(summary_rows, row.names = FALSE, digits = 3)
cat("\nSaved: results/metrics/maturation_h2_permutation_test_w6.csv\n")


cat("\nDone.\n")
