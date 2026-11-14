#!/usr/bin/env Rscript
#________________________________________
# - File Name : 02_per_animal_phenotypes.R
# - Purpose : Per-animal trajectory descriptors: 24 scalars across six metric classes
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
# Mon May 25 10:10:00 2026 ------------------------------

setwd(here::here())
suppressPackageStartupMessages({
  library(arrow)
  library(vegan)
})

# --- Load data ---------------------------------------------------------------

wide <- read_parquet("data/processed/mag_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
feat_cols <- setdiff(names(wide), "animal_id")
mat <- as.matrix(wide[, feat_cols])
rownames(mat) <- animal_ids

col_parts <- strsplit(feat_cols, "_(?=W[0-9]+$)", perl = TRUE)
feature_names <- unique(sapply(col_parts, function(x) x[1]))
tp_labels <- unique(sapply(col_parts, function(x) x[2]))
week_map <- c(W3 = 3, W6 = 6, W9 = 9, W12 = 12, W15 = 15,
              W17 = 17, W21 = 21, W24 = 24)

ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
line_lookup <- setNames(ped$line, ped$animal_id)

curves <- read.csv("results/metrics/maturation_curves_per_sample.csv",
                    stringsAsFactors = FALSE)
maz_animal <- read.csv("results/metrics/maz_per_animal.csv",
                        stringsAsFactors = FALSE)

cat("Loaded:", length(animal_ids), "animals,", length(feature_names), "features\n")

# --- Helper functions --------------------------------------------------------

get_profile <- function(animal, tp) {
  cols <- paste0(feature_names, "_", tp)
  v <- mat[animal, cols]
  if (all(is.na(v))) return(NULL)
  v
}

clr <- function(x) {
  x[x == 0] <- 0.5
  lx <- log(x)
  lx - mean(lx)
}

# --- W24 adult centroid (proportions) ----------------------------------------

w24_profiles <- list()
for (a in animal_ids) {
  p <- get_profile(a, "W24")
  if (!is.null(p)) w24_profiles[[a]] <- p
}
w24_mat <- do.call(rbind, w24_profiles)
w24_centroid <- colMeans(w24_mat)
w24_centroid_prop <- w24_centroid / sum(w24_centroid)

# --- Per-animal metrics ------------------------------------------------------

cat("\nComputing per-animal metrics...\n")
results <- list()

for (ai in seq_along(animal_ids)) {
  a <- animal_ids[ai]
  line <- line_lookup[[a]]

  # Identify observed timepoints
  obs_tps <- character(0)
  for (tp in tp_labels) {
    if (!is.na(mat[a, paste0(feature_names[1], "_", tp)])) {
      obs_tps <- c(obs_tps, tp)
    }
  }
  obs_weeks <- week_map[obs_tps]
  n_obs <- length(obs_tps)
  if (n_obs < 3) next

  # -- Class 1: Age prediction (MAZ) --
  maz_row <- maz_animal[maz_animal$animal_id == a, ]
  maz_mean <- if (nrow(maz_row) == 1) maz_row$mean_maz else NA
  maz_slope <- if (nrow(maz_row) == 1) maz_row$maz_slope else NA

  # -- Class 2: Distance to adult --
  bc_vals <- numeric(n_obs)
  aitchison_vals <- numeric(n_obs)
  for (i in seq_along(obs_tps)) {
    p <- get_profile(a, obs_tps[i])
    p_prop <- p / sum(p)
    bc_vals[i] <- sum(abs(p_prop - w24_centroid_prop)) / 2
    aitchison_vals[i] <- sqrt(sum((clr(p) - clr(w24_centroid))^2))
  }

  bc_to_adult_mean <- mean(bc_vals)
  bc_to_adult_slope <- coef(lm(bc_vals ~ obs_weeks))[2]
  bc_to_adult_w24 <- bc_vals[which(obs_tps == "W24")]
  if (length(bc_to_adult_w24) == 0) bc_to_adult_w24 <- NA

  aitchison_to_adult_mean <- mean(aitchison_vals)
  aitchison_to_adult_slope <- coef(lm(aitchison_vals ~ obs_weeks))[2]

  bc_first <- bc_vals[1]
  rel_maturity_w24 <- if (!is.na(bc_to_adult_w24) && bc_first > 0) {
    1 - bc_to_adult_w24 / bc_first
  } else NA

  # -- Class 3: Growth curve parameters (preliminary Gompertz/logistic) --
  # Fit to 1 - BC (maturity, increasing toward 1)
  maturity_vals <- 1 - bc_vals

  gompertz_k <- NA; gompertz_ti <- NA; gompertz_A <- NA; gompertz_converged <- FALSE
  tryCatch({
    fit <- nls(y ~ A * exp(-exp(-k * (t - ti))),
               data = data.frame(y = maturity_vals, t = obs_weeks),
               start = list(A = max(maturity_vals) + 0.05, k = 0.15, ti = 10),
               control = nls.control(maxiter = 200, warnOnly = TRUE))
    gompertz_k <- coef(fit)["k"]
    gompertz_ti <- coef(fit)["ti"]
    gompertz_A <- coef(fit)["A"]
    gompertz_converged <- fit$convInfo$isConv
  }, error = function(e) {})

  logistic_k <- NA; logistic_t0 <- NA; logistic_L <- NA
  tryCatch({
    fit_log <- nls(y ~ L / (1 + exp(-k * (t - t0))),
                   data = data.frame(y = maturity_vals, t = obs_weeks),
                   start = list(L = max(maturity_vals) + 0.05, k = 0.15, t0 = 10),
                   control = nls.control(maxiter = 200, warnOnly = TRUE))
    logistic_k <- coef(fit_log)["k"]
    logistic_t0 <- coef(fit_log)["t0"]
    logistic_L <- coef(fit_log)["L"]
  }, error = function(e) {})

  # -- Class 4: Stability / volatility --
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

  # Intra-individual BC variance (all pairwise)
  all_profiles <- list()
  for (tp in obs_tps) {
    p <- get_profile(a, tp)
    all_profiles[[tp]] <- p / sum(p)
  }
  pairwise_bc <- c()
  for (i in seq_len(n_obs - 1)) {
    for (j in (i + 1):n_obs) {
      pairwise_bc <- c(pairwise_bc,
        sum(abs(all_profiles[[obs_tps[i]]] - all_profiles[[obs_tps[j]]])) / 2)
    }
  }
  intra_bc_var <- var(pairwise_bc)

  animal_curves <- curves[curves$animal_id == a, ]
  shannon_vals <- animal_curves$shannon[match(obs_tps, animal_curves$tp)]
  shannon_cv <- sd(shannon_vals, na.rm = TRUE) / mean(shannon_vals, na.rm = TRUE)

  # -- Class 5: Diversity trajectory --
  shannon_w24 <- shannon_vals[which(obs_tps == "W24")]
  if (length(shannon_w24) == 0) shannon_w24 <- NA
  shannon_slope <- coef(lm(shannon_vals ~ obs_weeks))[2]
  shannon_mean <- mean(shannon_vals, na.rm = TRUE)

  richness_vals <- animal_curves$richness[match(obs_tps, animal_curves$tp)]
  richness_slope <- coef(lm(richness_vals ~ obs_weeks))[2]
  richness_w24 <- richness_vals[which(obs_tps == "W24")]
  if (length(richness_w24) == 0) richness_w24 <- NA

  results[[ai]] <- data.frame(
    animal_id = a, line = line, n_obs = n_obs,
    maz_mean = maz_mean, maz_slope = maz_slope,
    bc_to_adult_mean = bc_to_adult_mean,
    bc_to_adult_slope = bc_to_adult_slope,
    bc_to_adult_w24 = bc_to_adult_w24,
    aitchison_to_adult_mean = aitchison_to_adult_mean,
    aitchison_to_adult_slope = aitchison_to_adult_slope,
    rel_maturity_w24 = rel_maturity_w24,
    gompertz_k = gompertz_k, gompertz_ti = gompertz_ti, gompertz_A = gompertz_A,
    gompertz_converged = gompertz_converged,
    logistic_k = logistic_k, logistic_t0 = logistic_t0, logistic_L = logistic_L,
    volatility_mean = volatility_mean,
    volatility_slope = volatility_slope,
    intra_bc_var = intra_bc_var,
    shannon_cv = shannon_cv,
    shannon_mean = shannon_mean, shannon_w24 = shannon_w24,
    shannon_slope = shannon_slope,
    richness_w24 = richness_w24, richness_slope = richness_slope,
    stringsAsFactors = FALSE
  )

  if (ai %% 20 == 0) cat(sprintf("  %d/%d animals\n", ai, length(animal_ids)))
}

df <- do.call(rbind, results)
cat(sprintf("\nComputed metrics for %d animals\n", nrow(df)))

# --- Merge FPCA scores (Class 6) --------------------------------------------

bc_fpc <- read.csv("results/metrics/maturation_bc_fpc_scores.csv",
                    stringsAsFactors = FALSE)
shan_fpc <- read.csv("results/metrics/maturation_fpc_scores.csv",
                      stringsAsFactors = FALSE)

bc_fpc_cols <- grep("^BC_FPC", names(bc_fpc), value = TRUE)
df <- merge(df, bc_fpc[, c("animal_id", bc_fpc_cols)],
            by = "animal_id", all.x = TRUE)

shan_fpc_cols <- grep("^FPC", names(shan_fpc), value = TRUE)
shan_rename <- setNames(paste0("shan_", shan_fpc_cols), shan_fpc_cols)
for (old in names(shan_rename)) {
  names(shan_fpc)[names(shan_fpc) == old] <- shan_rename[[old]]
}
df <- merge(df, shan_fpc[, c("animal_id", unname(shan_rename))],
            by = "animal_id", all.x = TRUE)

# --- Save --------------------------------------------------------------------

out_file <- "results/metrics/maturation_all_metrics_per_animal.csv"
write.csv(df, out_file, row.names = FALSE)
cat(sprintf("Saved: %s (%d animals, %d columns)\n", out_file, nrow(df), ncol(df)))

# --- Summary by line ---------------------------------------------------------

cat("\n=== SEN vs TOL summary ===\n")
numeric_cols <- names(df)[sapply(df, is.numeric) & !names(df) %in% "n_obs"]
for (col in numeric_cols) {
  vals_sen <- df[[col]][df$line == "SEN"]
  vals_tol <- df[[col]][df$line == "TOL"]
  if (all(is.na(vals_sen)) || all(is.na(vals_tol))) next
  tt <- tryCatch(t.test(vals_sen, vals_tol), error = function(e) NULL)
  if (is.null(tt)) next
  cat(sprintf("  %-30s  SEN=%.4f  TOL=%.4f  p=%.4f%s\n",
              col,
              mean(vals_sen, na.rm = TRUE), mean(vals_tol, na.rm = TRUE),
              tt$p.value,
              ifelse(tt$p.value < 0.05, " *", ifelse(tt$p.value < 0.1, " .", ""))))
}

cat(sprintf("\nGompertz converged: %d/%d animals (%.0f%%)\n",
            sum(df$gompertz_converged, na.rm = TRUE), nrow(df),
            100 * mean(df$gompertz_converged, na.rm = TRUE)))

# --- Cross-class correlation matrix ------------------------------------------

cat("\n=== Cross-class Spearman correlations ===\n")
key_metrics <- c("maz_mean", "maz_slope", "bc_to_adult_slope",
                 "aitchison_to_adult_slope", "rel_maturity_w24",
                 "gompertz_k", "gompertz_ti", "volatility_mean",
                 "volatility_slope", "shannon_slope", "shannon_cv")
cor_mat <- cor(df[, key_metrics], use = "pairwise.complete.obs", method = "spearman")
for (i in seq_len(nrow(cor_mat))) {
  cat(sprintf("  %-25s %s\n", rownames(cor_mat)[i],
              paste(sprintf("%6.2f", cor_mat[i, ]), collapse = " ")))
}

cat("\nDone.\n")

# PART 2: Time-to-threshold descriptors (t55, t60, t65, t70, t75)
setwd(here::here())
suppressPackageStartupMessages({
  library(arrow); library(asreml)
})

# --- Load PACE-imputed wide matrix -----------------------------------------

wide <- read_parquet("data/processed/mag_matrix_wide_pace_imputed.parquet")
animal_ids <- as.character(wide$animal_id)
feat_cols <- setdiff(names(wide), "animal_id")
mat <- as.matrix(wide[, feat_cols])
rownames(mat) <- animal_ids
feature_names <- unique(sapply(strsplit(feat_cols, "_(?=W[0-9]+$)", perl=TRUE),
                               function(x) x[1]))
tp_labels <- c("W3","W6","W9","W12","W15","W17","W21","W24")
weeks <- c(W3=3, W6=6, W9=9, W12=12, W15=15, W17=17, W21=21, W24=24)

# --- Build maturity = 1 - BC(profile, W24-centroid) ------------------------

w24 <- colMeans(mat[, paste0(feature_names, "_W24")], na.rm=TRUE)
w24_prop <- w24 / sum(w24)

bc_df <- list()
for (a in animal_ids) for (tp in tp_labels) {
  cols <- paste0(feature_names, "_", tp)
  v <- mat[a, cols]
  if (all(is.na(v))) next
  p <- v / sum(v)
  bc_df[[length(bc_df)+1]] <- data.frame(
    animal_id=a, w=weeks[[tp]],
    maturity = 1 - sum(abs(p - w24_prop))/2,
    stringsAsFactors=FALSE)
}
bc_df <- do.call(rbind, bc_df)

# --- Per-animal time-to-threshold via linear interpolation -----------------

threshes <- c(0.55, 0.60, 0.65, 0.70, 0.75)
results <- data.frame(animal_id=unique(bc_df$animal_id), stringsAsFactors=FALSE)

cat("Computing time-to-threshold descriptors for", nrow(results), "animals...\n")
for (a in unique(bc_df$animal_id)) {
  d <- bc_df[bc_df$animal_id==a, ]
  d <- d[order(d$w), ]
  if (nrow(d) < 5) next
  # Use full W3-W24 range (PACE-imputed trajectory is smooth)
  t_grid <- seq(min(d$w), max(d$w), by=0.1)
  m_interp <- approx(d$w, d$maturity, xout=t_grid, rule=2)$y
  for (thr in threshes) {
    idx <- which(m_interp >= thr)[1]
    val <- if (is.na(idx)) NA else t_grid[idx]
    results[results$animal_id==a, sprintf("t%02d", round(thr*100))] <- val
  }
}

# Attach line and litter
ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors=FALSE)
ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors=FALSE)
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
litter_lookup <- setNames(as.character(ped$litter_id), ped$animal_id)
line_lookup <- setNames(ped$line, ped$animal_id)
results$litter_id <- litter_lookup[results$animal_id]
results$line <- line_lookup[results$animal_id]

write.csv(results, "results/metrics/maturation_threshold_metrics_imputed.csv",
          row.names=FALSE)
cat("Saved: results/metrics/maturation_threshold_metrics_imputed.csv\n\n")

# --- Animal model h² per threshold -----------------------------------------

results$animal_id_f <- factor(results$animal_id, levels=ped3$id)
results$litter_id_f <- factor(results$litter_id)

fit_h2 <- function(y, label) {
  d <- data.frame(animal_id=results$animal_id_f, litter_id=results$litter_id_f, y=y)
  d <- d[!is.na(d$y), ]
  if (nrow(d) < 10 || sd(d$y) < 1e-8) {
    return(data.frame(metric=label, n=nrow(d), h2=NA, h2_se=NA, c2=NA,
                      lrt_p=NA, stringsAsFactors=FALSE))
  }
  d$y <- scale(d$y)[,1]
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data=d, na.action=na.method(y="include"),
           trace=FALSE, maxit=100))), error=function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(metric=label, n=nrow(d), h2=NA, h2_se=NA, c2=NA,
                      lrt_p=NA, stringsAsFactors=FALSE))
  }
  vc <- summary(fit)$varcomp
  pat_a <- "vm\\(animal_id"
  s2_a <- max(vc[grep(pat_a, rownames(vc))[1], "component"], 0)
  s2_c <- max(vc[which(rownames(vc)=="litter_id"), "component"], 0)
  s2_e <- vc[grep("units", rownames(vc))[1], "component"]
  h2 <- s2_a / (s2_a+s2_c+s2_e)
  c2 <- s2_c / (s2_a+s2_c+s2_e)
  h2_se <- tryCatch(vpredict(fit, h2 ~ V2/(V1+V2+V3))$SE[1], error=function(e) NA)
  fit0 <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id, residual = ~ idv(units), data=d,
           na.action=na.method(y="include"), trace=FALSE, maxit=100))),
    error=function(e) NULL)
  lrt_p <- if (is.null(fit0)) NA else
    0.5 * pchisq(max(2*(fit$loglik-fit0$loglik), 0), df=1, lower.tail=FALSE)
  data.frame(metric=label, n=nrow(d), h2=h2, h2_se=h2_se, c2=c2,
             lrt_p=lrt_p, stringsAsFactors=FALSE)
}

cat("=== h² per threshold ===\n")
h2_rows <- list()
for (col in paste0("t", c("55","60","65","70","75"))) {
  r <- fit_h2(results[[col]], col)
  h2_rows[[col]] <- r
  cat(sprintf("%-5s n=%-3d h2=%.3f (%.3f)  c2=%.3f  LRT_p=%.3f\n",
              r$metric, r$n, r$h2, r$h2_se, r$c2, r$lrt_p))
}
h2_df <- do.call(rbind, h2_rows)
write.csv(h2_df, "results/metrics/maturation_h2_threshold.csv", row.names=FALSE)
cat("\nSaved: results/metrics/maturation_h2_threshold.csv\n\n")

# --- SEN/TOL line tests -----------------------------------------------------

cat("=== SEN/TOL line tests ===\n")
for (col in paste0("t", c("55","60","65","70","75"))) {
  s <- results[[col]][results$line=="SEN"]; t <- results[[col]][results$line=="TOL"]
  s <- s[!is.na(s)]; t <- t[!is.na(t)]
  if (length(s) < 5 || length(t) < 5 || sd(s) < 1e-8 || sd(t) < 1e-8) {
    cat(sprintf("%-5s skip (constant or low n)\n", col)); next
  }
  cat(sprintf("%-5s SEN=%5.2f (n=%d, sd=%.2f)  TOL=%5.2f (n=%d, sd=%.2f)  p=%.4f\n",
              col, mean(s), length(s), sd(s), mean(t), length(t), sd(t),
              t.test(s, t)$p.value))
}

cat("\nDone.\n")

# PART 3: Extended alpha-diversity metrics (Shannon, Pielou, Simpson, invSimpson, Faith's PD)
setwd(here::here())
suppressPackageStartupMessages({
  library(arrow)
  library(vegan)     # diversity, specnumber
  library(picante)   # pd (Faith's phylogenetic diversity)
  library(ape)       # read.tree
  library(asreml)
})

# --- Load matrices ----------------------------------------------------------

wide <- read_parquet("data/processed/mag_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
cols <- setdiff(names(wide), "animal_id")
mat_wide <- as.matrix(wide[, cols])
rownames(mat_wide) <- animal_ids

# Convert wide (animals × MAG_W) to long (sample × MAG)
parts <- strsplit(cols, "_(?=W[0-9]+$)", perl = TRUE)
feat_per_col <- sapply(parts, function(x) x[1])
tp_per_col   <- sapply(parts, function(x) x[2])
features  <- unique(feat_per_col)
tps       <- sort(unique(tp_per_col))

# Build sample × feature matrix; sample id = "<animal>_<tp>"
samp_mat <- matrix(NA_real_, nrow = length(animal_ids) * length(tps),
                   ncol = length(features),
                   dimnames = list(NULL, features))
samp_meta <- expand.grid(animal_id = animal_ids, tp = tps,
                         stringsAsFactors = FALSE)
samp_meta$sample_id <- paste0(samp_meta$animal_id, "_", samp_meta$tp)
rownames(samp_mat) <- samp_meta$sample_id

for (i in seq_along(cols)) {
  feat <- feat_per_col[i]; tp <- tp_per_col[i]
  for (a in animal_ids) {
    samp_mat[paste0(a, "_", tp), feat] <- mat_wide[a, cols[i]]
  }
}

# Drop samples that are entirely NA (animal not sampled at that timepoint)
keep <- rowSums(is.na(samp_mat)) < ncol(samp_mat)
samp_mat <- samp_mat[keep, , drop = FALSE]
samp_meta <- samp_meta[keep, ]
cat(sprintf("Sample matrix: %d samples × %d features\n",
            nrow(samp_mat), ncol(samp_mat)))

# Replace NAs (drop-sample) with 0 for diversity calc on remaining cells
samp_mat[is.na(samp_mat)] <- 0

# --- Phylo tree for Faith's PD ----------------------------------------------

tree <- tryCatch(read.tree("data/raw/gtdbtk.bac120.nwk"), error = function(e) NULL)
have_tree <- !is.null(tree)
if (have_tree) {
  # Match tree tips to feature names
  common <- intersect(tree$tip.label, colnames(samp_mat))
  cat(sprintf("Tree tips: %d, intersect with features: %d\n",
              length(tree$tip.label), length(common)))
  if (length(common) >= 50) {
    tree_pruned <- ape::keep.tip(tree, common)
    samp_for_pd <- samp_mat[, common, drop = FALSE]
  } else {
    have_tree <- FALSE
    cat("Tree intersection too small; skipping Faith's PD.\n")
  }
} else {
  cat("Phylo tree missing; skipping Faith's PD.\n")
}

# --- Compute alpha-diversity metrics per sample -----------------------------

# Convert to relative abundance for Simpson/Shannon (counts ok for richness/Pielou)
n_obs <- ncol(samp_mat); n_samp <- nrow(samp_mat)

samp_rel <- samp_mat / pmax(rowSums(samp_mat), 1)
shannon  <- diversity(samp_rel, index = "shannon")
simpson  <- diversity(samp_rel, index = "simpson")  # 1 - sum(p^2)
invsimp  <- diversity(samp_rel, index = "invsimpson")
richness <- specnumber(samp_mat)
pielou   <- ifelse(richness > 1, shannon / log(richness), NA)

if (have_tree) {
  pd_obj <- pd(samp_for_pd, tree_pruned, include.root = FALSE)
  faith_pd <- pd_obj$PD
  faith_sr <- pd_obj$SR
} else {
  faith_pd <- NA; faith_sr <- NA
}

per_sample <- data.frame(
  sample_id  = samp_meta$sample_id,
  animal_id  = samp_meta$animal_id,
  tp         = samp_meta$tp,
  shannon    = shannon,
  simpson    = simpson,
  invsimpson = invsimp,
  richness   = richness,
  pielou     = pielou,
  faith_pd   = faith_pd,
  stringsAsFactors = FALSE
)
per_sample$week <- as.integer(sub("^W", "", per_sample$tp))

write.csv(per_sample,
          "results/metrics/alpha_diversity_extended_per_sample.csv",
          row.names = FALSE)
cat("Wrote per-sample metrics.\n")

# --- Per-animal summaries (mean, slope, W24, CV) ----------------------------

per_animal_summary <- function(df, col) {
  fit <- tryCatch(stats::lm(df[[col]] ~ df$week), error = function(e) NULL)
  list(
    mean_v = mean(df[[col]], na.rm = TRUE),
    slope  = if (!is.null(fit) && nrow(coef(summary(fit))) >= 2) coef(fit)[2] else NA,
    cv     = sd(df[[col]], na.rm = TRUE) / max(abs(mean(df[[col]], na.rm = TRUE)), 1e-9),
    w24    = if (any(df$tp == "W24")) df[[col]][df$tp == "W24"][1] else NA
  )
}

metrics_to_summarise <- c("shannon", "simpson", "invsimpson", "richness", "pielou", "faith_pd")
metrics_to_summarise <- metrics_to_summarise[metrics_to_summarise %in% names(per_sample)]

per_animal <- data.frame(animal_id = unique(per_sample$animal_id))
for (m in metrics_to_summarise) {
  vals_mean <- vals_slope <- vals_cv <- vals_w24 <- numeric(nrow(per_animal))
  for (i in seq_len(nrow(per_animal))) {
    sub <- per_sample[per_sample$animal_id == per_animal$animal_id[i], ]
    s <- per_animal_summary(sub, m)
    vals_mean[i]  <- s$mean_v
    vals_slope[i] <- s$slope
    vals_cv[i]    <- s$cv
    vals_w24[i]   <- s$w24
  }
  per_animal[[paste0(m, "_mean")]]  <- vals_mean
  per_animal[[paste0(m, "_slope")]] <- vals_slope
  per_animal[[paste0(m, "_cv")]]    <- vals_cv
  per_animal[[paste0(m, "_w24")]]   <- vals_w24
}

write.csv(per_animal,
          "results/metrics/alpha_diversity_extended_per_animal.csv",
          row.names = FALSE)
cat(sprintf("Wrote per-animal metrics: %d animals × %d cols\n",
            nrow(per_animal), ncol(per_animal)))

# --- h² per per-animal alpha-diversity summary (animal model) --------------

ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

per_animal <- merge(per_animal, ped[, c("animal_id", "litter_id", "line")],
                    by = "animal_id", all.x = TRUE)
per_animal$animal_id <- factor(per_animal$animal_id, levels = ped3$id)
per_animal$litter_id <- factor(per_animal$litter_id)

extract_pieces <- function(fit) {
  vc <- summary(fit)$varcomp
  rn <- rownames(vc)
  s2_A <- {a <- grep("vm\\(animal_id", rn); if (length(a)) max(vc[a[1], "component"], 0) else 0}
  s2_c <- {c <- which(rn == "litter_id");   if (length(c)) max(vc[c, "component"], 0) else 0}
  s2_e <- {e <- grep("units", rn);          if (length(e)) vc[e[1], "component"] else NA}
  total <- s2_A + s2_c + s2_e
  h2vp <- tryCatch(vpredict(fit, h2 ~ V2 / (V1 + V2 + V3)), error = function(e) NULL)
  c2vp <- tryCatch(vpredict(fit, c2 ~ V1 / (V1 + V2 + V3)), error = function(e) NULL)
  list(
    s2_A = s2_A, s2_c = s2_c, s2_e = s2_e, total = total,
    h2 = if (!is.null(h2vp)) h2vp$Estimate[1] else s2_A / total,
    h2_se = if (!is.null(h2vp)) h2vp$SE[1] else NA,
    c2 = if (!is.null(c2vp)) c2vp$Estimate[1] else s2_c / total,
    c2_se = if (!is.null(c2vp)) c2vp$SE[1] else NA,
    loglik = fit$loglik
  )
}

h2_results <- data.frame()
target_cols <- setdiff(names(per_animal),
                       c("animal_id", "line", "litter_id"))

for (col in target_cols) {
  y_raw <- per_animal[[col]]
  n_obs <- sum(!is.na(y_raw))
  if (n_obs < 30) next
  if (sd(y_raw, na.rm = TRUE) < 1e-9) next

  df_fit <- data.frame(
    y = as.numeric(scale(y_raw)),
    animal_id = per_animal$animal_id,
    litter_id = per_animal$litter_id
  )
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = df_fit, trace = FALSE, maxit = 100,
           na.action = na.method(y = "include"))
  )), error = function(e) NULL)
  if (is.null(fit)) next
  p <- extract_pieces(fit)

  h2_results <- rbind(h2_results, data.frame(
    metric = col, n = n_obs,
    s2_A = p$s2_A, s2_c = p$s2_c, s2_e = p$s2_e,
    h2 = p$h2, h2_se = p$h2_se,
    c2 = p$c2, c2_se = p$c2_se,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  %-25s n=%d  h²=%.3f±%.3f  c²=%.3f±%.3f\n",
              col, n_obs, p$h2, p$h2_se, p$c2, p$c2_se))
}

write.csv(h2_results,
          "results/metrics/alpha_diversity_extended_h2.csv",
          row.names = FALSE)
cat("\n=== outputs ===\n")
cat("results/metrics/alpha_diversity_extended_per_sample.csv\n")
cat("results/metrics/alpha_diversity_extended_per_animal.csv\n")
cat("results/metrics/alpha_diversity_extended_h2.csv\n")
