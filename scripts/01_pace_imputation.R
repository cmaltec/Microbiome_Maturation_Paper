#!/usr/bin/env Rscript
#________________________________________
# - File Name : 01_pace_imputation.R
# - Purpose : PACE (fdapace FPCA) imputation of per-sample MAG and Shannon trajectories
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
suppressPackageStartupMessages({library(fdapace); library(arrow)})

# === Load data ===
wide <- read_parquet("data/processed/mag_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
feat_cols <- setdiff(names(wide), "animal_id")
mat <- as.matrix(wide[, feat_cols]); rownames(mat) <- animal_ids
log_mat <- log1p(mat)

col_parts <- strsplit(feat_cols, "_(?=W[0-9]+$)", perl = TRUE)
feature_names <- unique(sapply(col_parts, function(x) x[1]))
tp_labels <- unique(sapply(col_parts, function(x) x[2]))
week_map <- c(W3=3, W6=6, W9=9, W12=12, W15=15, W17=17, W21=21, W24=24)
all_weeks <- as.numeric(week_map[tp_labels])

ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

n_animals <- length(animal_ids)
n_features <- length(feature_names)
n_tp <- length(tp_labels)
n_pcs <- 20

cat(sprintf("Data: %d animals, %d features, %d timepoints\n", n_animals, n_features, n_tp))

# === Build sample-level data + CLR + PCA ===
rows <- list()
for (a in animal_ids) for (tp in tp_labels) {
  if (!is.na(log_mat[a, paste0(feature_names[1], "_", tp)])) {
    rows[[length(rows) + 1]] <- data.frame(
      animal_id = a, tp = tp, week = week_map[[tp]],
      stringsAsFactors = FALSE)
  }
}
sample_df <- do.call(rbind, rows)
n_samples <- nrow(sample_df)

count_mat <- matrix(NA_real_, n_samples, n_features)
for (fi in seq_along(feature_names)) for (ri in seq_len(n_samples)) {
  count_mat[ri, fi] <- log_mat[sample_df$animal_id[ri],
                                paste0(feature_names[fi], "_", sample_df$tp[ri])]
}
clr_mat <- count_mat - rowMeans(count_mat, na.rm = TRUE)
pca_fit <- prcomp(clr_mat, center = TRUE, scale. = FALSE, rank. = n_pcs)
pc_scores <- pca_fit$x

cat(sprintf("Observed samples: %d, PCs: %d\n", n_samples, n_pcs))
cat(sprintf("Missingness: %.1f%% of animal×tp cells\n",
            100 * (1 - n_samples / (n_animals * n_tp))))

# Precompute row-mean offsets for back-transform
sample_row_mean <- rowMeans(count_mat, na.rm = TRUE)
tp_offset_mat <- matrix(NA_real_, n_animals, n_tp,
                        dimnames = list(animal_ids, tp_labels))
for (si in seq_len(n_samples)) {
  tp_offset_mat[sample_df$animal_id[si], sample_df$tp[si]] <- sample_row_mean[si]
}

# === Fit FPCA per PC on ALL observed data, predict at ALL timepoints ===
cat("\nFitting FPCA on full data...\n")

# Output: predicted PC scores for every animal × every timepoint
pred_pc_mat <- array(NA_real_, dim = c(n_animals, n_tp, n_pcs),
                     dimnames = list(animal_ids, tp_labels, paste0("PC", 1:n_pcs)))

for (pci in seq_len(n_pcs)) {
  # Build Ly/Lt per animal (observed data only)
  Ly <- list(); Lt <- list(); animal_order <- character(0)
  for (a in animal_ids) {
    a_rows <- which(sample_df$animal_id == a)
    if (length(a_rows) == 0) next
    Ly[[length(Ly) + 1]] <- pc_scores[a_rows, pci]
    Lt[[length(Lt) + 1]] <- sample_df$week[a_rows]
    animal_order <- c(animal_order, a)
  }

  # Fit FPCA
  fpca_res <- tryCatch(
    suppressWarnings(FPCA(Ly, Lt,
      optns = list(dataType = "Sparse", methodMuCovEst = "smooth",
                   methodSelectK = 3, useBinnedData = "OFF",
                   error = TRUE, verbose = FALSE))),
    error = function(e) { cat(sprintf("  PC%d FPCA failed: %s\n", pci, e$message)); NULL })

  if (is.null(fpca_res)) next

  # Predict complete trajectories for all animals
  pred <- tryCatch(
    predict(fpca_res, newLy = Ly, newLt = Lt,
            K = min(3, length(fpca_res$lambda)),
            xiMethod = "CE")$predCurves,
    error = function(e) { cat(sprintf("  PC%d predict failed: %s\n", pci, e$message)); NULL })

  if (is.null(pred)) next

  # pred is n_curves × length(workGrid); interpolate to our 8 timepoints
  wg <- fpca_res$workGrid
  for (ai in seq_along(animal_order)) {
    a <- animal_order[ai]
    for (ti in seq_along(tp_labels)) {
      wi <- which.min(abs(wg - all_weeks[ti]))
      pred_pc_mat[a, tp_labels[ti], pci] <- pred[ai, wi]
    }
  }

  # Animals with no observed data: use population mean function
  missing_animals <- setdiff(animal_ids, animal_order)
  if (length(missing_animals) > 0) {
    for (a in missing_animals) {
      for (ti in seq_along(tp_labels)) {
        wi <- which.min(abs(wg - all_weeks[ti]))
        pred_pc_mat[a, tp_labels[ti], pci] <- fpca_res$mu[wi]
      }
    }
  }

  if (pci %% 5 == 0 || pci == 1) cat(sprintf("  PC%d done\n", pci))
}

cat("FPCA prediction complete.\n")

# === Back-transform: PC space → CLR → log space ===
cat("Back-transforming to log-count space...\n")

imputed_wide <- matrix(NA_real_, n_animals, length(feat_cols),
                       dimnames = list(animal_ids, feat_cols))

for (ai in seq_along(animal_ids)) {
  a <- animal_ids[ai]
  for (ti in seq_along(tp_labels)) {
    tp <- tp_labels[ti]
    pred_pcs <- pred_pc_mat[a, tp, ]
    if (any(is.na(pred_pcs))) next

    # PC → CLR
    imputed_clr <- as.numeric(pred_pcs %*% t(pca_fit$rotation)) + pca_fit$center

    # Offset: use animal's own mean if observed, population mean at tp otherwise
    offset <- tp_offset_mat[a, tp]
    if (is.na(offset)) {
      offset <- mean(tp_offset_mat[, tp], na.rm = TRUE)
    }

    imputed_log <- imputed_clr + offset
    cols <- paste0(feature_names, "_", tp)
    imputed_wide[a, cols] <- imputed_log
  }
}

# Convert back to count space
imputed_counts <- expm1(imputed_wide)
imputed_counts[imputed_counts < 0] <- 0

# === Diagnostics ===
n_total_cells <- n_animals * n_tp
n_filled <- sum(!is.na(imputed_wide)) / n_features
cat(sprintf("\nImputation summary:\n"))
cat(sprintf("  Original: %d/%d animal×tp cells observed (%.1f%%)\n",
            n_samples, n_total_cells, 100 * n_samples / n_total_cells))
cat(sprintf("  Imputed:  %.0f/%d animal×tp cells filled (%.1f%%)\n",
            n_filled, n_total_cells, 100 * n_filled / n_total_cells))
cat(sprintf("  Remaining NA: %.0f cells\n", sum(is.na(imputed_wide)) / n_features))

# Compare observed vs PACE-smoothed values
obs_vs_pace <- data.frame()
for (ri in seq_len(n_samples)) {
  a <- sample_df$animal_id[ri]; tp <- sample_df$tp[ri]
  cols <- paste0(feature_names, "_", tp)
  obs <- log_mat[a, cols]
  imp <- imputed_wide[a, cols]
  if (any(is.na(imp))) next
  obs_vs_pace <- rbind(obs_vs_pace, data.frame(
    animal_id = a, tp = tp,
    cor_pearson = cor(obs, imp, use = "complete.obs"),
    rmse = sqrt(mean((obs - imp)^2, na.rm = TRUE)),
    stringsAsFactors = FALSE))
}
cat(sprintf("  Observed vs PACE-smoothed: median r=%.4f, median RMSE=%.4f\n",
            median(obs_vs_pace$cor_pearson, na.rm = TRUE),
            median(obs_vs_pace$rmse, na.rm = TRUE)))

# === Save ===
out_df <- data.frame(animal_id = animal_ids, imputed_counts, check.names = FALSE)
out_file <- "data/processed/mag_matrix_wide_pace_imputed.parquet"
write_parquet(out_df, out_file)
cat(sprintf("\nSaved: %s\n", out_file))

write.csv(obs_vs_pace, "results/metrics/pace_full_fit_summary.csv", row.names = FALSE)
cat("Saved: results/metrics/pace_full_fit_summary.csv\n")

cat("\nDone.\n")
