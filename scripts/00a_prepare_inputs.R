#!/usr/bin/env Rscript
#________________________________________
# - File Name : 00a_prepare_inputs.R
# - Purpose : Reshape the Zenodo deposit into the layout the pipeline expects
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
# Sat Aug 01 13:44:36 2026 ------------------------------

# Reads the Zenodo files as deposited (flat, under data/) and writes the
# per-animal pedigree and the two wide feature matrices the pipeline reads.
# Column order is load-bearing: scripts derive the timepoint sequence from it,
# so timepoints run chronologically (outer) and features ascending (inner).

setwd(here::here())
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr)
})

dir.create("data/metadata",  recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
# 10_feature_level_h2.R reads the KO -> module map here; without it, it fetches
# the current KEGG release, which differs from the deposited one.
stopifnot(file.copy("data/ko_module_map.csv", "data/processed/ko_module_map.csv", overwrite = TRUE))
# Derived inputs deposited under data/derived/ (Zenodo >= 1.1.0) are read by the
# scripts from results/metrics/.
dir.create("results/metrics", recursive = TRUE, showWarnings = FALSE)
derived <- list.files("data/derived", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(derived) == 12,
          file.copy(derived, "results/metrics", overwrite = TRUE))

TPS <- c("W3", "W6", "W9", "W12", "W15", "W17", "W21", "W24")

# --- Per-animal pedigree from the per-sample metadata ------------------------
meta <- read.csv("data/sample_metadata.csv", stringsAsFactors = FALSE)
ped <- meta %>%
  select(animal_id, litter_id, sire_id, age_days_w3, line, pedigree_inferred) %>%
  distinct() %>%
  arrange(animal_id)
stopifnot(!anyDuplicated(ped$animal_id))
write.csv(ped, "data/metadata/pedigree.csv", row.names = FALSE)
cat(sprintf("Wrote data/metadata/pedigree.csv (%d animals)\n", nrow(ped)))

# --- Long -> wide feature matrices -------------------------------------------
to_wide <- function(long_path, feature_col, out_path) {
  lg <- read.csv(long_path, stringsAsFactors = FALSE)
  stopifnot(all(lg$timepoint %in% TPS))
  lg$col <- paste0(lg[[feature_col]], "_", lg$timepoint)
  feats <- sort(unique(lg[[feature_col]]))
  col_order <- as.vector(outer(feats, TPS, function(f, t) paste0(f, "_", t)))

  w <- lg %>%
    select(animal_id, col, count) %>%
    pivot_wider(names_from = col, values_from = count) %>%
    arrange(animal_id)

  missing <- setdiff(col_order, names(w))
  if (length(missing)) w[missing] <- NA_real_
  w <- w[, c(col_order, "animal_id")]

  write_parquet(w, out_path)
  cat(sprintf("Wrote %s (%d animals x %d feature-timepoint columns)\n",
              out_path, nrow(w), length(col_order)))
}

to_wide("data/mag_abundance_long.csv.gz", "mag_id",
        "data/processed/mag_matrix_wide.parquet")
to_wide("data/ko_coverage_long.csv.gz", "ko_id",
        "data/processed/ko_matrix_wide.parquet")
