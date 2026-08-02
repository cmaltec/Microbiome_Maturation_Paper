#!/usr/bin/env Rscript
#________________________________________
# - File Name : 09_selection_response.R
# - Purpose : Expected per-generation response to mass and sire-BLUP selection
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
# Mon May 25 10:09:00 2026 ------------------------------

setwd(here::here())
suppressPackageStartupMessages({ library(dplyr) })

h2_full <- read.csv("results/metrics/maturation_h2_animal_model.csv",
                    stringsAsFactors = FALSE)
h2_w6   <- read.csv("results/metrics/maturation_h2_w6_scalars.csv",
                    stringsAsFactors = FALSE)
h2_thr  <- read.csv("results/metrics/maturation_h2_threshold.csv",
                    stringsAsFactors = FALSE)
metrics <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                    stringsAsFactors = FALSE)
thr     <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv",
                    stringsAsFactors = FALSE)
acc     <- read.csv("results/metrics/maturation_sire_blup_accuracy.csv",
                    stringsAsFactors = FALSE)

# Trait-level mean PEV-based accuracy across the 7 sires (script 08 output).
r_by_trait <- acc %>%
  group_by(metric) %>%
  summarise(r_mean = mean(accuracy, na.rm = TRUE),
            r_min  = suppressWarnings(min(accuracy, na.rm = TRUE)),
            r_max  = suppressWarnings(max(accuracy, na.rm = TRUE)),
            .groups = "drop")

get_r_sire <- function(metric) {
  rec <- r_by_trait[r_by_trait$metric == metric, ]
  if (nrow(rec) == 0 || !is.finite(rec$r_mean)) return(NA_real_)
  rec$r_mean
}

# Build the per-trait input table. For the W3-robust subset we prefer the
# W6+ refit (more conservative SE); for the threshold metrics we pull from
# maturation_h2_threshold.csv.
# Focal trait set (the W3-robust headlines plus the threshold descriptors)
focal <- list(
  list(metric = "maz_mean",            source = "obs full",       h2_src = h2_full,  sigma_src = metrics,
       sigma_col = "maz_mean"),
  list(metric = "bc_to_adult_slope",   source = "obs W6+ refit",  h2_src = h2_w6,    sigma_src = metrics,
       sigma_col = "bc_to_adult_slope"),
  list(metric = "aitchison_slope",     source = "obs full (imp)", h2_src = h2_full,  sigma_src = metrics,
       sigma_col = "aitchison_to_adult_slope"),
  list(metric = "shannon_slope",       source = "obs W6+ refit",  h2_src = h2_w6,    sigma_src = metrics,
       sigma_col = "shannon_slope"),
  list(metric = "intra_bc_var",        source = "obs W6+ refit",  h2_src = h2_w6,    sigma_src = metrics,
       sigma_col = "intra_bc_var"),
  list(metric = "BC_FPC2",             source = "obs full",       h2_src = h2_full,  sigma_src = metrics,
       sigma_col = "BC_FPC2"),
  list(metric = "t70",                 source = "PACE-imputed",   h2_src = h2_thr,   sigma_src = thr,
       sigma_col = "t70"),
  list(metric = "t75",                 source = "PACE-imputed",   h2_src = h2_thr,   sigma_src = thr,
       sigma_col = "t75")
)

# selection intensities
i_20 <- 1.40   # top 20%
i_05 <- 2.06   # top 5%

rows <- data.frame()
for (rec in focal) {
  h2_row <- rec$h2_src[rec$h2_src$metric == rec$metric, ]
  if (nrow(h2_row) == 0) {
    # h2_w6 uses h2_w6 column instead of h2; try fallback
    if ("h2_w6" %in% names(rec$h2_src)) {
      h2_row <- rec$h2_src[rec$h2_src$metric == rec$metric, ]
    }
  }
  if (nrow(h2_row) == 0) next
  if ("h2_w6" %in% names(h2_row) && !"h2" %in% names(h2_row)) {
    h2_val <- h2_row$h2_w6
    h2_se_val <- h2_row$h2_w6_se
  } else {
    h2_val <- h2_row$h2
    h2_se_val <- h2_row$h2_se
  }
  h2_val <- pmax(h2_val, 0)
  sigma_col <- rec$sigma_col
  sigma_P <- if (sigma_col %in% names(rec$sigma_src)) {
    sd(rec$sigma_src[[sigma_col]], na.rm = TRUE)
  } else NA_real_

  r_sire_trait <- get_r_sire(rec$metric)

  dG_mass_20 <- i_20 * h2_val
  dG_mass_05 <- i_05 * h2_val
  dG_sire_20 <- if (!is.na(r_sire_trait)) i_20 * r_sire_trait * sqrt(h2_val) else NA_real_
  dG_sire_05 <- if (!is.na(r_sire_trait)) i_05 * r_sire_trait * sqrt(h2_val) else NA_real_

  se_mass_20 <- i_20 * h2_se_val
  se_mass_05 <- i_05 * h2_se_val
  se_sire_20 <- if (h2_val > 1e-6 && !is.na(r_sire_trait)) i_20 * r_sire_trait * h2_se_val / (2 * sqrt(h2_val)) else NA_real_
  se_sire_05 <- if (h2_val > 1e-6 && !is.na(r_sire_trait)) i_05 * r_sire_trait * h2_se_val / (2 * sqrt(h2_val)) else NA_real_

  rows <- rbind(rows, data.frame(
    metric = rec$metric, source = rec$source,
    h2 = h2_val, h2_se = h2_se_val, sigma_P = sigma_P, r_sire_blup = r_sire_trait,
    dG_mass_i14   = dG_mass_20, dG_mass_i14_lo = pmax(0, dG_mass_20 - 1.96 * se_mass_20),
                                  dG_mass_i14_hi = dG_mass_20 + 1.96 * se_mass_20,
    dG_mass_i21   = dG_mass_05, dG_mass_i21_lo = pmax(0, dG_mass_05 - 1.96 * se_mass_05),
                                  dG_mass_i21_hi = dG_mass_05 + 1.96 * se_mass_05,
    dG_sire_i14   = dG_sire_20, dG_sire_i14_lo = pmax(0, dG_sire_20 - 1.96 * se_sire_20),
                                  dG_sire_i14_hi = dG_sire_20 + 1.96 * se_sire_20,
    dG_sire_i21   = dG_sire_05, dG_sire_i21_lo = pmax(0, dG_sire_05 - 1.96 * se_sire_05),
                                  dG_sire_i21_hi = dG_sire_05 + 1.96 * se_sire_05,
    stringsAsFactors = FALSE))
}

write.csv(rows, "results/metrics/maturation_selection_response.csv",
          row.names = FALSE)
cat(sprintf("Wrote results/metrics/maturation_selection_response.csv  (%d traits)\n",
            nrow(rows)))
cat("Per-trait mean sire-BLUP accuracy (PEV-based, from script 08):\n")
print(r_by_trait, row.names = FALSE)
print(rows[, c("metric", "source", "h2", "h2_se", "sigma_P", "r_sire_blup",
               "dG_mass_i14", "dG_sire_i14")], row.names = FALSE)
