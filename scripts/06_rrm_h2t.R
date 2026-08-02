#!/usr/bin/env Rscript
#________________________________________
# - File Name : 06_rrm_h2t.R
# - Purpose : Random-regression animal model for community-level h2(t)
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
suppressPackageStartupMessages({
  library(asreml); library(dplyr)
})

cur  <- read.csv("results/metrics/maturation_curves_per_sample.csv",
                 stringsAsFactors = FALSE)
maz  <- read.csv("results/metrics/maz_per_sample.csv", stringsAsFactors = FALSE)
ped  <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)

cur$week <- as.integer(sub("^W", "", cur$tp))
cur <- merge(cur, ped[, c("animal_id", "litter_id")], by = "animal_id", all.x = TRUE)
cur <- merge(cur, maz[, c("animal_id", "tp", "maz")], by = c("animal_id", "tp"), all.x = TRUE)
schedule_mean <- mean(c(3, 6, 9, 12, 15, 17, 21, 24))
cur$week_c <- cur$week - schedule_mean
cur$animal_id <- factor(as.character(cur$animal_id), levels = ped3$id)
cur$litter_id <- factor(as.character(cur$litter_id))

reload_ainv <- function() {
  A <- readRDS("data/processed/A_matrix_animals.rds")
  Ai <- solve(A)
  attr(Ai, "rowNames") <- rownames(A)
  attr(Ai, "colNames") <- colnames(A)
  Ai
}

fit_rrm <- function(trait_col) {
  Ainv <- reload_ainv()
  df <- cur[!is.na(cur[[trait_col]]), ]
  df$y <- as.numeric(scale(df[[trait_col]]))
  suppressWarnings(suppressMessages(
    asreml(
      fixed    = y ~ 1 + week_c,
      random   = ~ str( ~ vm(animal_id, Ainv) + week_c:vm(animal_id, Ainv),
                        ~ us(2):id(78)) + litter_id,
      residual = ~ idv(units),
      data     = df, trace = FALSE, maxit = 200,
      na.action = na.method(y = "include"))
  ))
}

extract_components <- function(fit) {
  vc <- summary(fit)$varcomp
  pick <- function(pat) {
    idx <- grep(pat, rownames(vc))
    if (length(idx) == 0) return(c(NA, NA))
    c(vc[idx[1], "component"], vc[idx[1], "std.error"])
  }
  list(
    v_int  = pick("us\\(2\\)_1:1$"),
    v_slo  = pick("us\\(2\\)_2:2$"),
    v_cov  = pick("us\\(2\\)_2:1$"),
    v_litt = pick("^litter_id"),
    v_res  = pick("units!units$"))
}

h2t_from_components <- function(comp,
                                weeks = c(3, 6, 9, 12, 15, 17, 21, 24)) {
  v_int  <- as.numeric(comp$v_int[1])
  v_slo  <- as.numeric(comp$v_slo[1])
  v_cov  <- as.numeric(comp$v_cov[1])
  v_litt <- as.numeric(comp$v_litt[1]); if (is.na(v_litt)) v_litt <- 0
  v_res  <- as.numeric(comp$v_res[1])
  do.call(rbind, lapply(weeks, function(w) {
    t <- w - schedule_mean
    s2A <- max(v_int + 2 * v_cov * t + v_slo * t * t, 0)
    data.frame(week = w, s2_A = s2A, s2_c = v_litt, s2_e = v_res,
               h2 = s2A / (s2A + v_litt + v_res))
  }))
}

# Delta-method SE on σ²_A(t) and h²(t) at scheduled weeks, treating the three
# additive-genetic random-regression components (v_int, v_cov, v_slope) and the
# environmental components (v_litter, v_residual) as approximately independent.
# This understates SE relative to the full AI-matrix delta method (which would
# pick up the typically negative correlation between v_int and v_slope) but
# avoids the ASReml-R `vpredict()` V-index ambiguity for str() composed random
# terms, and yields ribbons that bracket reasonable uncertainty for the figure.
#   σ²_A(t) = v_int + 2 t v_cov + t² v_slope
#   Var(σ²_A(t)) ≈ SE²(v_int) + (2t)² SE²(v_cov) + (t²)² SE²(v_slope)
#   σ²_P(t)  = σ²_A(t) + v_litter + v_residual  (held at point estimate in denom)
#   SE(h²(t)) ≈ SE(σ²_A(t)) / σ²_P(t)
h2t_se_per_week <- function(comp, weeks = c(3, 6, 9, 12, 15, 17, 21, 24)) {
  v_int  <- as.numeric(comp$v_int[1]);  se_int  <- as.numeric(comp$v_int[2])
  v_slo  <- as.numeric(comp$v_slo[1]);  se_slo  <- as.numeric(comp$v_slo[2])
  v_cov  <- as.numeric(comp$v_cov[1]);  se_cov  <- as.numeric(comp$v_cov[2])
  v_litt <- as.numeric(comp$v_litt[1]); if (is.na(v_litt)) v_litt <- 0
  v_res  <- as.numeric(comp$v_res[1])
  do.call(rbind, lapply(weeks, function(w) {
    t   <- w - schedule_mean
    s2A <- v_int + 2 * v_cov * t + v_slo * t * t
    s2P <- max(s2A, 0) + v_litt + v_res
    var_s2A <- se_int^2 + (2 * t)^2 * se_cov^2 + (t * t)^2 * se_slo^2
    se_s2A  <- sqrt(max(var_s2A, 0))
    h2     <- pmax(0, s2A) / s2P
    h2_se  <- se_s2A / s2P
    data.frame(week = w, h2 = h2, h2_se = h2_se)
  }))
}

traits <- c("bc_to_adult", "shannon", "maz")
trait_label <- c(bc_to_adult = "BC-to-adult", shannon = "Shannon", maz = "MAZ")

all_h2t <- data.frame()
all_vc  <- data.frame()

for (tr in traits) {
  cat(sprintf("\n=== RRM fit: %s ===\n", trait_label[tr]))
  fit <- tryCatch(fit_rrm(tr),
                  error = function(e) { cat("  failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) { cat("  skipping\n"); next }
  comp <- extract_components(fit)
  cat(sprintf("  v_int   = %.4g  (SE %.4g)\n", comp$v_int[1], comp$v_int[2]))
  cat(sprintf("  v_slope = %.4g  (SE %.4g)\n", comp$v_slo[1], comp$v_slo[2]))
  cat(sprintf("  v_cov   = %.4g  (SE %.4g)\n", comp$v_cov[1], comp$v_cov[2]))
  cat(sprintf("  v_litter= %.4g  (SE %.4g)\n", comp$v_litt[1], comp$v_litt[2]))
  cat(sprintf("  v_res   = %.4g  (SE %.4g)\n", comp$v_res[1], comp$v_res[2]))
  h2t <- h2t_from_components(comp); h2t$trait <- trait_label[tr]
  # Delta-method SE on h²(t) per week, components-independent approximation
  h2t_se <- h2t_se_per_week(comp)
  h2t$h2_se <- h2t_se$h2_se
  cat("  h²(t) ± SE per week:\n")
  for (i in seq_len(nrow(h2t))) {
    cat(sprintf("    W%-2d  h²=%6.3f  SE=%6.3f\n",
                h2t$week[i], h2t$h2[i], h2t$h2_se[i]))
  }
  all_h2t <- rbind(all_h2t, h2t[, c("trait", "week", "s2_A", "s2_c", "s2_e", "h2", "h2_se")])
  all_vc <- rbind(all_vc, data.frame(
    trait = trait_label[tr],
    component = c("v_int", "v_slope", "v_cov", "v_litter", "v_residual"),
    estimate  = c(comp$v_int[1], comp$v_slo[1], comp$v_cov[1], comp$v_litt[1], comp$v_res[1]),
    se        = c(comp$v_int[2], comp$v_slo[2], comp$v_cov[2], comp$v_litt[2], comp$v_res[2]),
    stringsAsFactors = FALSE))
}

write.csv(all_h2t, "results/metrics/maturation_rrm_h2_by_week.csv",
          row.names = FALSE)
write.csv(all_vc,  "results/metrics/maturation_rrm_varcomp.csv",
          row.names = FALSE)

cat("\n=== Wrote ===\n")
cat("  results/metrics/maturation_rrm_h2_by_week.csv\n")
cat("  results/metrics/maturation_rrm_varcomp.csv\n\n")
cat("=== h²(t) ===\n"); print(all_h2t, row.names = FALSE)
cat("\n=== VC + SE ===\n"); print(all_vc, row.names = FALSE)
