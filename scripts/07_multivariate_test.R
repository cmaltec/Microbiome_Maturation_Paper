#!/usr/bin/env Rscript
#________________________________________
# - File Name : 07_multivariate_test.R
# - Purpose : Hotelling's T2 and PERMANOVA on the six-class metric vector
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
  library(vegan)
})

# --- Load per-animal metrics (per-metric source rule per §S2) ---------------

obs <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                stringsAsFactors = FALSE)
imp <- read.csv("results/metrics/maturation_all_metrics_per_animal_imputed.csv",
                stringsAsFactors = FALSE)
thr <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv",
                stringsAsFactors = FALSE)
ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

# Per-class representative metric (highest-h² from each class; per-metric
# source rule per Supp §S2: stability=obs, distance=imp, threshold=imp,
# others=obs).
metric_panel <- list(
  list(name = "maz_slope",          class = "1-Age",         source = "obs", col = "maz_slope"),
  list(name = "aitchison_slope",    class = "2-Distance",    source = "imp", col = "aitchison_to_adult_slope"),
  list(name = "t75",                class = "3-Growth-AltC", source = "imp", col = "t75"),
  list(name = "intra_bc_var",       class = "4-Stability",   source = "obs", col = "intra_bc_var"),
  list(name = "shannon_slope",      class = "5-Diversity",   source = "obs", col = "shannon_slope"),
  list(name = "BC_FPC2",            class = "6-FPCA",        source = "obs", col = "BC_FPC2"))

# Build the 6-vector matrix per animal using per-metric source -------------
build_vector <- function() {
  out <- data.frame(animal_id = obs$animal_id)
  for (m in metric_panel) {
    src <- if (m$source == "obs") obs else if (m$source == "imp") imp else thr
    if (m$col %in% names(src)) {
      idx <- match(out$animal_id, src$animal_id)
      out[[m$name]] <- src[[m$col]][idx]
    } else if (m$col %in% names(thr) && m$source == "imp") {
      idx <- match(out$animal_id, thr$animal_id)
      out[[m$name]] <- thr[[m$col]][idx]
    } else {
      stop(sprintf("col '%s' not found in %s source", m$col, m$source))
    }
  }
  out
}
mat_df <- build_vector()
mat_df <- merge(mat_df, ped[, c("animal_id", "line", "sire_id", "litter_id")],
                by = "animal_id", all.x = TRUE)

# Drop animals with any NA in the 6-vector
metric_cols <- sapply(metric_panel, function(m) m$name)
keep <- complete.cases(mat_df[, metric_cols])
mat_df <- mat_df[keep, ]
cat(sprintf("Animals with complete 6-vector: %d\n", nrow(mat_df)))
cat(sprintf("Per-class composition (line):\n"))
print(table(mat_df$line))

# Z-score columns ------------------------------------------------------------
Y <- as.matrix(mat_df[, metric_cols])
Y_z <- scale(Y)

# === Test 1: Hotelling's T^2 — SEN vs TOL ==================================
hotelling_T2 <- function(Y_z, group) {
  groups <- unique(group)
  if (length(groups) != 2) stop("need exactly 2 groups for Hotelling's T^2")
  g1 <- group == groups[1]; g2 <- group == groups[2]
  n1 <- sum(g1); n2 <- sum(g2)
  m1 <- colMeans(Y_z[g1, , drop = FALSE])
  m2 <- colMeans(Y_z[g2, , drop = FALSE])
  S1 <- cov(Y_z[g1, , drop = FALSE])
  S2 <- cov(Y_z[g2, , drop = FALSE])
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  d  <- m1 - m2
  T2 <- (n1 * n2 / (n1 + n2)) * as.numeric(t(d) %*% MASS::ginv(Sp) %*% d)
  p <- ncol(Y_z)
  F_stat <- T2 * (n1 + n2 - p - 1) / (p * (n1 + n2 - 2))
  pval <- pf(F_stat, df1 = p, df2 = n1 + n2 - p - 1, lower.tail = FALSE)
  list(T2 = T2, F = F_stat, df1 = p, df2 = n1 + n2 - p - 1, pval = pval,
       n1 = n1, n2 = n2, group_labels = groups)
}

cat("\n=== Hotelling's T^2: SEN vs TOL on 6-class metric vector ===\n")
ht_line <- hotelling_T2(Y_z, mat_df$line)
cat(sprintf("  T^2 = %.3f, F(%d, %d) = %.3f, p = %.4g\n",
            ht_line$T2, ht_line$df1, ht_line$df2, ht_line$F, ht_line$pval))
cat(sprintf("  groups: %s (n=%d) vs %s (n=%d)\n",
            ht_line$group_labels[1], ht_line$n1,
            ht_line$group_labels[2], ht_line$n2))

# Mahalanobis-distance permutation alternative (more robust at n=78) -------
mahal_perm <- function(Y_z, group, n_perm = 9999) {
  groups <- unique(group)
  m1 <- colMeans(Y_z[group == groups[1], , drop = FALSE])
  m2 <- colMeans(Y_z[group == groups[2], , drop = FALSE])
  Sp <- cov(Y_z)
  d_obs <- as.numeric(sqrt(t(m1 - m2) %*% MASS::ginv(Sp) %*% (m1 - m2)))
  perm_d <- replicate(n_perm, {
    perm_g <- sample(group)
    pm1 <- colMeans(Y_z[perm_g == groups[1], , drop = FALSE])
    pm2 <- colMeans(Y_z[perm_g == groups[2], , drop = FALSE])
    as.numeric(sqrt(t(pm1 - pm2) %*% MASS::ginv(Sp) %*% (pm1 - pm2)))
  })
  list(d_mahal = d_obs, p_perm = mean(perm_d >= d_obs))
}
set.seed(42)
mp <- mahal_perm(Y_z, mat_df$line)
cat(sprintf("  Mahalanobis distance: %.3f, permutation p = %.4g (9999 perms)\n",
            mp$d_mahal, mp$p_perm))

# === Test 2: PERMANOVA — variance by line, sire, litter ====================
cat("\n=== PERMANOVA on 6-class metric vector ===\n")
mat_df$sire_id_f <- factor(mat_df$sire_id)
mat_df$litter_id_f <- factor(mat_df$litter_id)
mat_df$line_f <- factor(mat_df$line)

D <- dist(Y_z, method = "euclidean")
# Line is fully nested within sire (each sire is in one line) — so we test
# them in separate models rather than jointly with marginal effects.
adonis_line <- adonis2(D ~ line_f, data = mat_df, permutations = 9999)
cat("PERMANOVA — line alone:\n")
print(adonis_line)
adonis_sire <- adonis2(D ~ sire_id_f, data = mat_df, permutations = 9999)
cat("\nPERMANOVA — sire alone:\n")
print(adonis_sire)
adonis_litter <- adonis2(D ~ litter_id_f, data = mat_df, permutations = 9999)
cat("\nPERMANOVA — litter alone:\n")
print(adonis_litter)

# === Save outputs ===========================================================
ht_out <- data.frame(
  test = c("Hotelling_T2_line", "Mahalanobis_perm_line"),
  statistic = c(ht_line$T2, mp$d_mahal),
  F_or_d = c(ht_line$F, mp$d_mahal),
  df1 = c(ht_line$df1, NA),
  df2 = c(ht_line$df2, NA),
  pval = c(ht_line$pval, mp$p_perm),
  n_total = c(ht_line$n1 + ht_line$n2, nrow(mat_df)),
  metrics_in_vector = paste(metric_cols, collapse = ";")
)
write.csv(ht_out,
          "results/metrics/maturation_multivariate_tests.csv",
          row.names = FALSE)

stack_adonis <- function(a, label) {
  data.frame(
    model      = label,
    term       = rownames(a),
    Df         = a$Df,
    SumOfSqs   = a$SumOfSqs,
    R2         = a$R2,
    F          = a$F,
    pval       = a$`Pr(>F)`,
    stringsAsFactors = FALSE
  )
}
permanova_out <- rbind(
  stack_adonis(adonis_line,   "line_only"),
  stack_adonis(adonis_sire,   "sire_only"),
  stack_adonis(adonis_litter, "litter_only")
)
write.csv(permanova_out,
          "results/metrics/maturation_multivariate_permanova.csv",
          row.names = FALSE)

cat("\n=== outputs written ===\n")
cat("  results/metrics/maturation_multivariate_tests.csv\n")
cat("  results/metrics/maturation_multivariate_permanova.csv\n")
