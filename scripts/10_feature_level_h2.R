#!/usr/bin/env Rscript
#________________________________________
# - File Name : 10_feature_level_h2.R
# - Purpose : Per-MAG, per-KO and per-KEGG-module heritability panels
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
suppressPackageStartupMessages({
  library(arrow); library(asreml); library(dplyr)
})

# ---- Pedigree / A-matrix ----
ped  <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
# Standardize column name to 'id' if needed
if (!"id" %in% names(ped3) && "animal_id" %in% names(ped3)) ped3$id <- ped3$animal_id
Ainv_path <- "data/processed/A_inverse_3gen.rds"
if (!file.exists(Ainv_path)) {
  A <- readRDS("data/processed/A_matrix_animals.rds")
  Ainv <- solve(A)
  attr(Ainv, "rowNames") <- rownames(A); attr(Ainv, "colNames") <- colnames(A)
} else {
  Ainv <- readRDS(Ainv_path)
}
litter_lookup <- setNames(as.character(ped$litter_id), ped$animal_id)
animal_ids <- as.character(ped$animal_id)

# Reusable: fit per-animal scalar through the canonical animal model
fit_per_animal_h2 <- function(scalar_per_animal) {
  if (sum(!is.na(scalar_per_animal)) < 30) return(NULL)
  if (sd(scalar_per_animal, na.rm = TRUE) < 1e-9) return(NULL)
  df <- data.frame(
    y = as.numeric(scale(scalar_per_animal)),
    animal_id = factor(animal_ids, levels = ped3$id),
    litter_id = factor(as.character(litter_lookup[animal_ids]))
  )
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1,
           random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units),
           data = df, trace = FALSE, maxit = 50,
           na.action = na.method(y = "include"))
  )), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  vc <- summary(fit)$varcomp
  rn <- rownames(vc)
  pick <- function(pat) {
    idx <- grep(pat, rn)
    if (length(idx)) max(vc[idx[1], "component"], 0) else 0
  }
  s2A <- pick("vm\\(animal_id")
  s2c <- pick("^litter_id")
  s2e <- pick("units")
  total <- s2A + s2c + s2e
  h2_vp <- tryCatch(suppressWarnings(vpredict(fit, h2 ~ V2/(V1+V2+V3))),
                    error = function(e) NULL)
  list(sigma2_a = s2A, sigma2_litter = s2c, sigma2_e = s2e,
       sigma2_total = total,
       h2 = if (total > 0) s2A / total else NA_real_,
       h2_se = if (!is.null(h2_vp)) h2_vp$SE[1] else NA_real_,
       c2 = if (total > 0) s2c / total else NA_real_,
       converged = TRUE)
}

# Helper: load a wide matrix (animals × feature_tp) and return per-animal mean
# of log1p(feature) across timepoints, for each feature.
build_per_animal_means <- function(wide_path) {
  m <- as.matrix(read_parquet(wide_path) %>% as.data.frame())
  # rownames may be in the parquet — re-read with arrow to preserve
  raw <- read_parquet(wide_path)
  raw <- as.data.frame(raw)
  if ("animal_id" %in% names(raw)) {
    rn <- as.character(raw$animal_id)
    raw$animal_id <- NULL
  } else {
    rn <- rownames(raw)
  }
  feat_cols <- colnames(raw)
  # Each column is feature_tp; split at the last '_'
  split_col <- regmatches(feat_cols, regexpr("_W[0-9]+$", feat_cols))
  if (length(split_col) != length(feat_cols)) {
    stop(sprintf("Could not parse feature_tp from columns; sample col '%s'", feat_cols[1]))
  }
  tp <- sub("^_", "", split_col)
  feature <- substr(feat_cols, 1, nchar(feat_cols) - nchar(split_col))
  unique_features <- unique(feature)
  cat(sprintf("  features: %d, timepoints: %s\n",
              length(unique_features),
              paste(sort(unique(tp)), collapse = ", ")))
  per_animal_means <- matrix(NA_real_, nrow = nrow(raw),
                             ncol = length(unique_features),
                             dimnames = list(rn, unique_features))
  for (f in unique_features) {
    cols <- which(feature == f)
    sub <- as.matrix(raw[, cols, drop = FALSE])
    sub <- log1p(pmax(sub, 0))
    per_animal_means[, f] <- rowMeans(sub, na.rm = TRUE)
  }
  per_animal_means
}

# Per-MAG
cat("\n=== Per-MAG h² ===\n")
mag_means <- build_per_animal_means("data/processed/mag_matrix_wide.parquet")
cat(sprintf("  matrix: %d animals × %d MAGs\n", nrow(mag_means), ncol(mag_means)))

# Reorder rows to match pedigree (canonical animal_ids)
mag_means <- mag_means[match(animal_ids, rownames(mag_means)), , drop = FALSE]

mag_results <- data.frame()
n_features <- ncol(mag_means)
for (i in seq_len(n_features)) {
  f <- colnames(mag_means)[i]
  out <- fit_per_animal_h2(mag_means[, f])
  if (is.null(out)) next
  mag_results <- rbind(mag_results, data.frame(
    feature = f,
    sigma2_a = out$sigma2_a, sigma2_litter = out$sigma2_litter,
    sigma2_e = out$sigma2_e, sigma2_total = out$sigma2_total,
    h2 = out$h2, h2_se = out$h2_se, c2 = out$c2, converged = out$converged,
    stringsAsFactors = FALSE))
  if (i %% 50 == 0) cat(sprintf("  ... %d / %d MAGs fit\n", i, n_features))
}
write.csv(mag_results, "results/metrics/mag_per_feature_h2.csv", row.names = FALSE)
cat(sprintf("\nWrote results/metrics/mag_per_feature_h2.csv  (%d features)\n",
            nrow(mag_results)))
cat(sprintf("  median h² = %.4f, frac h² > 0.2 = %.1f%%\n",
            median(mag_results$h2, na.rm = TRUE),
            100 * mean(mag_results$h2 > 0.2, na.rm = TRUE)))

# Per-KO
cat("\n=== Per-KO h² ===\n")
ko_means <- build_per_animal_means("data/processed/ko_matrix_wide.parquet")
cat(sprintf("  matrix: %d animals × %d KOs\n", nrow(ko_means), ncol(ko_means)))
ko_means <- ko_means[match(animal_ids, rownames(ko_means)), , drop = FALSE]

ko_results <- data.frame()
n_features <- ncol(ko_means)
for (i in seq_len(n_features)) {
  f <- colnames(ko_means)[i]
  out <- fit_per_animal_h2(ko_means[, f])
  if (is.null(out)) next
  ko_results <- rbind(ko_results, data.frame(
    feature = f,
    sigma2_a = out$sigma2_a, sigma2_litter = out$sigma2_litter,
    sigma2_e = out$sigma2_e, sigma2_total = out$sigma2_total,
    h2 = out$h2, h2_se = out$h2_se, c2 = out$c2, converged = out$converged,
    stringsAsFactors = FALSE))
  if (i %% 200 == 0) cat(sprintf("  ... %d / %d KOs fit\n", i, n_features))
}
write.csv(ko_results, "results/metrics/ko_per_feature_h2.csv", row.names = FALSE)
cat(sprintf("\nWrote results/metrics/ko_per_feature_h2.csv  (%d features)\n",
            nrow(ko_results)))
cat(sprintf("  median h² = %.4f, frac h² > 0.2 = %.1f%%\n",
            median(ko_results$h2, na.rm = TRUE),
            100 * mean(ko_results$h2 > 0.2, na.rm = TRUE)))

cat("\nDone.\n")

# PART 2: per-KEGG-module heritability (KO → module via KEGGREST)
setwd(here::here())
suppressPackageStartupMessages({
  library(arrow)
  library(asreml)
  library(KEGGREST)
})

cache_path <- "data/processed/ko_module_map.csv"

# --- KO → module mapping (cached) -------------------------------------------

if (file.exists(cache_path)) {
  cat("Loading cached KO → module mapping from", cache_path, "\n")
  ko_module <- read.csv(cache_path, stringsAsFactors = FALSE)
} else {
  cat("Fetching KO → module mapping from KEGG REST...\n")
  ko_module_raw <- tryCatch(KEGGREST::keggLink("module", "ko"), error = function(e) NULL)
  if (is.null(ko_module_raw) || length(ko_module_raw) == 0) {
    stop("KEGGREST keggLink failed; check internet access")
  }
  # ko_module_raw: named char; names like "ko:K00001", values like "md:M00001"
  ko_module <- data.frame(
    ko     = sub("^ko:", "", names(ko_module_raw)),
    module = sub("^md:", "", as.character(ko_module_raw)),
    stringsAsFactors = FALSE
  )
  write.csv(ko_module, cache_path, row.names = FALSE)
  cat(sprintf("Cached %d KO ↔ module pairs to %s\n",
              nrow(ko_module), cache_path))
}

n_unique_modules <- length(unique(ko_module$module))
n_unique_kos     <- length(unique(ko_module$ko))
cat(sprintf("Mapping: %d unique KOs ↔ %d unique modules\n",
            n_unique_kos, n_unique_modules))

# --- Load KO matrix ---------------------------------------------------------

wide <- read_parquet("data/processed/ko_matrix_wide.parquet")
animal_ids <- as.character(wide$animal_id)
cols <- setdiff(names(wide), "animal_id")
mat_wide <- as.matrix(wide[, cols])
rownames(mat_wide) <- animal_ids

parts <- strsplit(cols, "_(?=W[0-9]+$)", perl = TRUE)
ko_per_col <- sapply(parts, function(x) x[1])
tp_per_col <- sapply(parts, function(x) x[2])
unique_kos <- unique(ko_per_col)
unique_tps <- sort(unique(tp_per_col))
cat(sprintf("KO matrix: %d animals × %d (KO × tp) columns; %d KOs × %d tps\n",
            length(animal_ids), length(cols),
            length(unique_kos), length(unique_tps)))

# --- KO → module aggregation ------------------------------------------------

# Restrict to KOs that map to at least one module
ko_in_map <- intersect(unique_kos, ko_module$ko)
cat(sprintf("KOs in module map: %d / %d (%.1f%%)\n",
            length(ko_in_map), length(unique_kos),
            100 * length(ko_in_map) / length(unique_kos)))

# Build module list (only modules with KOs in our data)
modules_to_use <- unique(ko_module$module[ko_module$ko %in% ko_in_map])
cat(sprintf("Modules in our data: %d\n", length(modules_to_use)))

# Build per-animal × per-tp module-level matrix by summing KO counts per module
# Strategy: for each module × tp, sum the (KO × tp) columns over KOs in module
agg_mat <- matrix(0, nrow = length(animal_ids),
                  ncol = length(modules_to_use) * length(unique_tps))
agg_cols <- character(0)
col_idx <- 1
for (m in modules_to_use) {
  kos_in_m <- ko_module$ko[ko_module$module == m]
  for (tp in unique_tps) {
    src_cols <- which(ko_per_col %in% kos_in_m & tp_per_col == tp)
    if (length(src_cols) == 0) {
      agg_mat[, col_idx] <- NA_real_
    } else if (length(src_cols) == 1) {
      agg_mat[, col_idx] <- mat_wide[, src_cols]
    } else {
      agg_mat[, col_idx] <- rowSums(mat_wide[, src_cols, drop = FALSE], na.rm = FALSE)
    }
    agg_cols <- c(agg_cols, paste0(m, "_", tp))
    col_idx <- col_idx + 1
  }
}
colnames(agg_mat) <- agg_cols
rownames(agg_mat) <- animal_ids
log_agg <- log1p(agg_mat)
cat(sprintf("Module matrix: %d animals × %d (module × tp)\n",
            nrow(agg_mat), ncol(agg_mat)))

# --- Per-module h² (animal model, summing across timepoints per module) ----

# Pedigree
ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
ped3 <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
Ainv <- readRDS("data/processed/A_inverse_3gen.rds")
litter_lookup <- setNames(as.character(ped$litter_id), ped$animal_id)

extract <- function(fit) {
  vc <- summary(fit)$varcomp
  rn <- rownames(vc)
  s2_A <- {a <- grep("vm\\(animal_id", rn); if (length(a)) max(vc[a[1], "component"], 0) else 0}
  s2_c <- {c <- which(rn == "litter_id");   if (length(c)) max(vc[c, "component"], 0) else 0}
  s2_e <- {e <- grep("units", rn);          if (length(e)) vc[e[1], "component"] else NA}
  total <- s2_A + s2_c + s2_e
  list(
    s2_A = s2_A, s2_c = s2_c, s2_e = s2_e, total = total,
    h2 = if (total > 0) s2_A / total else NA,
    c2 = if (total > 0) s2_c / total else NA
  )
}

# Strategy: per module, fit ONE animal model on the per-animal mean of log(module count) across
# timepoints (analogous to per-KO scalar). This gives one h² per module, comparable to per-KO h².

cat("\n=== Per-module h² (mean log-count across timepoints) ===\n")
results <- data.frame()

for (m in modules_to_use) {
  # collect this module's columns
  m_cols <- grep(paste0("^", m, "_W"), agg_cols, value = TRUE)
  if (length(m_cols) == 0) next
  vals_per_animal <- rowMeans(log_agg[, m_cols, drop = FALSE], na.rm = TRUE)
  if (sum(!is.na(vals_per_animal)) < 30) next
  if (sd(vals_per_animal, na.rm = TRUE) < 1e-6) next

  df_fit <- data.frame(
    y = as.numeric(scale(vals_per_animal)),
    animal_id = factor(animal_ids, levels = ped3$id),
    litter_id = factor(as.character(litter_lookup[animal_ids]))
  )
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(y ~ 1, random = ~ litter_id + vm(animal_id, Ainv),
           residual = ~ idv(units), data = df_fit, trace = FALSE, maxit = 50,
           na.action = na.method(y = "include"))
  )), error = function(e) NULL)
  if (is.null(fit)) next
  p <- extract(fit)
  h2_vp <- tryCatch(suppressWarnings(vpredict(fit, h2 ~ V2/(V1+V2+V3))),
                    error = function(e) NULL)
  h2_se <- if (!is.null(h2_vp)) h2_vp$SE[1] else NA_real_
  results <- rbind(results, data.frame(
    module = m, n_animals = sum(!is.na(vals_per_animal)),
    n_kos_in_module = sum(ko_module$module == m & ko_module$ko %in% ko_in_map),
    s2_A = p$s2_A, s2_c = p$s2_c, s2_e = p$s2_e,
    h2 = p$h2, h2_se = h2_se, c2 = p$c2,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("Computed h² for %d modules\n", nrow(results)))
cat(sprintf("Mean module h²:   %.4f\n", mean(results$h2, na.rm = TRUE)))
cat(sprintf("Median module h²: %.4f\n", median(results$h2, na.rm = TRUE)))
cat(sprintf("Mean module c²:   %.4f\n", mean(results$c2, na.rm = TRUE)))
cat(sprintf("Comparison reference (paper §S4): per-KO mean h² = 0.028\n"))

write.csv(results, "results/metrics/ko_per_module_h2.csv", row.names = FALSE)
cat("\n=== output ===\nresults/metrics/ko_per_module_h2.csv\n")
