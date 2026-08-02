#!/usr/bin/env Rscript
#________________________________________
# - File Name : 00b_build_a_matrix.R
# - Purpose : Two-generation pedigree A matrix and ASReml-form A-inverse
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
# Sat Aug 01 13:44:36 2026 ------------------------------

# Dam is unknown, so litter_id is the dam proxy: animals sharing a litter share
# a dam. Same-litter animals are therefore full sibs (r = 0.5) and
# same-sire-across-litters animals paternal half sibs (r = 0.25).

setwd(here::here())
suppressPackageStartupMessages({
  library(Matrix); library(asreml)
})

ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
stopifnot(all(c("animal_id", "litter_id", "sire_id") %in% names(ped)))

litter_sire <- unique(ped[, c("litter_id", "sire_id")])
ls_counts <- table(litter_sire$litter_id)
if (any(ls_counts > 1)) {
  stop("litter_id mapped to multiple sires, dam-proxy invalid: ",
       paste(names(ls_counts)[ls_counts > 1], collapse = ", "))
}

animal_ids <- ped$animal_id
sire_ids   <- as.character(unique(ped$sire_id))
dam_ids    <- paste0("DAM_", unique(ped$litter_id))
stopifnot(length(intersect(animal_ids, sire_ids)) == 0)
stopifnot(length(intersect(animal_ids, dam_ids))  == 0)
stopifnot(length(intersect(sire_ids,   dam_ids))  == 0)

ped2 <- rbind(
  data.frame(id = c(sire_ids, dam_ids), sire = NA_character_,
             dam = NA_character_, stringsAsFactors = FALSE),
  data.frame(id = ped$animal_id, sire = as.character(ped$sire_id),
             dam = paste0("DAM_", ped$litter_id), stringsAsFactors = FALSE)
)
n   <- nrow(ped2)
idx <- setNames(seq_len(n), ped2$id)

# Henderson tabular algorithm. Parents precede offspring.
A <- matrix(0, nrow = n, ncol = n)
for (i in seq_len(n)) {
  s <- ped2$sire[i]; d <- ped2$dam[i]
  si <- if (!is.na(s)) idx[[s]] else NA_integer_
  di <- if (!is.na(d)) idx[[d]] else NA_integer_

  if (is.na(si) && is.na(di)) {
    A[i, i] <- 1
  } else if (is.na(si) || is.na(di)) {
    known <- if (is.na(si)) di else si
    A[i, i] <- 1
    if (i > 1) {
      for (j in seq_len(i - 1)) {
        A[j, i] <- 0.5 * A[j, known]
        A[i, j] <- A[j, i]
      }
    }
  } else {
    A[i, i] <- 1 + 0.5 * A[si, di]
    if (i > 1) {
      for (j in seq_len(i - 1)) {
        A[j, i] <- 0.5 * (A[j, si] + A[j, di])
        A[i, j] <- A[j, i]
      }
    }
  }
}
rownames(A) <- ped2$id
colnames(A) <- ped2$id

A_an <- A[animal_ids, animal_ids]
A_an_sp <- as(A_an, "dgCMatrix")
saveRDS(A_an_sp, "data/processed/A_matrix_animals.rds")

# ASReml 3-column triplet form (Row, Column, Ainverse) that vm() expects.
ped_for_ai <- ped2
ped_for_ai$sire[is.na(ped_for_ai$sire)] <- 0
ped_for_ai$dam[is.na(ped_for_ai$dam)]   <- 0
Ainv <- ainverse(ped_for_ai)
saveRDS(Ainv, "data/processed/A_inverse_proper.rds")

# --- sanity checks -----------------------------------------------------------
ped_lookup <- ped[, c("animal_id", "litter_id", "sire_id")]
rownames(ped_lookup) <- ped_lookup$animal_id
pairs <- t(combn(animal_ids, 2))
rel <- character(nrow(pairs))
for (k in seq_len(nrow(pairs))) {
  a <- pairs[k, 1]; b <- pairs[k, 2]
  same_litter <- ped_lookup[a, "litter_id"] == ped_lookup[b, "litter_id"]
  same_sire   <- ped_lookup[a, "sire_id"]   == ped_lookup[b, "sire_id"]
  rel[k] <- if (same_litter && same_sire) "full_sib"
            else if (same_sire) "paternal_half_sib"
            else "unrelated"
}
obs <- A_an[lower.tri(A_an)]
by_rel <- tapply(obs, rel, function(x)
  c(n = length(x), mean = mean(x), min = min(x), max = max(x)))

sink("data/processed/A_matrix_summary.txt")
cat("A matrix build summary\n======================\n\n")
cat("n animals: ", length(animal_ids), "\n")
cat("n founders (sires + dam-proxies): ", length(sire_ids) + length(dam_ids), "\n")
cat("  sires: ", length(sire_ids), "\n")
cat("  dam-proxy litters: ", length(dam_ids), "\n\n")
cat("Diagonal (self-relationship):\n"); print(summary(diag(A_an)))
cat("\nOff-diagonal distribution:\n"); print(table(round(obs, 4)))
cat("\nBy expected relationship class:\n"); print(by_rel)
cat("\nExpected values: full_sib=0.5, paternal_half_sib=0.25, unrelated=0\n")
sink()

cat(sprintf("A (animals): %d x %d\n", nrow(A_an_sp), ncol(A_an_sp)))
cat(sprintf("Ainv: %d non-zero entries, %d ids, logdet = %.3f\n",
            nrow(Ainv), length(attr(Ainv, "rowNames")), attr(Ainv, "logdet")))
cat("Wrote: data/processed/A_matrix_animals.rds\n")
cat("Wrote: data/processed/A_inverse_proper.rds\n")
cat("Wrote: data/processed/A_matrix_summary.txt\n")
