#!/usr/bin/env Rscript
#________________________________________
# - File Name : 04_pedigree.R
# - Purpose : Three-generation pedigree and ASReml-form A and A-inverse
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
  library(dplyr); library(tidyr)
})

# Dam and litter share one identifier throughout, so litter_id serves as dam_id.
d <- read.csv("data/pedigree_3generation.csv", stringsAsFactors = FALSE)
names(d)[names(d) == "litter_id"] <- "dam_id"
cat(sprintf("Pedigree+growth table: %d animals (%d with growth data)\n",
            nrow(d), sum(!is.na(d$birthweight))))

# ---- Save growth phenotypes ------------------------------------------------
growth <- d %>% select(animal_id, line, sire_id, dam_id,
                       birthweight, d7_weight, d14_weight, weaning_weight,
                       pre_weaning_adg)
write.csv(growth, "data/metadata/growth_phenotypes.csv", row.names = FALSE)
cat(sprintf("Saved: data/metadata/growth_phenotypes.csv (%d animals)\n", nrow(growth)))

# ---- Build 3-generation pedigree -------------------------------------------
# Format: id, sire, dam (the 3-column format ASReml needs for ainverse())
# Include: 78 animals + 13 dams + 7 sires + 6 grandsires + 10 granddams = 114
# (two dual-paternity grandsire ids are dropped below, 8 unique -> 6)

# Map dam_id → grand_sire, grand_dam (each dam should have consistent grandparents)
dam_to_gparents <- d %>%
  select(dam_id, grand_sire, grand_dam) %>%
  distinct()
cat(sprintf("\nDam → grandparents map: %d dams, %d unique grand sires, %d unique grand dams\n",
            nrow(dam_to_gparents),
            length(unique(dam_to_gparents$grand_sire)),
            length(unique(dam_to_gparents$grand_dam))))

# Build the 3-column pedigree:
# Top generation: grand sires + grand dams (no parents)
# Middle:         sires + dams (parents are grand sire + grand dam where known)
# Bottom:         animals (parents are sire + dam)

# Grand parents (no parents)
grand_sires <- unique(na.omit(d$grand_sire))
grand_dams  <- unique(na.omit(d$grand_dam))

# Sires (also no parents — grandparents not provided for sires)
sires <- unique(d$sire_id)
dams  <- unique(d$dam_id)

# Build pedigree rows
ped3gen <- bind_rows(
  data.frame(id = as.character(grand_sires), sire = NA_character_, dam = NA_character_),
  data.frame(id = as.character(grand_dams),  sire = NA_character_, dam = NA_character_),
  data.frame(id = as.character(sires), sire = NA_character_, dam = NA_character_),
  data.frame(id = as.character(dams),
             sire = as.character(dam_to_gparents$grand_sire[match(dams, dam_to_gparents$dam_id)]),
             dam  = as.character(dam_to_gparents$grand_dam[ match(dams, dam_to_gparents$dam_id)])),
  data.frame(id = as.character(d$animal_id),
             sire = as.character(d$sire_id),
             dam  = as.character(d$dam_id))
) %>% distinct(id, .keep_all = TRUE)

# Drop ambiguous dual-paternity grandsire IDs (e.g. "2024571/1418638"); they encode
# uncertain pedigree links and would arbitrarily bias relationships if retained.
slash_ids <- ped3gen$id[grepl("/", ped3gen$id)]
if (length(slash_ids) > 0) {
  cat(sprintf("Dropping %d ambiguous dual-paternity ids: %s\n",
              length(slash_ids), paste(slash_ids, collapse = ", ")))
  ped3gen$sire[ped3gen$sire %in% slash_ids] <- NA
  ped3gen$dam[ped3gen$dam   %in% slash_ids] <- NA
  ped3gen <- ped3gen[!ped3gen$id %in% slash_ids, ]
}

# Verify
cat(sprintf("\n3-gen pedigree: %d total ids\n", nrow(ped3gen)))
cat(sprintf("  founders (NA parents): %d\n", sum(is.na(ped3gen$sire) & is.na(ped3gen$dam))))
cat(sprintf("  animals (both parents known): %d\n", sum(!is.na(ped3gen$sire) & !is.na(ped3gen$dam))))

# Sanity-check: every animal's sire and dam are in the id column
miss_sire <- setdiff(na.omit(ped3gen$sire), ped3gen$id)
miss_dam  <- setdiff(na.omit(ped3gen$dam),  ped3gen$id)
cat(sprintf("  sires missing from id column: %d\n", length(miss_sire)))
cat(sprintf("  dams missing from id column: %d\n", length(miss_dam)))
if (length(miss_sire) > 0) cat("    Missing sires:", paste(miss_sire, collapse=", "), "\n")
if (length(miss_dam)  > 0) cat("    Missing dams:",  paste(miss_dam,  collapse=", "), "\n")

# Save
write.csv(ped3gen, "data/metadata/pedigree_3gen_3col.csv", row.names = FALSE)
cat("\nSaved: data/metadata/pedigree_3gen_3col.csv\n")

# Also save extended pedigree with line info for animals
ped3gen_wide <- d %>% select(animal_id, line, sire_id, dam_id, grand_sire, grand_dam,
                              birthweight, d7_weight, d14_weight, weaning_weight,
                              pre_weaning_adg)
write.csv(ped3gen_wide, "data/metadata/pedigree_3gen.csv", row.names = FALSE)
cat("Saved: data/metadata/pedigree_3gen.csv\n")

# ---- Heritability and pre-weaning growth summary ----------------------------
cat(sprintf("\n=== Pre-weaning growth phenotypes (%d animals with data) ===\n",
            sum(!is.na(d$birthweight))))
for (col in c("birthweight", "d7_weight", "d14_weight", "weaning_weight", "pre_weaning_adg")) {
  v <- d[[col]]
  cat(sprintf("  %-18s n = %d, mean = %.3f, SD = %.3f, range = [%.2f, %.2f]\n",
              col, sum(!is.na(v)), mean(v, na.rm=TRUE), sd(v, na.rm=TRUE),
              min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
}

# PART 2: Build A (numerator relationship matrix) and Ainv from 3-gen pedigree
setwd(here::here())
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(asreml)
})

ped3gen <- read.csv("data/metadata/pedigree_3gen_3col.csv", stringsAsFactors = FALSE)
ped3gen$id   <- as.character(ped3gen$id)
ped3gen$sire <- as.character(ped3gen$sire)
ped3gen$dam  <- as.character(ped3gen$dam)
ped3gen$sire[ped3gen$sire %in% c("NA", "")] <- NA
ped3gen$dam[ped3gen$dam %in% c("NA", "")] <- NA

cat(sprintf("3-gen pedigree: %d ids\n", nrow(ped3gen)))

# Order pedigree so parents precede progeny.
# Founders first (no parents), then animals whose parents are already listed.
# Iterative topological order.
remaining <- ped3gen
ordered <- character(0)
while (nrow(remaining) > 0) {
  # ids whose parents are NA or already in ordered
  ready <- remaining$id[
    (is.na(remaining$sire) | remaining$sire %in% ordered) &
    (is.na(remaining$dam)  | remaining$dam  %in% ordered)
  ]
  if (length(ready) == 0) {
    cat("Topological order failed; remaining ids:\n")
    print(remaining)
    stop("Pedigree has a cycle or unresolved parents.")
  }
  ordered <- c(ordered, ready)
  remaining <- remaining[!remaining$id %in% ready, ]
}
ped3gen <- ped3gen[match(ordered, ped3gen$id), ]
cat(sprintf("Ordered: %d founders + %d animals with parents\n",
            sum(is.na(ped3gen$sire) & is.na(ped3gen$dam)),
            sum(!is.na(ped3gen$sire) | !is.na(ped3gen$dam))))

# Henderson's tabular algorithm
n <- nrow(ped3gen)
idx <- setNames(seq_len(n), ped3gen$id)
A <- matrix(0, nrow = n, ncol = n)
for (i in seq_len(n)) {
  s <- ped3gen$sire[i]; d <- ped3gen$dam[i]
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
rownames(A) <- ped3gen$id
colnames(A) <- ped3gen$id

# Animal block (just the 78 animals)
growth <- read.csv("data/metadata/growth_phenotypes.csv", stringsAsFactors = FALSE)
animal_ids <- as.character(growth$animal_id)
A_an <- A[animal_ids, animal_ids]
A_an_sp <- as(A_an, "dgCMatrix")
saveRDS(A_an_sp, "data/processed/A_matrix_3gen.rds")
cat(sprintf("\nA-matrix (animals only): %d x %d\n", nrow(A_an_sp), ncol(A_an_sp)))

# Full Ainv via ASReml's ainverse() — produces the 3-column triplet form
# (Row, Column, Ainverse) with rowNames + logdet attributes that vm() expects.
ped_for_ai <- ped3gen
ped_for_ai$sire[is.na(ped_for_ai$sire)] <- 0
ped_for_ai$dam[is.na(ped_for_ai$dam)]   <- 0
Ainv_3col <- ainverse(ped_for_ai)
saveRDS(Ainv_3col, "data/processed/A_inverse_3gen.rds")
cat(sprintf("Ainv (full pedigree, ASReml ainverse): %d non-zero entries, %d ids, logdet = %.3f\n",
            nrow(Ainv_3col), length(attr(Ainv_3col, "rowNames")),
            attr(Ainv_3col, "logdet")))

cat("\nSaved:\n")
cat("  data/processed/A_matrix_3gen.rds\n")
cat("  data/processed/A_inverse_3gen.rds\n")

# ---- Sanity checks ---------------------------------------------------------
# Full sibs (same sire + same dam → r ≥ 0.5 with possible inbreeding adjustment)
current_ped <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)
fs_pair <- current_ped %>%
  group_by(sire_id, litter_id) %>%
  filter(n() >= 2) %>%
  summarize(a1 = first(animal_id), a2 = nth(animal_id, 2), .groups = "drop") %>%
  slice(1)
if (nrow(fs_pair) > 0) {
  a1 <- fs_pair$a1[1]; a2 <- fs_pair$a2[1]
  if (a1 %in% rownames(A_an) && a2 %in% rownames(A_an)) {
    cat(sprintf("\nFull-sib check (%s × %s): A = %.4f (expected ~0.5)\n",
                a1, a2, A_an[a1, a2]))
  }
}

# Paternal half-sibs
phs_pair <- current_ped %>%
  group_by(sire_id) %>%
  filter(n_distinct(litter_id) >= 2) %>%
  do(head(., 2))
if (nrow(phs_pair) >= 2) {
  a1 <- phs_pair$animal_id[1]; a2 <- phs_pair$animal_id[2]
  if (a1 %in% rownames(A_an) && a2 %in% rownames(A_an)) {
    cat(sprintf("Paternal half-sib check (%s × %s): A = %.4f (expected ~0.25)\n",
                a1, a2, A_an[a1, a2]))
  }
}

# Diagonal range
cat(sprintf("\nDiagonal range: [%.4f, %.4f] (1.0 = no inbreeding; >1 = inbred)\n",
            min(diag(A_an)), max(diag(A_an))))

# Compare to current 2-gen A
A_old <- readRDS("data/processed/A_matrix_animals.rds")
common <- intersect(animal_ids, rownames(A_old))
cat(sprintf("\nComparison to 2-gen A on %d common animals:\n", length(common)))
A_old_sub <- as.matrix(A_old[common, common])
A_new_sub <- as.matrix(A_an[common, common])

# Off-diagonal
off_old <- A_old_sub[upper.tri(A_old_sub)]
off_new <- A_new_sub[upper.tri(A_new_sub)]
delta <- off_new - off_old
cat(sprintf("  Off-diagonal entries: %d\n", length(delta)))
cat(sprintf("  Mean change: %.5f, Max abs change: %.4f\n", mean(delta), max(abs(delta))))
cat(sprintf("  Entries changed > 0.01: %d\n", sum(abs(delta) > 0.01)))
cat(sprintf("  Entries changed > 0.05: %d\n", sum(abs(delta) > 0.05)))

# Show a few examples of changed pairs
if (sum(abs(delta) > 0.01) > 0) {
  cat("\nTop 5 most-changed off-diagonal pairs:\n")
  changed_idx <- which(upper.tri(A_old_sub) & abs(A_new_sub - A_old_sub) > 0.01,
                       arr.ind = TRUE)
  if (nrow(changed_idx) > 0) {
    top5 <- changed_idx[order(-abs(A_new_sub[changed_idx] - A_old_sub[changed_idx]))[1:min(5, nrow(changed_idx))], ,
                        drop = FALSE]
    for (k in seq_len(nrow(top5))) {
      i <- top5[k, 1]; j <- top5[k, 2]
      a1 <- common[i]; a2 <- common[j]
      cat(sprintf("  %s × %s: old=%.4f → new=%.4f (Δ=%.4f)\n",
                  a1, a2, A_old_sub[i, j], A_new_sub[i, j],
                  A_new_sub[i, j] - A_old_sub[i, j]))
    }
  }
}
