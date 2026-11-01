#!/usr/bin/env Rscript
#________________________________________
# - File Name : 00c_revision_inputs.R
# - Purpose : 3-gen and 2-gen A-inverse, per-sample long table and per-animal traits for scripts 06, 13, 14
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
# Sun Nov 01 15:30:00 2026 ------------------------------

setwd(here::here())

suppressPackageStartupMessages({ library(asreml); library(dplyr) })

OUT <- "results/revision/inputs"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 3-generation pedigree (same construction as 04_pedigree.R) -------------
d <- read.csv("data/pedigree_3generation.csv", stringsAsFactors = FALSE)
names(d)[names(d) == "litter_id"] <- "dam_id"
gp <- distinct(d[, c("dam_id", "grand_sire", "grand_dam")])
stopifnot(!anyDuplicated(gp$dam_id))

ped3 <- bind_rows(
  data.frame(id = unique(na.omit(d$grand_sire)), sire = NA_character_, dam = NA_character_),
  data.frame(id = unique(na.omit(as.character(d$grand_dam))), sire = NA_character_, dam = NA_character_),
  data.frame(id = unique(d$sire_id), sire = NA_character_, dam = NA_character_),
  data.frame(id = gp$dam_id, sire = gp$grand_sire, dam = as.character(gp$grand_dam)),
  data.frame(id = d$animal_id, sire = d$sire_id, dam = d$dam_id)
) %>% distinct(id, .keep_all = TRUE)
slash <- ped3$id[grepl("/", ped3$id)]
ped3$sire[ped3$sire %in% slash] <- NA
ped3$dam[ped3$dam %in% slash]   <- NA
ped3 <- ped3[!ped3$id %in% slash, ]
cat(sprintf("3-gen pedigree: %d ids (dropped ambiguous: %s)\n",
            nrow(ped3), paste(slash, collapse = ", ")))
stopifnot(nrow(ped3) == 114)
write.csv(ped3, file.path(OUT, "pedigree_3gen_3col.csv"), row.names = FALSE)

p <- ped3; p$sire[is.na(p$sire)] <- 0; p$dam[is.na(p$dam)] <- 0
Ainv3 <- ainverse(p)
saveRDS(Ainv3, file.path(OUT, "Ainv_3gen.rds"))

# ---- 2-generation pedigree (as 00b_build_a_matrix.R) ------------------------
ped2 <- rbind(
  data.frame(id = c(unique(d$sire_id), paste0("DAM_", unique(d$dam_id))),
             sire = 0, dam = 0),
  data.frame(id = d$animal_id, sire = d$sire_id, dam = paste0("DAM_", d$dam_id)))
Ainv2 <- ainverse(ped2)
saveRDS(Ainv2, file.path(OUT, "Ainv_2gen.rds"))

# ---- animal metadata --------------------------------------------------------
meta <- read.csv("data/sample_metadata.csv", stringsAsFactors = FALSE)
anim <- distinct(meta[, c("animal_id", "litter_id", "sire_id", "line")])
stopifnot(nrow(anim) == 78, !anyDuplicated(anim$animal_id))
ld <- distinct(data.frame(litter = anim$litter_id,
                          dam = d$dam_id[match(anim$animal_id, d$animal_id)]))
stopifnot(all(ld$litter == ld$dam))
write.csv(anim, file.path(OUT, "animals.csv"), row.names = FALSE)

# ---- long per-sample table for random regression ---------------------------
cur <- read.csv("results/metrics/maturation_curves_per_sample.csv",
                stringsAsFactors = FALSE)
maz <- read.csv("results/metrics/maz_per_sample.csv", stringsAsFactors = FALSE)
long <- cur[, c("animal_id", "week", "tp", "bc_to_adult", "shannon")] %>%
  left_join(maz[, c("animal_id", "tp", "maz")], by = c("animal_id", "tp")) %>%
  left_join(anim, by = "animal_id") %>%
  arrange(animal_id, week)
cat(sprintf("Long table: %d records, %d animals; per week:\n",
            nrow(long), length(unique(long$animal_id))))
print(table(long$week))
stopifnot(nrow(long) == 500)
write.csv(long, file.path(OUT, "per_sample_long.csv"), row.names = FALSE)

# ---- per-animal traits (observed-data scalars used by 05) -------------------
met <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv",
                stringsAsFactors = FALSE)
met <- inner_join(met, anim[, c("animal_id", "litter_id", "sire_id")], by = "animal_id")
write.csv(met, file.path(OUT, "per_animal_traits.csv"), row.names = FALSE)
cat(sprintf("Per-animal traits: %d animals\n", nrow(met)))
