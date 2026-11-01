#!/usr/bin/env Rscript
#________________________________________
# - File Name : 15_sequencing_depth.R
# - Purpose : Per-sample sequencing depth after QC and reads mapped to the 633 core MAGs (Methods)
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

OUT <- "results/revision/sequencing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
lg <- read.csv(gzfile("data/mag_abundance_long.csv.gz"))
dep <- aggregate(count ~ sample_id + animal_id + timepoint, data = lg, FUN = sum)
stopifnot(nrow(dep) == 500)
write.csv(dep, file.path(OUT, "mapped_reads_per_sample.csv"), row.names = FALSE)
s <- data.frame(n = nrow(dep), mean = mean(dep$count), sd = sd(dep$count),
                median = median(dep$count), min = min(dep$count), max = max(dep$count))
write.csv(s, file.path(OUT, "mapped_reads_summary.csv"), row.names = FALSE)
print(s)

rc   <- read.delim("data/atlas_read_counts.tsv", stringsAsFactors = FALSE)
meta <- read.csv("data/sample_metadata.csv", stringsAsFactors = FALSE)
map  <- read.csv("results/revision/sequencing/mapped_reads_per_sample.csv", stringsAsFactors = FALSE)
stopifnot(all(meta$sample_id %in% rc$Sample), nrow(meta) == 500)
rc <- rc[rc$Sample %in% meta$sample_id, ]
# Total_Bases / Total_Reads ~ 2 x read length => counts are read pairs
rc$bp_per_unit <- rc$Total_Bases / rc$Total_Reads
summ <- do.call(rbind, lapply(split(rc, rc$Step), function(d) data.frame(
  step = d$Step[1], n = nrow(d), mean_pairs = mean(d$Total_Reads), sd_pairs = sd(d$Total_Reads),
  min_pairs = min(d$Total_Reads), max_pairs = max(d$Total_Reads),
  mean_bp_per_pair = mean(d$bp_per_unit))))
qc <- rc[rc$Step == "QC", c("Sample", "Total_Reads")]
m  <- merge(map, qc, by.x = "sample_id", by.y = "Sample")
m$ratio_mapped_core_to_qc_pairs <- m$count / m$Total_Reads
summ_map <- data.frame(n = nrow(m), mean_ratio = mean(m$ratio_mapped_core_to_qc_pairs),
                       sd_ratio = sd(m$ratio_mapped_core_to_qc_pairs),
                       min_ratio = min(m$ratio_mapped_core_to_qc_pairs),
                       max_ratio = max(m$ratio_mapped_core_to_qc_pairs))
write.csv(summ, "results/revision/sequencing/atlas_read_depth_by_step.csv", row.names = FALSE)
write.csv(m, "results/revision/sequencing/mapped_vs_qc_per_sample.csv", row.names = FALSE)
write.csv(summ_map, "results/revision/sequencing/mapped_vs_qc_summary.csv", row.names = FALSE)
print(summ, digits = 4); print(summ_map, digits = 4)
