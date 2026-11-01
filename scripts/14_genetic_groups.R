#!/usr/bin/env Rscript
#________________________________________
# - File Name : 14_genetic_groups.R
# - Purpose : Animal model with genetic groups for unknown parents vs selection line as fixed (Table S2b footnote)
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

IN  <- "results/revision/inputs"
OUT <- "results/revision/genetic_groups"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

Ainv3 <- readRDS(file.path(IN, "Ainv_3gen.rds"))
Ainv2 <- readRDS(file.path(IN, "Ainv_2gen.rds"))
ped3  <- read.csv(file.path(IN, "pedigree_3gen_3col.csv"), stringsAsFactors = FALSE)
anim  <- read.csv(file.path(IN, "animals.csv"), stringsAsFactors = FALSE)
dat   <- read.csv(file.path(IN, "per_animal_traits.csv"), stringsAsFactors = FALSE)
p3raw <- read.csv("data/pedigree_3generation.csv", stringsAsFactors = FALSE)

sire_line <- distinct(anim[, c("sire_id", "line")])
dam_line  <- distinct(data.frame(dam = anim$litter_id, line = anim$line))
stopifnot(!anyDuplicated(sire_line$sire_id), !anyDuplicated(dam_line$dam))

# Q matrix by recursion: founders take the mean of their two phantom-parent
# groups; others take the mean of their parents' Q rows.
q_matrix <- function(ped, founder_group) {
  groups <- sort(unique(unlist(founder_group)))
  Q <- matrix(0, nrow(ped), length(groups), dimnames = list(ped$id, groups))
  for (i in seq_len(nrow(ped))) {
    s <- ped$sire[i]; d <- ped$dam[i]
    qs <- if (is.na(s)) setNames(as.numeric(groups == founder_group[[ped$id[i]]][1]), groups) else Q[s, ]
    qd <- if (is.na(d)) setNames(as.numeric(groups == founder_group[[ped$id[i]]][2]), groups) else Q[d, ]
    Q[i, ] <- 0.5 * qs + 0.5 * qd
  }
  Q
}

# GG3
fg3 <- list()
for (id in ped3$id) {
  if (id %in% sire_line$sire_id) {
    g <- sire_line$line[sire_line$sire_id == id]; fg3[[id]] <- c(g, g)
  } else fg3[[id]] <- c("BASE", "BASE")
}
Q3 <- q_matrix(ped3, fg3)

# GG2
ped2 <- rbind(data.frame(id = c(sire_line$sire_id, dam_line$dam), sire = NA, dam = NA),
              data.frame(id = anim$animal_id, sire = anim$sire_id, dam = anim$litter_id))
fg2 <- list()
for (id in ped2$id) {
  g <- if (id %in% sire_line$sire_id) sire_line$line[sire_line$sire_id == id] else
    if (id %in% dam_line$dam) dam_line$line[dam_line$dam == id] else NA
  fg2[[id]] <- c(g, g)
}
Q2 <- q_matrix(ped2, fg2)

qa3 <- Q3[anim$animal_id, , drop = FALSE]; qa2 <- Q2[anim$animal_id, , drop = FALSE]
write.csv(data.frame(animal_id = anim$animal_id, line = anim$line, qa3), file.path(OUT, "Q_GG3.csv"), row.names = FALSE)
write.csv(data.frame(animal_id = anim$animal_id, line = anim$line, qa2), file.path(OUT, "Q_GG2.csv"), row.names = FALSE)

# Column-space check: is [1, Q] the same space as [1, line]?
X_line <- cbind(1, as.numeric(anim$line == "TOL"))
span_eq <- function(Q) {
  X <- cbind(1, Q)
  c(rank_Q = qr(X)$rank, rank_joint = qr(cbind(X, X_line))$rank, rank_line = qr(X_line)$rank)
}
sp <- rbind(GG3 = span_eq(qa3), GG2 = span_eq(qa2))
cat("Rank checks (equal ranks => group fixed effects span the line contrast):\n"); print(sp)
cat("\nGG3 Q rows by line (unique):\n"); print(unique(data.frame(line = anim$line, round(qa3, 4))))
cat("\nGG2 Q rows by line (unique):\n"); print(unique(data.frame(line = anim$line, round(qa2, 4))))
write.csv(data.frame(model = rownames(sp), sp), file.path(OUT, "group_rank_check.csv"), row.names = FALSE)

traits <- c("maz_mean", "maz_slope", "bc_to_adult_mean", "bc_to_adult_slope",
            "aitchison_to_adult_slope", "rel_maturity_w24", "volatility_mean",
            "volatility_slope", "intra_bc_var", "shannon_cv", "shannon_mean",
            "shannon_w24", "shannon_slope", "BC_FPC1", "BC_FPC2", "BC_FPC3",
            "shan_FPC1", "shan_FPC2", "shan_FPC3")

dat <- inner_join(dat, data.frame(animal_id = anim$animal_id,
                                  GG3_SEN = qa3[, "SEN"], GG3_TOL = qa3[, "TOL"],
                                  GG2_TOL = qa2[, "TOL"]), by = "animal_id")
dat$line <- factor(dat$line); dat$litter_id <- factor(dat$litter_id)

fit_h2 <- function(tr, fixed_rhs, Ainv) {
  df <- dat
  df$y <- df[[tr]]
  df$animal_id <- factor(df$animal_id, levels = attr(Ainv, "rowNames"))
  assign("Ainv_cur", Ainv, envir = .GlobalEnv)
  fit <- tryCatch(suppressWarnings(suppressMessages(asreml(
    fixed = as.formula(paste("y ~", fixed_rhs)),
    random = ~ litter_id + vm(animal_id, Ainv_cur), residual = ~ idv(units),
    data = df, na.action = na.method(y = "include"), trace = FALSE, maxit = 200))),
    error = function(e) NULL)
  if (is.null(fit) || !fit$converge) return(c(h2 = NA, se = NA, logLik = NA))
  vc <- summary(fit)$varcomp; rn <- rownames(vc)
  ia <- grep("vm\\(animal_id", rn); ic <- which(rn == "litter_id"); ie <- grep("units!R|units!units", rn)[1]
  sA <- max(vc$component[ia], 0); tot <- sA + max(vc$component[ic], 0) + vc$component[ie]
  vp <- tryCatch(vpredict(fit, as.formula(sprintf("h2 ~ V%d/(V%d+V%d+V%d)", ia, ic, ia, ie))),
                 error = function(e) data.frame(SE = NA))
  c(h2 = sA / tot, se = vp$SE[1], logLik = fit$loglik)
}

res <- do.call(rbind, lapply(traits, function(tr) {
  b  <- fit_h2(tr, "1", Ainv3)
  l  <- fit_h2(tr, "1 + line", Ainv3)
  g3 <- fit_h2(tr, "1 + GG3_SEN + GG3_TOL", Ainv3)
  g2 <- fit_h2(tr, "1 + GG2_TOL", Ainv2)
  l2 <- fit_h2(tr, "1 + line", Ainv2)
  data.frame(trait = tr,
             h2_primary = b["h2"], se_primary = b["se"],
             h2_line_fixed = l["h2"], se_line_fixed = l["se"],
             h2_GG3 = g3["h2"], se_GG3 = g3["se"],
             logLik_line = l["logLik"], logLik_GG3 = g3["logLik"],
             h2_GG2_A2 = g2["h2"], h2_line_fixed_A2 = l2["h2"])
}))
rownames(res) <- NULL
write.csv(res, file.path(OUT, "h2_genetic_groups.csv"), row.names = FALSE)
cat("\n"); print(res, digits = 3)
