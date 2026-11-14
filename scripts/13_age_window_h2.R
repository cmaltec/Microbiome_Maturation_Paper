#!/usr/bin/env Rscript
#________________________________________
# - File Name : 13_age_window_h2.R
# - Purpose : Heritability of the per-animal descriptors in early (W3-W12) and late (W15-W24) windows (Fig 4A, Table S15)
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
# Sun Nov 01 15:30:00 2026 ------------------------------

setwd(here::here())

suppressPackageStartupMessages({ library(asreml); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 1999L
set.seed(20261101)

IN  <- "results/revision/inputs"
OUT <- "results/revision/window_h2"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

Ainv <- readRDS(file.path(IN, "Ainv_3gen.rds")); ped_ids <- attr(Ainv, "rowNames")
anim <- read.csv(file.path(IN, "animals.csv"), stringsAsFactors = FALSE)
lg   <- read.csv(gzfile("data/mag_abundance_long.csv.gz"), stringsAsFactors = FALSE)
curves <- read.csv("results/metrics/maturation_curves_per_sample.csv", stringsAsFactors = FALSE)
maz_s  <- read.csv("results/metrics/maz_per_sample.csv", stringsAsFactors = FALSE)
canon  <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv", stringsAsFactors = FALSE)

WEEK <- c(W3 = 3, W6 = 6, W9 = 9, W12 = 12, W15 = 15, W17 = 17, W21 = 21, W24 = 24)
mags <- sort(unique(lg$mag_id))
prof <- split(lg, paste(lg$animal_id, lg$timepoint))
get_profile <- function(a, tp) {
  x <- prof[[paste(a, tp)]]
  if (is.null(x)) return(NULL)
  setNames(x$count, x$mag_id)[mags]
}
clr <- function(x) { x[x == 0] <- 0.5; lx <- log(x); lx - mean(lx) }

w24 <- do.call(rbind, lapply(anim$animal_id, get_profile, tp = "W24"))
w24_c <- colMeans(w24); w24_p <- w24_c / sum(w24_c)

compute_metrics <- function(wmin, wmax, excl = numeric(0)) {
  out <- list()
  for (a in sort(anim$animal_id)) {
    tps <- names(WEEK)[WEEK >= wmin & WEEK <= wmax & !(WEEK %in% excl)]
    tps <- tps[vapply(tps, function(tp) !is.null(prof[[paste(a, tp)]]), logical(1))]
    n <- length(tps); if (n < 3) next
    wk <- WEEK[tps]
    mz <- maz_s[maz_s$animal_id == a & maz_s$week >= wmin & maz_s$week <= wmax & !(maz_s$week %in% excl), ]
    P  <- lapply(tps, function(tp) get_profile(a, tp))
    Pp <- lapply(P, function(p) p / sum(p))
    bc <- vapply(Pp, function(p) sum(abs(p - w24_p)) / 2, numeric(1))
    ai <- vapply(P, function(p) sqrt(sum((clr(p) - clr(w24_c))^2)), numeric(1))
    cb <- vapply(seq_len(n - 1), function(i) sum(abs(Pp[[i]] - Pp[[i + 1]])) / 2, numeric(1))
    vs <- if (length(cb) >= 3) coef(lm(cb ~ I((wk[-n] + wk[-1]) / 2)))[2] else NA
    pw <- c(); for (i in 1:(n - 1)) for (j in (i + 1):n) pw <- c(pw, sum(abs(Pp[[i]] - Pp[[j]])) / 2)
    sh <- curves$shannon[curves$animal_id == a][match(tps, curves$tp[curves$animal_id == a])]
    out[[a]] <- data.frame(
      animal_id = a, n_obs = n,
      maz_mean = if (nrow(mz) >= 2) mean(mz$maz) else NA,
      maz_slope = if (nrow(mz) >= 2) coef(lm(maz ~ week, mz))[2] else NA,
      bc_to_adult_mean = mean(bc), bc_to_adult_slope = coef(lm(bc ~ wk))[2],
      aitchison_to_adult_slope = coef(lm(ai ~ wk))[2],
      volatility_mean = mean(cb), volatility_slope = vs, intra_bc_var = var(pw),
      shannon_cv = sd(sh) / mean(sh), shannon_mean = mean(sh),
      shannon_slope = coef(lm(sh ~ wk))[2], row.names = NULL)
  }
  do.call(rbind, out)
}

traits <- c("maz_mean", "maz_slope", "bc_to_adult_mean", "bc_to_adult_slope",
            "aitchison_to_adult_slope", "volatility_mean", "volatility_slope",
            "intra_bc_var", "shannon_cv", "shannon_mean", "shannon_slope")

# Reimplementation check against the deposited full-window scalars
full <- compute_metrics(3, 24)
chk <- sapply(traits, function(t) {
  m <- merge(full[, c("animal_id", t)], canon[, c("animal_id", t)], by = "animal_id")
  c(n = nrow(m), max_abs_diff = max(abs(m[[2]] - m[[3]]), na.rm = TRUE),
    r = cor(m[[2]], m[[3]], use = "complete.obs"))
})
cat("Full-window reimplementation vs deposited scalars:\n"); print(t(chk), digits = 4)
write.csv(data.frame(trait = traits, t(chk)), file.path(OUT, "reimplementation_check.csv"), row.names = FALSE)

early <- compute_metrics(3, 12); late <- compute_metrics(15, 24)
write.csv(early, file.path(OUT, "scalars_W3_W12.csv"), row.names = FALSE)
write.csv(late,  file.path(OUT, "scalars_W15_W24.csv"), row.names = FALSE)
cat(sprintf("Animals: early %d, late %d\n", nrow(early), nrow(late)))

fit_vc <- function(df) {
  fit <- tryCatch(suppressWarnings(suppressMessages(asreml(
    fixed = y ~ 1, random = ~ litter_id + vm(animal_id, Ainv), residual = ~ idv(units),
    data = df, na.action = na.method(y = "include"), trace = FALSE, maxit = 200))),
    error = function(e) NULL)
  if (is.null(fit) || !fit$converge) return(c(h2 = NA, c2 = NA, se = NA))
  vc <- summary(fit)$varcomp; rn <- rownames(vc)
  ia <- grep("vm\\(animal_id", rn); ic <- which(rn == "litter_id"); ie <- grep("units!R", rn)
  s <- c(max(vc$component[ia], 0), max(vc$component[ic], 0), vc$component[ie]); tot <- sum(s)
  se <- tryCatch(vpredict(fit, as.formula(sprintf("h2 ~ V%d/(V%d+V%d+V%d)", ia, ic, ia, ie)))$SE,
                 error = function(e) NA)
  c(h2 = s[1] / tot, c2 = s[2] / tot, se = se)
}

run_window <- function(sc, label) {
  sc <- merge(sc, anim[, c("animal_id", "litter_id")], by = "animal_id")
  do.call(rbind, lapply(traits, function(tr) {
    df0 <- sc[!is.na(sc[[tr]]), c("animal_id", "litter_id", tr)]
    df0$y <- as.numeric(scale(df0[[tr]]))
    df0$litter_id <- factor(df0$litter_id)
    df <- df0; df$animal_id <- factor(df0$animal_id, levels = ped_ids)
    obs <- fit_vc(df)
    nul <- numeric(B)
    for (b in seq_len(B)) {
      dp <- df0; dp$animal_id <- factor(sample(df0$animal_id), levels = ped_ids)
      nul[b] <- fit_vc(dp)["h2"]
    }
    ok <- !is.na(nul)
    p <- (1 + sum(nul[ok] >= obs["h2"])) / (sum(ok) + 1)
    cat(sprintf("  %-8s %-26s n=%2d h2=%.3f (SE %.3f) p=%.4f null95=%.3f valid=%d\n",
                label, tr, nrow(df0), obs["h2"], obs["se"], p, quantile(nul[ok], 0.95), sum(ok)))
    data.frame(window = label, trait = tr, n = nrow(df0), h2 = obs["h2"], h2_se = obs["se"],
               c2 = obs["c2"], p_perm = p, null_q95 = quantile(nul[ok], 0.95),
               B_valid = sum(ok), row.names = NULL)
  }))
}

res <- rbind(run_window(early, "W3-W12"), run_window(late, "W15-W24"))
write.csv(res, file.path(OUT, "window_h2.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)

# Late window without the W21 visit
set.seed(20261101)
late_x <- compute_metrics(15, 24, excl = 21)
write.csv(late_x, file.path(OUT, "scalars_W15_W24_noW21.csv"), row.names = FALSE)
cat(sprintf("Animals: late without W21 %d\n", nrow(late_x)))
res2 <- run_window(late_x, "W15-W24 excl W21")
write.csv(res2, file.path(OUT, "window_h2_late_noW21.csv"), row.names = FALSE)
print(res2, digits = 3, row.names = FALSE)

# ---- Table S15 -----------------------------------------------------------------
w  <- read.csv("results/revision/window_h2/window_h2.csv", stringsAsFactors = FALSE)
w2 <- read.csv("results/revision/window_h2/window_h2_late_noW21.csv", stringsAsFactors = FALSE)
a  <- rbind(w, w2)
lab <- c(maz_mean = "MAZ mean", maz_slope = "MAZ slope", bc_to_adult_mean = "BC mean",
         bc_to_adult_slope = "BC slope", aitchison_to_adult_slope = "Aitchison slope",
         volatility_mean = "Volatility", volatility_slope = "Volatility slope",
         intra_bc_var = "Intra-BC variance", shannon_cv = "Shannon CV",
         shannon_mean = "Shannon mean", shannon_slope = "Shannon slope")
fmt <- function(r) {
  if (is.na(r$h2)) return(c("n.e.", "—", as.character(r$n)))
  se <- if (is.na(r$h2_se) || r$h2 < 1e-4) "—" else sprintf("%.3f", r$h2_se)
  p  <- sprintf("%.4f", r$p_perm)
  c(sprintf("%.3f (%s)", r$h2, se), p, as.character(r$n))
}
wins <- c("W3-W12", "W15-W24", "W15-W24 excl W21")
tab <- do.call(rbind, lapply(names(lab), function(t) {
  v <- unlist(lapply(wins, function(x) fmt(a[a$trait == t & a$window == x, ])))
  data.frame(Descriptor = lab[[t]], t(v), row.names = NULL)
}))
names(tab) <- c("Descriptor", as.vector(outer(c("h2 (SE)", "p", "n"), wins,
                                              function(a, b) paste(b, a))))
write.csv(tab, "results/revision/window_h2/table_S15.csv", row.names = FALSE)
print(tab, row.names = FALSE)
