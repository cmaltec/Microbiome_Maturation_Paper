#!/usr/bin/env Rscript
#________________________________________
# - File Name : 12_figures.R
# - Purpose : All manuscript and supplementary figures
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
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(png); library(grid); library(Matrix); library(reshape2); library(ggrepel)
  library(ggh4x); library(ggnewscale)
})
source("scripts/palette.R")

week_map <- c(W3=3, W6=6, W9=9, W12=12, W15=15, W17=17, W21=21, W24=24)

# ---- Core data ----
curves   <- read.csv("results/metrics/maturation_curves_per_sample.csv", stringsAsFactors=FALSE)
metrics  <- read.csv("results/metrics/maturation_all_metrics_per_animal.csv", stringsAsFactors=FALSE)
h2       <- read.csv("results/metrics/maturation_h2_animal_model.csv", stringsAsFactors=FALSE)
thr      <- read.csv("results/metrics/maturation_threshold_metrics_imputed.csv", stringsAsFactors=FALSE)
h2_thr   <- read.csv("results/metrics/maturation_h2_threshold.csv", stringsAsFactors=FALSE)
ped      <- read.csv("data/metadata/pedigree.csv", stringsAsFactors=FALSE)
samp_meta <- read.csv("data/sample_metadata.csv", stringsAsFactors=FALSE)

# ---- h² and h²(t) data for Figs 3 and 4 ----
win_h2   <- read.csv("results/revision/window_h2/window_h2.csv", stringsAsFactors=FALSE)
win_n    <- c(early = nrow(read.csv("results/revision/window_h2/scalars_W3_W12.csv")),
              late  = nrow(read.csv("results/revision/window_h2/scalars_W15_W24.csv")))
sire_blu <- read.csv("results/metrics/maturation_sire_blups.csv", stringsAsFactors=FALSE)
w6_h2    <- tryCatch(read.csv("results/metrics/maturation_h2_w6_scalars.csv", stringsAsFactors=FALSE),
                     error=function(e) NULL)
perm     <- read.csv("results/metrics/maturation_h2_permutation_test.csv", stringsAsFactors=FALSE)
perm_w6  <- read.csv("results/metrics/maturation_h2_permutation_test_w6.csv", stringsAsFactors=FALSE)

# ---- Fig 6 inputs: per-feature h² across aggregation levels ----
mag_h2   <- read.csv("results/metrics/mag_per_feature_h2.csv", stringsAsFactors=FALSE)
ko_h2    <- read.csv("results/metrics/ko_per_feature_h2.csv", stringsAsFactors=FALSE)
mod_h2   <- read.csv("results/metrics/ko_per_module_h2.csv", stringsAsFactors=FALSE)

curves$week <- week_map[curves$tp]
cat("Data loaded.\n")

# FIG 1 — Design and data scope
#   Panel A: external NC-Design-I mating schematic (fig1_mating_design.png)
#   Panel B: pedigree A-matrix heatmap (78×78, animals ordered line→sire→litter)
#   Panel C: per-animal distributions of the six class-representative metrics, by line
cat("Fig 1 (A bundled + B pedigree + C six-class violins)...\n")

# Panel A — programmatic NC-Design-I schematic built from ped + sample metadata
# Three tiers: 7 sires → 13 litters → 78 progeny, with actual progeny counts.
ped_full <- ped %>%
  arrange(line, sire_id, litter_id, animal_id) %>%
  mutate(x_progeny = dplyr::row_number())
litter_pos <- ped_full %>%
  group_by(line, sire_id, litter_id) %>%
  summarise(x_litter = mean(x_progeny),
            n_litter = dplyr::n(), .groups = "drop") %>%
  arrange(line, sire_id, litter_id)
sire_pos <- ped_full %>%
  group_by(line, sire_id) %>%
  summarise(x_sire = mean(x_progeny),
            n_sire = dplyr::n(),
            n_litters = dplyr::n_distinct(litter_id), .groups = "drop") %>%
  arrange(line, sire_id) %>%
  mutate(sire_lbl = sprintf("S%d", dplyr::row_number()),
         sire_lbl = paste0(sire_lbl, "\n(n=", n_sire, ")"))
# Build short sire labels per sire_id for the segment data

# Edge data for sire→litter and litter→progeny connectors
edges_sire_litter <- litter_pos %>%
  left_join(sire_pos %>% select(sire_id, x_sire), by = "sire_id") %>%
  mutate(line_grp = line)
edges_litter_progeny <- ped_full %>%
  left_join(litter_pos %>% select(litter_id, x_litter), by = "litter_id")

y_sire <- 3; y_litter <- 2; y_progeny <- 1
sire_pos$y    <- y_sire
litter_pos$y  <- y_litter
ped_full$y    <- y_progeny

f1a <- ggplot() +
  # Sire → litter connectors
  geom_segment(data = edges_sire_litter,
               aes(x = x_sire, xend = x_litter,
                   y = y_sire - 0.30, yend = y_litter + 0.18,
                   color = line),
               linewidth = 0.4, alpha = 0.55) +
  # Litter → progeny connectors
  geom_segment(data = edges_litter_progeny,
               aes(x = x_litter, xend = x_progeny,
                   y = y_litter - 0.18, yend = y_progeny + 0.10,
                   color = line),
               linewidth = 0.25, alpha = 0.30) +
  # Sire boxes
  geom_tile(data = sire_pos,
            aes(x = x_sire, y = y_sire, fill = line),
            width = 3.6, height = 0.55, color = "grey20", linewidth = 0.45) +
  geom_text(data = sire_pos,
            aes(x = x_sire, y = y_sire, label = sire_lbl),
            size = 3.0, color = "white", fontface = "bold", lineheight = 0.85) +
  # Litter circles (size = n_progeny in litter)
  geom_point(data = litter_pos,
             aes(x = x_litter, y = y_litter, fill = line, size = n_litter),
             shape = 21, color = "grey20", stroke = 0.4) +
  geom_text(data = litter_pos,
            aes(x = x_litter, y = y_litter, label = n_litter),
            size = 2.6, color = "white", fontface = "bold") +
  # Progeny dots
  geom_point(data = ped_full,
             aes(x = x_progeny, y = y_progeny, color = line),
             size = 1.3, alpha = 0.85) +
  # Row labels
  annotate("text", x = -1.6, y = y_sire,    label = "Sires (7)",
           size = 3.4, hjust = 1, fontface = "bold", color = "grey15") +
  annotate("text", x = -1.6, y = y_litter,  label = "Litters (13)",
           size = 3.4, hjust = 1, fontface = "bold", color = "grey15") +
  annotate("text", x = -1.6, y = y_progeny, label = "Progeny (78)",
           size = 3.4, hjust = 1, fontface = "bold", color = "grey15") +
  # Line block headers
  annotate("text",
           x = mean(sire_pos$x_sire[sire_pos$line == "TOL"]),
           y = y_sire + 0.55, label = "TOL (3 sires · 7 litters · 39 progeny)",
           size = 3.2, fontface = "bold", color = line_fills["TOL"]) +
  annotate("text",
           x = mean(sire_pos$x_sire[sire_pos$line == "SEN"]),
           y = y_sire + 0.55, label = "SEN (4 sires · 6 litters · 39 progeny)",
           size = 3.2, fontface = "bold", color = line_fills["SEN"]) +
  scale_fill_manual (values = line_fills, guide = "none") +
  scale_color_manual(values = line_fills, guide = "none") +
  scale_size_continuous(range = c(4.5, 9.5), guide = "none") +
  scale_x_continuous(limits = c(-10, 80), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 3.9)) +
  labs(title = "NC Design I — sires nested above litters nested above progeny") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0,
                                  margin = margin(b = 4)))

# Panel B — sampling completeness by line × week (the SEN-W3 gap that motivates
# PACE imputation and the W6–W24 balanced refit throughout the paper).
line_n <- ped %>% group_by(line) %>% summarise(n_line = n(), .groups="drop")
samp_meta$week <- week_map[samp_meta$timepoint]
obs_grid <- samp_meta %>%
  group_by(line, week) %>%
  summarise(n_obs = dplyr::n_distinct(animal_id), .groups="drop") %>%
  left_join(line_n, by="line") %>%
  mutate(pct_obs = 100 * n_obs / n_line,
         label   = sprintf("%d / %d", n_obs, n_line))

# Highlight the W6–W24 analytical window used for the per-trait line tests.
# Use position rank (1..8) on x so the rect annotation aligns with the dodged bars.
obs_grid$xpos <- as.integer(factor(obs_grid$week,
                                   levels=c(3,6,9,12,15,17,21,24)))
# Ghost layer: PACE-imputed coverage reaches 100% at every (line × week) cell.
# This shows what is filled in by imputation, on top of which the observed
# bars sit. Per-trait line tests use observed-only W6–W24 (cream window);
# trajectory + threshold analyses use the imputed full W3–W24 window (lavender).
ghost <- expand.grid(line = c("SEN", "TOL"),
                     xpos = 1:8, stringsAsFactors = FALSE)
ghost$pct_obs <- 100

f1b <- ggplot(obs_grid, aes(x=xpos, y=pct_obs, fill=line)) +
  # Outer band: PACE-imputed coverage spans the full W3–W24 range
  annotate("rect", xmin=0.5, xmax=8.5, ymin=0, ymax=100,
           fill="#ece4f5", alpha=0.55, color=NA) +
  # Inner band: observed-only line-test window
  annotate("rect", xmin=1.5, xmax=8.5, ymin=0, ymax=100,
           fill="#fff3b0", alpha=0.85,
           color="grey25", linewidth=0.45) +
  # Imputation bracket on top (full W3–W24)
  annotate("segment", x=0.5, xend=8.5, y=112, yend=112,
           color="#6a4a9c", linewidth=0.8) +
  annotate("segment", x=0.5, xend=0.5, y=109, yend=115,
           color="#6a4a9c", linewidth=0.8) +
  annotate("segment", x=8.5, xend=8.5, y=109, yend=115,
           color="#6a4a9c", linewidth=0.8) +
  annotate("text", x=4.5, y=119,
           label="PACE-imputed → 78/78 at every week (trajectory + threshold + RRM analyses)",
           size=3.0, color="#3a2566", fontface="bold") +
  # Line-test window bracket below the imputation bracket
  annotate("segment", x=1.5, xend=8.5, y=104, yend=104,
           color="grey20", linewidth=0.7) +
  annotate("segment", x=1.5, xend=1.5, y=102, yend=106,
           color="grey20", linewidth=0.7) +
  annotate("segment", x=8.5, xend=8.5, y=102, yend=106,
           color="grey20", linewidth=0.7) +
  annotate("text", x=5.0, y=108,
           label="Observed-only window — per-trait line tests (W6–W24)",
           size=2.9, color="grey15", fontface="bold") +
  # Ghost (imputed-augmentation) layer
  geom_col(data=ghost, aes(x=xpos, y=pct_obs, fill=line),
           position=position_dodge(0.78), width=0.72,
           alpha=0.25, color=NA, show.legend=FALSE) +
  # Observed coverage on top
  geom_col(position=position_dodge(0.78), width=0.72,
           color="white", linewidth=0.3) +
  geom_text(aes(label=label), position=position_dodge(0.78),
            vjust=-0.4, size=2.8, color="grey25") +
  scale_fill_manual(values=line_fills, name="Line") +
  scale_x_continuous(breaks=1:8,
                     labels=as.character(c(3,6,9,12,15,17,21,24)),
                     expand=expansion(add=0.5)) +
  scale_y_continuous(limits=c(0, 125), breaks=seq(0, 100, 25)) +
  labs(x="Week", y="% of line observed",
       title="Sampling completeness — observed bars; PACE fills the gaps",
       subtitle="500 of 624 cells observed (80.1%); faint extension = PACE imputation augmenting each cell to 78/78") +
  theme_paper2(base_size=10) +
  theme(plot.title=element_text(face="bold", size=11,
                                margin=margin(b=8)),
        plot.subtitle=element_text(size=8, margin=margin(t=2, b=6)),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.box.margin=margin(t=4))

# Panel C — per-animal distributions of the six class-representative metrics, by line
# Pull t70 from PACE-imputed threshold table; rest from observed per-animal metrics.
focal_cols <- list(
  "MAZ mean"       = list(df = metrics, col = "maz_mean"),
  "BC slope"       = list(df = metrics, col = "bc_to_adult_slope"),
  "t70 (wk)"       = list(df = thr,     col = "t70"),
  "Intra-BC var"   = list(df = metrics, col = "intra_bc_var"),
  "Shannon slope"  = list(df = metrics, col = "shannon_slope"),
  "BC FPC2"        = list(df = metrics, col = "BC_FPC2"))

panel_c_rows <- do.call(rbind, lapply(names(focal_cols), function(lab) {
  rec <- focal_cols[[lab]]
  if (!rec$col %in% names(rec$df)) return(NULL)
  d <- rec$df[, c("animal_id", rec$col)]
  names(d) <- c("animal_id", "value")
  d$line <- ped$line[match(d$animal_id, ped$animal_id)]
  d$metric <- factor(lab, levels=names(focal_cols))
  d
}))

f1c <- ggplot(panel_c_rows, aes(x=line, y=value, fill=line)) +
  geom_violin(scale="width", width=0.85, alpha=0.55, color=NA) +
  geom_point(position=position_jitter(width=0.12, seed=20260518),
             size=0.9, alpha=0.6, color="grey20") +
  facet_wrap(~ metric, scales="free_y", ncol=3) +
  scale_fill_manual(values=line_fills, guide="none") +
  labs(x="", y="Per-animal value",
       title="Distribution of the six class-representative metrics, by line") +
  theme_paper2(base_size=10) +
  theme(strip.text=element_text(face="bold", size=9),
        plot.title=element_text(face="bold", size=11),
        axis.text.x=element_text(size=9))

fig1 <- f1a / (f1b | f1c) + plot_layout(heights=c(1, 1.3)) +
  plot_annotation(tag_levels = "A")
ggsave("paper/fig1_panels_BC.png", fig1, width=14, height=12, dpi=300)

# FIG 2 — Maturation as a quantitative trait (Obj 1)
#   A: Population BC-to-adult arc with W3–W17 Gompertz + W3–W24 modified
#      Gompertz-with-dip (R²~0.96); per-animal spaghetti underneath.
#   B: Per-animal maturation rate (bc_to_adult slope) — the heritable scalar
#      decomposed in Figs 3–4.
#   C: t70 / t75 threshold-timing distributions by line — the SEN-vs-TOL
#      early-stage gap and the input to the §7 selection-response calc.
cat("Fig 2...\n")

pop_bc <- curves %>% group_by(tp, week) %>%
  summarise(mean=mean(bc_to_adult), se=sd(bc_to_adult)/sqrt(dplyr::n()), .groups="drop")
curves$maturity <- 1 - curves$bc_to_adult
pop_mat <- curves %>% group_by(week) %>%
  summarise(mean=mean(maturity), .groups="drop")

# Vanilla Gompertz on W3–W17 maturity
gomp_pop <- tryCatch(
  nls(mean ~ A*exp(-exp(-k*(week-ti))),
      data=pop_mat[pop_mat$week<=17, ],
      start=list(A=0.7, k=0.15, ti=5)),
  error=function(e) NULL)

# Modified Gompertz with explicit Gaussian dip at W21 — see §1, §2 of narrative.
# Anchor mu=21, fit A, k, ti, C, sigma to the 8-week mean trajectory.
mod_gomp <- tryCatch(
  nls(mean ~ A*exp(-exp(-k*(week-ti))) -
              C*exp(-(week-21)^2 / (2*sigma^2)),
      data=pop_mat,
      start=list(A=0.7, k=0.15, ti=-1, C=0.14, sigma=1.5)),
  error=function(e) NULL)

# Pre-compute the fit curves and R² so we can place them in the subtitle.
pred_mod <- NULL; r2_mod <- NA_real_
if (!is.null(mod_gomp)) {
  pred_mod <- data.frame(week=seq(3, 24, 0.25))
  pred_mod$y <- 1 - predict(mod_gomp, newdata=pred_mod)
  r2_mod <- 1 - sum(residuals(mod_gomp)^2) /
                sum((pop_mat$mean - mean(pop_mat$mean))^2)
}
pred_van_fit <- NULL; pred_van_ext <- NULL
if (!is.null(gomp_pop)) {
  pred_van_fit <- data.frame(week=seq(3, 17, 0.25))
  pred_van_fit$y <- 1 - predict(gomp_pop, newdata=pred_van_fit)
  pred_van_ext <- data.frame(week=seq(17, 24, 0.25))
  pred_van_ext$y <- 1 - predict(gomp_pop, newdata=pred_van_ext)
}

fit_color <- "#d97826"  # warm orange, distinct from SEN/TOL palette
van_color <- "#5b6a8a"  # muted slate-blue for the simpler model

f2a <- ggplot() +
  # 78 per-animal trajectories, very light grey
  geom_line(data=curves, aes(x=week, y=bc_to_adult, group=animal_id),
            color="grey75", alpha=0.30, linewidth=0.30) +
  # Vertical marker at W21 to anchor the dip
  geom_vline(xintercept=21, lty=3, color="grey55", linewidth=0.45) +
  # Population mean + SE
  geom_ribbon(data=pop_bc, aes(x=week, ymin=mean-se, ymax=mean+se),
              alpha=0.22, fill="grey45") +
  geom_line (data=pop_bc, aes(x=week, y=mean), linewidth=1.6, color="grey15") +
  geom_point(data=pop_bc, aes(x=week, y=mean), size=2.6, color="grey15") +
  scale_x_continuous(breaks=c(3,6,9,12,15,17,21,24)) +
  labs(x="Week", y="Bray-Curtis to adult centroid",
       title="Population maturation arc",
       subtitle=if (!is.na(r2_mod))
         sprintf("Gompertz W3–W17 (blue, dashed → projection); Gompertz-with-dip (orange, R² = %.3f)",
                 r2_mod) else
         "Per-animal trajectories (grey); mean ± SE (dark)") +
  theme_paper2()

# Overlay vanilla Gompertz (W3–W17 fit + W17–W24 projection)
if (!is.null(pred_van_fit)) {
  f2a <- f2a +
    geom_line(data=pred_van_fit, aes(x=week, y=y),
              linewidth=0.85, color=van_color, lty=2) +
    geom_line(data=pred_van_ext, aes(x=week, y=y),
              linewidth=0.85, color=van_color, lty=3)
}
# Overlay modified Gompertz-with-dip across the full W3–W24 range
if (!is.null(pred_mod)) {
  f2a <- f2a +
    geom_line(data=pred_mod, aes(x=week, y=y),
              linewidth=1.3, color=fit_color, alpha=0.95)
}

# Panel B — per-animal maturation rate distribution (heritable scalar)
rate_df <- metrics %>%
  transmute(animal_id, bc_slope = bc_to_adult_slope,
            line = ped$line[match(animal_id, ped$animal_id)]) %>%
  filter(!is.na(bc_slope))
rate_df <- rate_df[order(rate_df$bc_slope), ]
rate_df$rank <- seq_len(nrow(rate_df))

f2b <- ggplot(rate_df, aes(x=rank, y=bc_slope, color=line, fill=line)) +
  geom_hline(yintercept=0, lty=3, color=pal$ref_line) +
  geom_hline(yintercept=median(rate_df$bc_slope, na.rm=TRUE),
             lty=2, color="grey40", linewidth=0.4) +
  geom_segment(aes(xend=rank, yend=median(rate_df$bc_slope, na.rm=TRUE)),
               linewidth=0.35, alpha=0.55) +
  geom_point(size=1.7, alpha=0.9, shape=21, stroke=0.3, color="grey15") +
  scale_color_manual(values=line_fills, guide="none") +
  scale_fill_manual (values=line_fills, name="Line") +
  labs(x="Animal rank (fast → slow convergence)",
       y="BC-to-adult slope (ΔBC / week)",
       title="Per-animal maturation rate — heritable scalar",
       subtitle="h² = 0.28 ± 0.21 (W6+ refit 0.35); more negative slope = faster convergence") +
  theme_paper2() +
  theme(legend.position=c(0.85, 0.18))

# Panel C — t70 / t75 threshold-timing distributions by line
thr_lab <- c(t70="t70 — 70% maturity (wk)", t75="t75 — 75% maturity (wk)")
thr_long <- thr %>%
  select(animal_id, line, t70, t75) %>%
  pivot_longer(cols=c(t70, t75), names_to="threshold", values_to="week_at_thr") %>%
  filter(!is.na(week_at_thr)) %>%
  mutate(thr_label = factor(thr_lab[threshold], levels=thr_lab))

f2c <- ggplot(thr_long, aes(x=line, y=week_at_thr, fill=line)) +
  geom_violin(scale="width", width=0.85, alpha=0.55, color=NA) +
  geom_point(position=position_jitter(width=0.12, seed=20260518),
             size=1.2, alpha=0.7, color="grey20") +
  stat_summary(fun=median, geom="crossbar", width=0.35,
               color="grey15", linewidth=0.4) +
  facet_wrap(~ thr_label, ncol=2) +
  scale_fill_manual(values=line_fills, guide="none") +
  labs(x="", y="Week at threshold (PACE-imputed)",
       title="Threshold timing by line",
       subtitle="t70 / t75 — heritable rate descriptors (t70 h² = 0.26 ± 0.21, t75 h² = 0.27 ± 0.23)") +
  theme_paper2() +
  theme(strip.text=element_text(face="bold", size=9))

# Apply a shared title/subtitle theme to all three panels so the headers line
# up vertically and don't squeeze against the plot bodies. Flush left.
fig2_header_theme <- theme(
  plot.title    = element_text(face="bold", size=11, hjust=0,
                               margin=margin(t=0, b=6)),
  plot.subtitle = element_text(size=8.5, hjust=0,
                               margin=margin(t=0, b=8)),
  plot.title.position    = "plot",
  plot.caption.position  = "plot",
  plot.margin            = margin(t=6, r=10, b=6, l=10))
fig2 <- (f2a + f2b + f2c + plot_layout(ncol=3) +
         plot_annotation(tag_levels = "A")) & fig2_header_theme
ggsave("paper/fig2_maturation.png", fig2, width=17, height=5.8, dpi=300)

# FIG 3 — Per-trait heritability (Obj 2A LEAD)
cat("Fig 3...\n")

drop_metrics <- c("gompertz_k","gompertz_ti","gompertz_A","logistic_k","logistic_t0")
h2 <- h2[!h2$metric %in% drop_metrics, ]

# Append Alt-C threshold heritabilities — include early thresholds too so the
# "early-stage thresholds not heritable, late ones are" pattern is visible.
h2_thr_keep <- h2_thr[h2_thr$metric %in% c("t55","t60","t65","t70","t75"), ]
h2_thr_rows <- data.frame(
  metric=h2_thr_keep$metric, class="3-Time-to-threshold", n=h2_thr_keep$n,
  s2_A=h2_thr_keep$h2 / pmax(1 - h2_thr_keep$c2, 0.001),
  s2_c=h2_thr_keep$c2, s2_e=1 - h2_thr_keep$h2 - h2_thr_keep$c2, s2_total=1,
  h2=h2_thr_keep$h2, h2_se=h2_thr_keep$h2_se,
  c2=h2_thr_keep$c2, c2_se=NA, e2=1 - h2_thr_keep$h2 - h2_thr_keep$c2,
  lrt_animal_chisq=NA, lrt_animal_p=h2_thr_keep$lrt_p,
  lrt_litter_chisq=NA, lrt_litter_p=NA, converged=TRUE,
  lrt_animal_p_meff=NA, sig_animal=FALSE, stringsAsFactors=FALSE)
h2_thr_rows <- h2_thr_rows[, intersect(names(h2_thr_rows), names(h2))]
for (cc in setdiff(names(h2), names(h2_thr_rows))) h2_thr_rows[[cc]] <- NA
h2_thr_rows <- h2_thr_rows[, names(h2)]
h2 <- rbind(h2, h2_thr_rows)

h2_labels <- c(
  maz_mean="MAZ mean", maz_slope="MAZ slope", bc_to_adult_mean="BC mean",
  bc_to_adult_slope="BC slope", aitchison_slope="Aitchison slope",
  rel_maturity_w24="Rel. maturity W24",
  t55="t55", t60="t60", t65="t65", t70="t70", t75="t75",
  volatility_mean="Volatility", volatility_slope="Vol. slope",
  intra_bc_var="Intra-BC var", shannon_cv="Shannon CV",
  shannon_mean="Shannon mean", shannon_w24="Shannon W24",
  shannon_slope="Shannon slope",
  BC_FPC1="BC FPC1", BC_FPC2="BC FPC2", BC_FPC3="BC FPC3",
  shan_FPC1="Shannon FPC1", shan_FPC2="Shannon FPC2", shan_FPC3="Shannon FPC3")
h2_class <- c(
  maz_mean="1-Age", maz_slope="1-Age",
  bc_to_adult_mean="2-Distance", bc_to_adult_slope="2-Distance",
  aitchison_slope="2-Distance", rel_maturity_w24="2-Distance",
  t55="3-Time-to-threshold", t60="3-Time-to-threshold", t65="3-Time-to-threshold",
  t70="3-Time-to-threshold", t75="3-Time-to-threshold",
  volatility_mean="4-Stability", volatility_slope="4-Stability",
  intra_bc_var="4-Stability", shannon_cv="4-Stability",
  shannon_mean="5-Diversity", shannon_w24="5-Diversity", shannon_slope="5-Diversity",
  BC_FPC1="6-FPCA", BC_FPC2="6-FPCA", BC_FPC3="6-FPCA",
  shan_FPC1="6-FPCA", shan_FPC2="6-FPCA", shan_FPC3="6-FPCA")
h2$label <- h2_labels[h2$metric]
h2$mclass <- h2_class[h2$metric]
h2 <- h2[!is.na(h2$label), ]

# Restrict the Fig 3 forest plot to informative h² estimates: drop boundary
# REML h²≈0 estimates AND uninformative metrics where h²/SE < 0.5. The full
# panel (including all zero and uninformative metrics) is in Supplementary
# Table S2. This keeps the lead figure focused on metrics with both
# detectable additive genetic variance and informative precision.
h2_drop_mask <- h2$h2 < 0.01 | is.na(h2$h2) | (h2$h2 / pmax(h2$h2_se, 1e-6)) < 0.5
h2_dropped_metrics <- h2$metric[h2_drop_mask]
cat(sprintf("  Fig 3 forest plot: dropping %d uninformative metrics (h²≈0 or h²/SE < 0.5; see Table S2): %s\n",
            length(h2_dropped_metrics), paste(h2_dropped_metrics, collapse=", ")))
h2 <- h2[!h2_drop_mask, ]

class_short <- c("1-Age"="Age prediction", "2-Distance"="Distance to adult",
                 "3-Time-to-threshold"="Time-to-threshold", "4-Stability"="Stability",
                 "5-Diversity"="Diversity", "6-FPCA"="FPCA modes")
h2$class_label <- class_short[h2$mclass]
# Class facet ordering: FPCA at top, Age at bottom.
h2$class_label <- factor(h2$class_label, levels=rev(unique(class_short[sort(unique(h2$mclass))])))

# Sort rows so that within each class, h² descends (highest first).
h2 <- h2 %>% arrange(class_label, desc(h2))

# Set label factor levels EXPLICITLY in bottom-to-top order so that ggplot's
# default discrete-axis behavior (factor level 1 at bottom) places the
# highest-h² metric of each class at the top of its facet panel.
h2$label <- factor(h2$label, levels = rev(as.character(h2$label)))

# W6+ refit overlay table (paired with full-window h² where available)
w6_map <- c(maz_mean="maz_mean", maz_slope="maz_slope",
            intra_bc_var="intra_bc_var", volatility_mean="volatility_mean",
            shannon_slope="shannon_slope", bc_to_adult_slope="bc_to_adult_slope",
            aitchison_to_adult_slope="aitchison_slope")
w6_overlay <- NULL
if (!is.null(w6_h2)) {
  w6_h2$metric_h2 <- w6_map[w6_h2$metric]
  w6_overlay <- merge(w6_h2[, c("metric_h2", "h2_w6", "h2_w6_se")],
                      h2[, c("metric", "label", "class_label")],
                      by.x="metric_h2", by.y="metric")
  w6_overlay$label <- factor(w6_overlay$label, levels=levels(h2$label))
  w6_overlay <- w6_overlay[!is.na(w6_overlay$label), ]
}

# Panel A — forest plot of h² ± SE, paired with W6+ refit overlay.
# Per-trait significance is shown by fill color of each marker:
#   filled circle  = full-window h², fill = full-window permutation p-value tier
#   filled diamond = W6+ refit h²,   fill = W6+    permutation p-value tier
full_color  <- "#234a72"   # deep navy-blue (Sidak-significant)
w6_border   <- "#3d3d3d"   # dark grey border for the W6+ diamonds
# Šidák at the M_eff = 9 effective independent tests reported in Methods.
sidak_alpha    <- 1 - (1 - 0.05)^(1 / 9)
sidak_alpha_w6 <- sidak_alpha

# Map metric names in the h² CSV to the full-window permutation CSV.
# (The perm panel uses "t70_imp" / "t75_imp" / "aitchison_slope"; h² CSV
# reports the same threshold metrics from the threshold-specific fit appended
# into h2 under "t70" / "t75" / "aitchison_to_adult_slope".)
perm_lookup <- setNames(perm$p_h2_one_tail, perm$trait)
perm_lookup["t70"] <- perm_lookup["t70_imp"]
perm_lookup["t75"] <- perm_lookup["t75_imp"]
perm_lookup["aitchison_to_adult_slope"] <- perm_lookup["aitchison_slope"]
h2$p_perm_A <- perm_lookup[h2$metric]
h2$perm_tier <- ifelse(is.na(h2$p_perm_A), "untested",
                ifelse(h2$p_perm_A <= sidak_alpha, "sidak",
                ifelse(h2$p_perm_A <= 0.05, "nominal",
                ifelse(h2$p_perm_A <= 0.10, "suggestive", "ns"))))

# W6+ permutation lookup. CSV uses "_w6" suffix on every name.
perm_lookup_w6 <- setNames(perm_w6$p_h2_one_tail,
                           sub("_w6$", "", perm_w6$trait))
# w6_map renames aitchison_to_adult_slope -> aitchison_slope, so the overlay is
# keyed on the h2-CSV name while this lookup is keyed on the perm-CSV name.
perm_lookup_w6["aitchison_slope"] <- perm_lookup_w6["aitchison_to_adult_slope"]
# Attach W6+ permutation p-values and tier to the overlay rows.
if (!is.null(w6_overlay) && nrow(w6_overlay) > 0) {
  w6_overlay$p_perm_A_w6 <- perm_lookup_w6[w6_overlay$metric_h2]
  w6_overlay$perm_tier_w6 <- ifelse(is.na(w6_overlay$p_perm_A_w6), "untested",
                             ifelse(w6_overlay$p_perm_A_w6 <= sidak_alpha_w6, "sidak",
                             ifelse(w6_overlay$p_perm_A_w6 <= 0.05, "nominal",
                             ifelse(w6_overlay$p_perm_A_w6 <= 0.10, "suggestive", "ns"))))
}
tier_fill <- c(sidak      = "#234a72",  # deep navy — strongest signal
               nominal    = "#5b8fbc",  # muted mid-blue
               suggestive = "#9bbdd9",  # lighter blue (0.05 < p ≤ 0.10)
               ns         = "#cfd8e3",  # very pale blue-grey
               untested   = "#e8e8e8")

# Class-specific shading. Each class gets a soft pastel band that runs across
# the panel AND fills the strip on the left, so the eye groups rows by class.
class_colors <- c(
  "FPCA modes"        = "#ece7f4",  # lavender
  "Diversity"         = "#e6f1e6",  # mint
  "Stability"         = "#fbf3e0",  # cream
  "Time-to-threshold" = "#f5e6e6",  # rose
  "Distance to adult" = "#e6eff5",  # ice blue
  "Age prediction"    = "#fbece0")  # peach


# Strip-background colors in the same facet order as the y-axis (top→bottom).
# levels(h2$class_label) is reversed by construction; match strip order to that.
facet_order   <- levels(h2$class_label)
strip_themed_y <- ggh4x::strip_themed(
  background_y = lapply(class_colors[facet_order], function(cc)
    element_rect(fill = cc, color = NA)))

# Connector lines linking full-window point to W6+ point per trait
conn <- NULL
if (!is.null(w6_overlay) && nrow(w6_overlay) > 0) {
  conn <- merge(h2[, c("label","h2","class_label")],
                w6_overlay[, c("label","h2_w6","h2_w6_se")], by="label")
  # merge() returns alphabetically by key; re-order to match h2 row order so
  # the y-axis is not re-derived from conn's row order under ggh4x facet_grid2.
  conn$label <- factor(conn$label, levels=levels(h2$label))
  conn$class_label <- factor(conn$class_label, levels=levels(h2$class_label))
  conn <- conn[order(conn$label, decreasing=TRUE), ]
}

f3a <- ggplot(h2, aes(y=label, x=h2)) +
  # Class grouping is provided by the strip backgrounds on the left
  # (strip_themed_y); the per-panel pastel rectangles were redundant.
  geom_vline(xintercept=0, color="grey75", linewidth=0.35)

# W6+ → full connector — disabled (geom_segment(data=conn) with aes(y=label,
# yend=label, ...) triggers a y-axis reordering bug under ggh4x facet_grid2
# space="free_y", placing the W6+ diamond next to its full-window circle is
# sufficient visual linkage without the connector line).
# if (!is.null(conn) && nrow(conn) > 0) {
#   f3a <- f3a +
#     geom_segment(data=conn,
#                  aes(y=label, yend=label, x=h2, xend=h2_w6),
#                  inherit.aes=FALSE, color="grey55",
#                  linewidth=0.5, alpha=0.6)
# }

f3a <- f3a +
  # SE whiskers — semi-transparent so the marker (focal point) stays dominant
  geom_errorbar(aes(xmin=pmax(0, h2 - h2_se), xmax=h2 + h2_se),
                width=0.32, color=full_color, linewidth=0.8, alpha=0.55,
                orientation="y") +
  geom_point(aes(x=h2, fill=perm_tier), shape=21, size=4.0,
             color="white", stroke=0.7) +
  scale_fill_manual(values=tier_fill,
                    breaks=c("sidak", "nominal", "suggestive", "ns"),
                    labels=c("Sidak (q < 0.05)",
                             "BH FDR q < 0.10",
                             "BH FDR q < 0.15",
                             "n.s."),
                    name="Full window")

if (!is.null(w6_overlay) && nrow(w6_overlay) > 0) {
  f3a <- f3a +
    # Reset the fill scale so the W6+ diamonds get their own legend, separate
    # from the full-window circle legend above. Both palettes share the same
    # three tier colours; what differs is which permutation panel they refer to.
    ggnewscale::new_scale_fill() +
    geom_point(data=w6_overlay, aes(y=label, x=h2_w6, fill=perm_tier_w6),
               inherit.aes=FALSE,
               shape=23, size=3.6, color=w6_border, stroke=0.9) +
    scale_fill_manual(values=tier_fill,
                      breaks=c("sidak", "nominal", "suggestive", "ns"),
                      labels=c("Sidak (q < 0.05)",
                               "BH FDR q < 0.10",
                               "BH FDR q < 0.15",
                               "n.s."),
                      name="W6–W24 refit")
}

f3a <- f3a +
  ggh4x::facet_grid2(class_label ~ ., scales="free_y", space="free_y",
                     switch="y", strip = strip_themed_y) +
  scale_x_continuous(limits=c(-0.02, 0.80),
                     breaks=seq(0, 0.8, 0.1), expand=expansion(add=0.01),
                     oob=scales::oob_keep) +
  labs(y="", x=expression(h^2 ~ "(animal model, REML)"),
       title="Per-trait heritability",
       subtitle="Circles = full window ± 1 SE; diamonds = W6–W24 balanced refit. Marker fill encodes the permutation p-value tier in that window. h²=0 boundary metrics: Table S2.") +
  theme_paper2() +
  theme(strip.placement="outside",
        strip.text.y.left=element_text(angle=0, hjust=1, face="bold",
                                       size=10.5, color="grey10",
                                       margin=margin(r=10, l=4)),
        panel.background=element_rect(fill=NA, color=NA),
        panel.border=element_blank(),
        panel.spacing=unit(1.0, "lines"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(color="white", linewidth=0.6),
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=9.5, color="grey15"),
        axis.ticks.y=element_blank(),
        plot.subtitle=element_text(size=8.5, color="grey25",
                                   margin=margin(b=10)))

# Panel B — obs-vs-W6+ h² scatter (the sensitivity test the paper relies on)
if (!is.null(w6_overlay) && nrow(w6_overlay) > 0) {
  ow <- merge(h2[, c("metric","label","h2","h2_se")],
              w6_overlay[, c("label","h2_w6","h2_w6_se")], by="label")
  ow$delta <- ow$h2_w6 - ow$h2
  ow$direction <- ifelse(ow$delta > 0.03, "rises (W6+ > full)",
                  ifelse(ow$delta < -0.03, "falls (W6+ < full)", "stable"))
  ow$direction <- factor(ow$direction,
                         levels=c("rises (W6+ > full)", "stable", "falls (W6+ < full)"))
  dir_colors <- c("rises (W6+ > full)" = pal$TOL,
                  "stable"             = "grey55",
                  "falls (W6+ < full)" = pal$SEN)

  f3b <- ggplot(ow, aes(x=h2, y=h2_w6, color=direction)) +
    geom_abline(slope=1, intercept=0, lty=3, color="grey50") +
    geom_segment(aes(x=h2 - h2_se, xend=h2 + h2_se,
                     y=h2_w6, yend=h2_w6), alpha=0.45, linewidth=0.5) +
    geom_segment(aes(x=h2, xend=h2,
                     y=pmax(0, h2_w6 - h2_w6_se),
                     yend=h2_w6 + h2_w6_se), alpha=0.45, linewidth=0.5) +
    geom_point(size=3.2, alpha=0.9) +
    geom_text_repel(aes(label=label), size=2.9, color="grey15",
                    seed=42, box.padding=0.4, point.padding=0.3,
                    min.segment.length=0.2, segment.color="grey60",
                    segment.size=0.3, max.overlaps=12) +
    scale_color_manual(values=dir_colors, name="") +
    coord_cartesian(xlim=c(-0.05, 0.7), ylim=c(-0.05, 0.7)) +
    labs(x=expression(h^2~"(full window)"),
         y=expression(h^2~"(W6+ refit)"),
         title="W6+ sensitivity — which traits survive the balanced refit",
         subtitle="Points above 1:1 line strengthen on W6+; below, attenuate") +
    theme_paper2() +
    theme(legend.position="bottom",
          legend.direction="horizontal",
          plot.subtitle=element_text(size=8.5, margin=margin(b=8)))
} else {
  f3b <- ggplot() + annotate("text", x=1, y=1, label="(W6+ CSV not found)") +
    theme_void()
}

# Shared header theme for visual consistency with Fig 2
fig3_header_theme <- theme(
  plot.title    = element_text(face="bold", size=11, hjust=0,
                               margin=margin(t=0, b=4)),
  plot.subtitle = element_text(size=8.5, hjust=0,
                               margin=margin(t=0, b=8)),
  plot.title.position = "plot",
  plot.margin   = margin(t=6, r=10, b=6, l=10))

fig3 <- (f3a + f3b + plot_layout(widths=c(1.7, 1)) +
         plot_annotation(tag_levels = "A")) & fig3_header_theme
ggsave("paper/fig3_h2_per_trait.png", fig3, width=15, height=9, dpi=300)

# FIG 4 — h² by age window + genetic correlations (Obj 2B)
#   Panel A: h² of the per-animal descriptors in the early (W3–W12) and late
#            (W15–W24) windows (13_age_window_h2.R).
#   Panel B: Calo r_EBV heatmap among per-animal trajectory traits, ordered
#            by hierarchical clustering on 1 − |r|.
cat("Fig 4...\n")

win_lab <- c(maz_mean = "MAZ mean", maz_slope = "MAZ slope",
             bc_to_adult_mean = "BC mean", bc_to_adult_slope = "BC slope",
             aitchison_to_adult_slope = "Aitchison slope",
             volatility_mean = "Volatility", volatility_slope = "Volatility slope",
             intra_bc_var = "Intra-BC var", shannon_cv = "Shannon CV",
             shannon_mean = "Shannon mean", shannon_slope = "Shannon slope")
win_cls <- c(maz_mean = "Age prediction", maz_slope = "Age prediction",
             bc_to_adult_mean = "Distance to adult", bc_to_adult_slope = "Distance to adult",
             aitchison_to_adult_slope = "Distance to adult",
             volatility_mean = "Stability", volatility_slope = "Stability",
             intra_bc_var = "Stability", shannon_cv = "Stability",
             shannon_mean = "Diversity", shannon_slope = "Diversity")
win_cls_col <- c("Diversity" = "#3f7a3f", "Stability" = "#a07000",
                 "Distance to adult" = "#2e5b85", "Age prediction" = "#b85c20")
win_levels <- c(sprintf("Weeks 3–12 (n = %d)", win_n[["early"]]),
                sprintf("Weeks 15–24 (n = %d)", win_n[["late"]]))
ord <- names(win_lab)
w <- win_h2
w$label  <- factor(win_lab[w$trait], levels = rev(win_lab[ord]))
w$window <- factor(w$window, levels = c("W3-W12", "W15-W24"), labels = win_levels)
w$se  <- ifelse(is.na(w$h2_se) | w$h2 < 1e-4, 0, w$h2_se)
w$lo  <- pmax(0, w$h2 - w$se); w$hi <- pmin(1, w$h2 + w$se)
w$sig <- w$p_perm <= 0.05
w$sidak <- w$p_perm <= 1 - 0.95^(1/9)
w$fillc <- ifelse(w$sig, as.character(w$window), "none")
w$dgrp <- factor(w$window, levels = rev(levels(w$window)))   # early window above late

win_cols <- setNames(c(pal$h2, pal$c2), win_levels)
y_cols <- win_cls_col[win_cls[ord]][rev(seq_along(ord))]
dodge  <- position_dodge(width = 0.55)

f4a <- ggplot(w, aes(x = h2, y = label, group = dgrp)) +
  geom_errorbar(aes(xmin = lo, xmax = hi, color = window), width = 0,
                orientation = "y", linewidth = 0.7, position = dodge) +
  geom_point(aes(shape = window, fill = fillc, color = window), size = 3.1,
             stroke = 0.9, position = dodge) +
  geom_text(data = transform(w[w$sidak, ],
                             yy = as.numeric(label) + ifelse(window == win_levels[1], 0.1375, -0.1375)),
            aes(x = hi + 0.015, y = yy, label = "**"), inherit.aes = FALSE,
            size = 3.6, color = "grey15", vjust = 0.75, hjust = 0) +
  scale_color_manual(values = win_cols, name = NULL) +
  scale_shape_manual(values = setNames(c(21, 24), win_levels), name = NULL) +
  scale_fill_manual(values = c(win_cols, none = "white"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1), expand = expansion(add = c(0.005, 0))) +
  coord_cartesian(xlim = c(0, 0.66)) +
  labs(x = expression(h^2 ~ "(± 1 SE)"), y = NULL,
       title = "Heritability by age window",
       subtitle = "Filled: permutation p ≤ 0.05; open: p > 0.05; ** p ≤ 0.0057 (family-wise threshold).") +
  theme_paper2() +
  theme(axis.text.y = element_text(size = 10, face = "bold", color = y_cols),
        legend.position = c(0.70, 0.10),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.height = unit(0.9, "lines")) +
  guides(color = guide_legend(override.aes = list(fill = unname(win_cols), shape = c(21, 24), linetype = 0)),
         shape = "none")

# Panel B — Calo r_EBV heatmap, ordered by hierarchical clustering on 1 − |r|
sire_traits <- c("maz_mean","maz_slope","bc_to_adult_mean","bc_to_adult_slope",
                 "aitchison_slope","rel_maturity_w24","volatility_mean",
                 "intra_bc_var","shannon_cv","shannon_slope","BC_FPC1","BC_FPC2")
sire_lab <- c(maz_mean="MAZ mean", maz_slope="MAZ slope",
              bc_to_adult_mean="BC mean", bc_to_adult_slope="BC slope",
              aitchison_slope="Aitchison slope", rel_maturity_w24="Rel. mat. W24",
              volatility_mean="Volatility", intra_bc_var="Intra-BC var",
              shannon_cv="Shannon CV", shannon_slope="Shannon slope",
              BC_FPC1="BC FPC1", BC_FPC2="BC FPC2")
sire_mat <- as.matrix(sire_blu[, sire_traits])
rownames(sire_mat) <- sire_blu$sire_id
r_ebv <- cor(sire_mat, use="pairwise.complete.obs", method="pearson")
rownames(r_ebv) <- sire_lab[sire_traits]
colnames(r_ebv) <- sire_lab[sire_traits]

# Hierarchical clustering ordering: distance = 1 − |r|, complete linkage
hc <- hclust(as.dist(1 - abs(r_ebv)), method = "complete")
trait_order <- hc$labels[hc$order]
n_traits <- length(trait_order)

r_long <- as.data.frame(as.table(r_ebv))
names(r_long) <- c("V1", "V2", "r")
# Drop the first trait from y-axis (its row has zero lower-triangle cells)
# and the last trait from x-axis (its column has zero lower-triangle cells).
# This removes the empty axis stubs that otherwise pad the matrix.
y_levels <- rev(trait_order)[-n_traits]    # drops trait_order[1] from y
x_levels <- trait_order[-n_traits]          # drops trait_order[n] from x
r_long$V1 <- factor(r_long$V1, levels=y_levels)
r_long$V2 <- factor(r_long$V2, levels=x_levels)
r_long    <- r_long[!is.na(r_long$V1) & !is.na(r_long$V2), ]

# Geometric coordinates: y_pos counts from bottom (1=bottom), x_pos from left.
# Trimming the axes (dropping BC slope from y_levels and BC mean from x_levels)
# removes one row and one column but does NOT reindex the remaining positions —
# so the diagonal condition stays at y_pos + x_pos = n_traits + 1 = 13.
y_pos <- match(as.character(r_long$V1), y_levels)
x_pos <- match(as.character(r_long$V2), x_levels)
r_long$is_diag  <- (y_pos + x_pos) == (n_traits + 1)
r_long$is_lower <- (y_pos + x_pos) <  (n_traits + 1)
r_long$abs_r    <- abs(r_long$r)
r_lower <- r_long[r_long$is_lower, ]

# Strong correlations: white label INSIDE the dot. Threshold chosen so the
# label fits comfortably within the dot diameter at max_size=13.
r_strong <- r_lower[r_lower$abs_r >= 0.55, ]
r_strong$label_txt <- sprintf("%+.2f", r_strong$r)

# Map sire-corrgram traits to the same class scheme used in Fig 3,
# then darker-saturated colors so labels remain readable as text.
trait_class_map <- c(
  "maz_mean"="Age prediction", "maz_slope"="Age prediction",
  "bc_to_adult_mean"="Distance to adult", "bc_to_adult_slope"="Distance to adult",
  "aitchison_slope"="Distance to adult", "rel_maturity_w24"="Distance to adult",
  "volatility_mean"="Stability", "intra_bc_var"="Stability", "shannon_cv"="Stability",
  "shannon_slope"="Diversity",
  "BC_FPC1"="FPCA modes", "BC_FPC2"="FPCA modes")
class_text_colors <- c(
  "FPCA modes"="#5e4789", "Diversity"="#3f7a3f", "Stability"="#a07000",
  "Distance to adult"="#2e5b85", "Age prediction"="#b85c20")

# Inverse: human-readable label → trait class → darker hex color
label_to_class <- setNames(trait_class_map[sire_traits], sire_lab[sire_traits])
label_to_color <- setNames(class_text_colors[label_to_class], names(label_to_class))
x_label_colors <- label_to_color[x_levels]   # left → right
y_label_colors <- label_to_color[y_levels]   # bottom → top

f4b <- ggplot(r_lower, aes(x=V2, y=V1)) +
  geom_point(aes(size=abs_r, fill=r),
             shape=21, color="grey30", stroke=0.35) +
  scale_size_area(max_size=13, limits=c(0, 1),
                  breaks=c(0.2, 0.5, 0.8),
                  name="|r|",
                  guide=guide_legend(override.aes=list(fill="grey85", color="grey30"))) +
  scale_fill_gradient2(low=pal$SEN, mid="white", high=pal$TOL,
                       midpoint=0, limits=c(-1, 1), name=expression(r[EBV])) +
  labs(x="", y="",
       title="Calo r_EBV — sire-product genetic correlation",
       subtitle="Lower triangle; dot area scales with |r|, fill = sign; labels for |r| ≥ 0.55. Axis-label colors mark trait class (see Fig 3).") +
  theme_paper2() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold",
                                 color=x_label_colors),
        axis.text.y=element_text(size=10, face="bold",
                                 color=y_label_colors),
        panel.grid.major=element_line(color="grey92", linewidth=0.3),
        panel.grid.minor=element_blank(),
        legend.position="right",
        legend.box="vertical",
        legend.spacing.y=unit(0.4, "lines"))

if (nrow(r_strong) > 0) {
  f4b <- f4b + geom_text(data=r_strong, aes(label=label_txt),
                         size=2.4, color="white", fontface="bold")
}

fig4_header_theme <- theme(
  plot.title    = element_text(face="bold", size=11, hjust=0,
                               margin=margin(t=0, b=4)),
  plot.subtitle = element_text(size=8.5, color="grey25", hjust=0,
                               margin=margin(t=0, b=8)),
  plot.title.position = "plot",
  plot.margin   = margin(t=6, r=10, b=6, l=10))

fig4 <- (f4a + f4b + plot_layout(widths=c(1, 1.35)) +
         plot_annotation(tag_levels = "A")) & fig4_header_theme
ggsave("paper/fig4_window_h2_and_correlations.png", fig4, width=15.5, height=8, dpi=300)

# FIG 6 — Cross-level feature heritability (Obj 3)
cat("Fig 6...\n")

mag_h2_df <- data.frame(h2=mag_h2$h2[!is.na(mag_h2$h2) & mag_h2$converged],
                        level="MAG (n=633)")
ko_h2_df  <- data.frame(h2=ko_h2$h2[!is.na(ko_h2$h2)],
                        level=sprintf("KO (n=%s)",
                                      format(sum(!is.na(ko_h2$h2)), big.mark=",")))
mod_h2_df <- data.frame(h2=mod_h2$h2[!is.na(mod_h2$h2)],
                        level=sprintf("KEGG module (n=%s)",
                                      format(sum(!is.na(mod_h2$h2)), big.mark=",")))
levels_lab <- c(unique(mag_h2_df$level), unique(ko_h2_df$level), unique(mod_h2_df$level))

mag_med <- median(mag_h2_df$h2); ko_med <- median(ko_h2_df$h2); mod_med <- median(mod_h2_df$h2)
mag_frac <- mean(mag_h2_df$h2 > 0.2); ko_frac <- mean(ko_h2_df$h2 > 0.2)
mod_frac <- mean(mod_h2_df$h2 > 0.2)
# Share of each panel pinned at the REML variance-component boundary h2 = 0.
mag_zero <- mean(mag_h2_df$h2 < 1e-6); ko_zero <- mean(ko_h2_df$h2 < 1e-6)
mod_zero <- mean(mod_h2_df$h2 < 1e-6)

cat(sprintf("  per-MAG h²: median %.3f, frac h²>0.2 = %.1f%%, at boundary %.1f%%\n",
            mag_med, 100*mag_frac, 100*mag_zero))
cat(sprintf("  per-KO  h²: median %.3f, frac h²>0.2 = %.1f%%, at boundary %.1f%%\n",
            ko_med, 100*ko_frac, 100*ko_zero))
cat(sprintf("  per-module h²: median %.3f, frac h²>0.2 = %.1f%%, at boundary %.1f%%\n",
            mod_med, 100*mod_frac, 100*mod_zero))

all_levels <- rbind(mag_h2_df, ko_h2_df, mod_h2_df)
all_levels$level <- factor(all_levels$level, levels = levels_lab)
level_colors <- setNames(c("#4F7A55", "#8A4A62", "#3F3A72"), levels_lab)

# Panel A — empirical cumulative distribution.
# Roughly a third of the MAG and KO panels sit at exactly h2 = 0 (REML
# variance-component boundary). A kernel density cannot represent a point mass:
# it smears the atom into a smooth hump whose width is a bandwidth artefact.
# The ECDF shows it as the true vertical rise at x = 0, and makes both reported
# statistics directly legible — the median at y = 0.5, the fraction above 0.2
# as 1 - F(0.2).
f6a <- ggplot(all_levels, aes(x = h2, color = level)) +
  geom_vline(xintercept = 0.2, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  stat_ecdf(geom = "step", linewidth = 0.7, pad = FALSE) +
  annotate("text", x = 0.215, y = 0.03, label = "h² = 0.2", hjust = 0,
           size = 3.0, color = "grey35") +
  annotate("text", x = 0.99, y = 0.535, label = "median", hjust = 1,
           size = 3.0, color = "grey35") +
  scale_color_manual(values = level_colors, name = "Aggregation level") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.025, 0.01))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.015, 0.04))) +
  labs(x = expression(h^2), y = "Cumulative fraction of features",
       title = "Distribution of per-feature heritability",
       subtitle = sprintf(paste0(
         "Median h² (fraction above 0.2): MAG %.3f (%.1f%%), KO %.3f (%.1f%%), ",
         "KEGG module %.3f (%.1f%%).\nThe rise at h² = 0 is the variance-component ",
         "boundary mass: %.1f%% (MAG), %.1f%% (KO), %.1f%% (module)."),
         mag_med, 100*mag_frac, ko_med, 100*ko_frac, mod_med, 100*mod_frac,
         100*mag_zero, 100*ko_zero, 100*mod_zero)) +
  theme_paper2() +
  theme(legend.position   = "right",
        legend.title      = element_text(face = "bold", size = 9.5, hjust = 0),
        legend.text       = element_text(size = 9),
        legend.key.height = unit(1.1, "lines"),
        legend.background = element_blank())

# Panel B — aggregation specificity. Pooling taxa up the GTDB ranks to the same
# breadth does not reproduce the module elevation. NULL_H2 is the module h2
# expected if modules were random KO sets: the variance-weighted mean per-KO h2
# taken over the KOs that actually constitute the module panel.
ko_map   <- read.csv("data/processed/ko_module_map.csv", stringsAsFactors = FALSE)
ko_const <- ko_h2[ko_h2$feature %in%
                    unique(ko_map$ko[ko_map$module %in% mod_h2$module]), ]
NULL_H2  <- sum(ko_const$sigma2_a, na.rm = TRUE) /
              sum(ko_const$sigma2_total, na.rm = TRUE)

lad <- read.csv("results/metrics/taxonomic_ladder_h2.csv", stringsAsFactors = FALSE)
lad_sum <- lad %>% filter(!is.na(h2)) %>%
  group_by(rank) %>%
  summarise(m = median(h2), n = dplyr::n(), sz = mean(n_mags), .groups = "drop")
gv <- function(r, what) lad_sum[[what]][lad_sum$rank == r]
mod_sz <- mean(mod_h2$n_kos_in_module[!is.na(mod_h2$h2)])

agg_df <- data.frame(
  key = c("MAG", "genus", "family", "phylum", "KO", "module"),
  grp = c(rep("Taxonomic pooling", 4), rep("Functional pooling", 2)),
  m   = c(gv("MAG","m"), gv("genus","m"), gv("family","m"), gv("phylum","m"),
          ko_med, mod_med),
  stringsAsFactors = FALSE)
# Row labels carry the pooling breadth, so the reader can see which taxonomic
# ranks are breadth-matched to KEGG modules and which are far coarser.
agg_df$lab <- c(
  sprintf("MAG\n%s features", format(as.integer(gv("MAG","n")), big.mark=",")),
  sprintf("Genus\n%d groups · %.1f MAGs", as.integer(gv("genus","n")), gv("genus","sz")),
  sprintf("Family\n%d groups · %.1f MAGs", as.integer(gv("family","n")), gv("family","sz")),
  sprintf("Phylum\n%d groups · %.1f MAGs", as.integer(gv("phylum","n")), gv("phylum","sz")),
  sprintf("KO\n%s features", format(nrow(ko_h2_df), big.mark=",")),
  sprintf("KEGG module\n%d groups · %.1f KOs", nrow(mod_h2_df), mod_sz))
agg_df$lab <- factor(agg_df$lab, levels = rev(agg_df$lab))
agg_df$grp <- factor(agg_df$grp, levels = c("Taxonomic pooling", "Functional pooling"))
agg_cols <- setNames(c("#4F7A55", "#9A8F82", "#9A8F82", "#9A8F82", "#8A4A62", "#3F3A72"),
                     levels(agg_df$lab)[match(agg_df$lab, levels(agg_df$lab))])
agg_cols <- setNames(c("#4F7A55","#9A8F82","#9A8F82","#9A8F82","#8A4A62","#3F3A72"),
                     as.character(agg_df$lab))

f6b <- ggplot(agg_df, aes(x = m, y = lab, color = lab)) +
  geom_vline(xintercept = NULL_H2, linetype = "dashed",
             color = "grey55", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = m, yend = lab), linewidth = 0.4, color = "grey80") +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", m)), hjust = -0.42, size = 3.0,
            fontface = "bold", show.legend = FALSE) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(values = agg_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 0.46), breaks = seq(0, 0.4, 0.1),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = expression(median~h^2), y = NULL,
       title = "Aggregation specificity",
       subtitle = sprintf(paste0("Median h² by pooling level, with the number of source ",
                                 "features pooled into each group.\nDashed line: the value ",
                                 "expected if modules were random KO sets (%.2f)."), NULL_H2)) +
  theme_paper2() +
  theme(strip.placement    = "outside",
        strip.text.y.left  = element_text(angle = 90, face = "bold", size = 9),
        axis.text.y        = element_text(size = 8.2, lineheight = 1.05),
        panel.grid.major.y = element_blank())

fig6 <- (f6a | f6b) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13)))

ggsave("paper/fig6_cross_level_h2.png", fig6, width = 13, height = 5.6, dpi = 300)

# Supplementary outputs (demoted from main per 2026-05-15 reorg)
cat("Supplementary figs...\n")

# S1 — W3 attenuation: line tests under full vs W6+ window
# Dumbbell forest plot of -log10(p) for SEN vs TOL contrasts under the
# full schedule (W3–W24) and the W6+ refit (W6–W24). Arrows show the
# direction of attenuation. Establishes that several line effects are
# driven by the structural W3 missingness in SEN.
cat("Supplementary S1 (W3 attenuation line tests)...\n")
mat_lt <- read.csv("results/metrics/maturation_w3_artifact_line_tests.csv",
                   stringsAsFactors = FALSE)
alpha_lt <- read.csv("results/metrics/alpha_diversity_w3_artifact_line_tests.csv",
                     stringsAsFactors = FALSE)
mat_lt$family   <- "Maturation / observed metrics"
alpha_lt$family <- "Alpha diversity"
common <- c("metric", "p_full", "p_w6", "family")
s3_df <- rbind(mat_lt[, common], alpha_lt[, common])

s3_labels <- c(
  maz_slope = "MAZ slope",
  bc_to_adult_slope = "BC slope",
  aitchison_to_adult_slope = "Aitchison slope",
  volatility_mean = "Volatility (mean)",
  volatility_slope = "Volatility slope",
  intra_bc_var = "Intra-BC variance",
  shannon_slope = "Shannon slope",
  shannon_cv = "Shannon CV",
  shannon_mean = "Shannon mean",
  shannon_w24 = "Shannon W24",
  simpson_mean = "Simpson mean",
  simpson_slope = "Simpson slope",
  simpson_cv = "Simpson CV",
  simpson_w24 = "Simpson W24",
  invsimpson_mean = "InvSimpson mean",
  invsimpson_slope = "InvSimpson slope",
  invsimpson_cv = "InvSimpson CV",
  invsimpson_w24 = "InvSimpson W24",
  pielou_mean = "Pielou mean",
  pielou_slope = "Pielou slope",
  pielou_cv = "Pielou CV",
  pielou_w24 = "Pielou W24"
)
s3_df$label <- s3_labels[s3_df$metric]
s3_df$label[is.na(s3_df$label)] <- s3_df$metric[is.na(s3_df$label)]

s3_df$mlp_full <- -log10(s3_df$p_full)
s3_df$mlp_w6   <- -log10(s3_df$p_w6)
s3_df$status <- with(s3_df, ifelse(p_full < 0.05 & p_w6 >= 0.05,
                                   "lost on W6+ (W3 artifact)",
                            ifelse(p_full < 0.05 & p_w6 < 0.05,
                                   "robust (sig under both)",
                            ifelse(p_full >= 0.05 & p_w6 < 0.05,
                                   "emerges on W6+",
                                   "n.s. under both"))))
s3_df$status <- factor(s3_df$status,
                       levels = c("robust (sig under both)",
                                  "lost on W6+ (W3 artifact)",
                                  "emerges on W6+",
                                  "n.s. under both"))
status_cols <- c("robust (sig under both)" = pal$SEN,
                 "lost on W6+ (W3 artifact)" = pal$c2,
                 "emerges on W6+" = pal$TOL,
                 "n.s. under both" = "grey65")
s3_df$family <- factor(s3_df$family,
                       levels = c("Maturation / observed metrics",
                                  "Alpha diversity"))
s3_df <- s3_df[order(s3_df$family, s3_df$p_full), ]
# Shannon slope and Shannon CV are tested in both families and carry different
# p-values in each, so the row identity must be family-scoped; a bare label
# factor collapses them to one level positioned by the first family.
row_key <- paste0(as.integer(s3_df$family), "_", s3_df$label)
s3_df$row_id <- factor(row_key, levels = rev(unique(row_key)))

figS1 <- ggplot(s3_df) +
  geom_vline(xintercept = -log10(0.05), lty = 2,
             color = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = -log10(0.10), lty = 3,
             color = "grey75", linewidth = 0.4) +
  geom_segment(aes(x = mlp_full, xend = mlp_w6,
                   y = row_id, yend = row_id, color = status),
               arrow = arrow(length = unit(0.10, "inches"),
                             type = "closed"),
               linewidth = 0.55, alpha = 0.85) +
  geom_point(aes(x = mlp_full, y = row_id, color = status),
             shape = 16, size = 2.4, alpha = 0.9) +
  geom_point(aes(x = mlp_w6, y = row_id, color = status),
             shape = 21, fill = "white", stroke = 1.1, size = 2.4) +
  scale_y_discrete(labels = function(x) sub("^[0-9]+_", "", x)) +
  scale_color_manual(values = status_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = c(0, 1, -log10(0.05), 2, 3),
                     labels = c("0", "1", "1.3", "2", "3")) +
  facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = expression(-log[10](italic(p))),
       y = NULL,
       title = "Line tests attenuate when W3 is removed",
       subtitle = paste("Filled = full-window p; open = W6+ refit p; arrow shows attenuation.",
                        "Dashed = p=0.05, dotted = p=0.10.")) +
  theme_paper2() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, face = "bold"),
        panel.spacing.y = unit(8, "pt"),
        axis.text.y = element_text(size = 8.5))
ggsave("paper/figS1_w3_attenuation_line_tests.png", figS1,
       width = 9, height = 8.5, dpi = 300)

# S2 — Permutation-test null distributions for σ²_A across the metric panel
# Each panel shows the empirical null distribution of ĥ² from B = 1,999
# random permutations of the animal-id → phenotype mapping, with the
# observed ĥ² and one-tailed p-value marked.
cat("Supplementary S2 (permutation null distributions)...\n")
null_dist <- read.csv("results/metrics/maturation_h2_permutation_null_dist.csv",
                      stringsAsFactors = FALSE)
perm_sum  <- read.csv("results/metrics/maturation_h2_permutation_test.csv",
                      stringsAsFactors = FALSE)

# Label map for display + facet ordering by h²_obs (descending)
perm_sum$label <- h2_labels[perm_sum$trait]
perm_sum$label[is.na(perm_sum$label)] <- perm_sum$trait[is.na(perm_sum$label)]
perm_sum <- perm_sum[order(perm_sum$h2_obs, decreasing = TRUE), ]
trait_order  <- perm_sum$trait
label_order  <- perm_sum$label
null_dist$label <- factor(perm_sum$label[match(null_dist$trait, perm_sum$trait)],
                          levels = label_order)
perm_sum$label <- factor(perm_sum$label, levels = label_order)

# Tier coloring matches Fig 3
perm_sum$tier <- ifelse(perm_sum$p_h2_one_tail <= sidak_alpha, "sidak",
                 ifelse(perm_sum$p_h2_one_tail <= 0.05, "nominal",
                 ifelse(perm_sum$p_h2_one_tail <= 0.10, "suggestive", "ns")))

figS2 <- ggplot(null_dist, aes(x = h2_perm_A)) +
  geom_histogram(bins = 30, fill = "grey80", color = "grey55",
                 linewidth = 0.2) +
  geom_vline(data = perm_sum, aes(xintercept = h2_obs, color = tier),
             linewidth = 0.9) +
  geom_text(data = perm_sum,
            aes(x = h2_obs, y = Inf,
                label = sprintf("p = %.3f", p_h2_one_tail)),
            hjust = -0.05, vjust = 1.4, size = 2.7, color = "grey15") +
  scale_color_manual(values = tier_fill,
                     breaks = c("sidak", "nominal", "suggestive", "ns"),
                     labels = c("Sidak (q < 0.05)",
                                "BH FDR q < 0.10",
                                "BH FDR q < 0.15",
                                "n.s."),
                     name = "Permutation tier") +
  facet_wrap(~ label, ncol = 4, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = expression(hat(h)^2 ~ "(permutation null)"),
       y = "Count (B = 1,999 permutations)",
       title = "Empirical null distributions of σ²_A from pedigree-position permutation",
       subtitle = paste("Vertical line = observed ĥ² in each metric, colored by",
                        "one-tailed permutation p-value tier. Panels ordered by",
                        "observed ĥ² (descending).")) +
  theme_paper2() +
  theme(strip.text = element_text(size = 8.5, face = "bold"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        legend.position = "bottom",
        panel.spacing = unit(6, "pt"))

ggsave("paper/figS2_animal_model_perm.png", figS2,
       width = 11, height = 10, dpi = 300)

# S3 — Extended line context across 4 maturation metrics
# Per-line population mean ± SE by week, with W3 (SEN n=2) flagged.
# Supersedes the single-metric figS_sen_tol.
cat("Supplementary S3 (extended SEN/TOL line context)...\n")
maz_in <- read.csv("results/metrics/maz_per_sample.csv", stringsAsFactors = FALSE)
maz_in$animal_id <- as.character(maz_in$animal_id)
maz_in$line <- ped$line[match(maz_in$animal_id, ped$animal_id)]

line_summ <- function(df, val) {
  df$val <- df[[val]]
  df %>% group_by(line, week) %>%
    summarise(mean = mean(val, na.rm = TRUE),
              se   = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
              n    = sum(!is.na(val)),
              .groups = "drop")
}

s3_bc  <- line_summ(curves, "bc_to_adult"); s3_bc$metric  <- "Bray–Curtis to adult"
s3_maz <- line_summ(maz_in, "maz");          s3_maz$metric <- "MAZ"
s3_shn <- line_summ(curves, "shannon");      s3_shn$metric <- "Shannon"

vol_s3 <- curves %>% group_by(animal_id, line) %>% arrange(week) %>%
  mutate(consec_bc = c(NA, abs(diff(bc_to_adult)))) %>%
  filter(!is.na(consec_bc)) %>% ungroup()
vol_s3$week <- vol_s3$week - 1.5
s3_vol <- vol_s3 %>% group_by(line, week) %>%
  summarise(mean = mean(consec_bc, na.rm = TRUE),
            se   = sd(consec_bc, na.rm = TRUE) / sqrt(n()),
            n    = n(), .groups = "drop")
s3_vol$metric <- "Consecutive BC (volatility)"

s3_all <- rbind(s3_bc, s3_maz, s3_shn, s3_vol)
s3_all$metric <- factor(s3_all$metric,
                        levels = c("Bray–Curtis to adult", "MAZ",
                                   "Shannon", "Consecutive BC (volatility)"))

# W3 SEN n flag for annotation
sen_w3 <- s3_all[s3_all$line == "SEN" &
                 abs(s3_all$week - 3) < 0.01 & !is.na(s3_all$n), ]

figS3 <- ggplot(s3_all, aes(x = week, y = mean, color = line, fill = line)) +
  geom_vline(xintercept = 3, lty = 2, color = "grey78", linewidth = 0.35) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.15, alpha = 0.9, lineend = "round") +
  geom_point(shape = 21, fill = "white", stroke = 1.1, size = 2.3) +
  scale_color_manual(values = line_colors, name = NULL) +
  scale_fill_manual(values = line_fills, name = NULL) +
  scale_x_continuous(breaks = c(3, 6, 9, 12, 15, 17, 21, 24)) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(x = "Week", y = "Population mean ± SE",
       title = "SEN vs TOL line context across maturation metrics",
       subtitle = "Vertical dashed line marks W3, where SEN has n=2 — line tests at W3 are W3-artifact-prone (see Fig S1).") +
  theme_paper2() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10))
ggsave("paper/figS3_line_context.png", figS3,
       width = 11, height = 7.5, dpi = 300)

# S4 — Per-litter trajectories for BC-to-adult and MAZ
# Per-litter mean per week (thin lines) over the line population mean
# (bold lines). Lets readers see that litter spread is small relative to
# line-level trajectories, supporting the c² ≈ 0 finding in S4 Panel B.
cat("Supplementary S4 (per-litter trajectories)...\n")
curves$litter_id <- ped$litter_id[match(curves$animal_id, ped$animal_id)]
maz_in$litter_id <- ped$litter_id[match(maz_in$animal_id, ped$animal_id)]

litter_summ <- function(df, val) {
  df$val <- df[[val]]
  df %>% group_by(line, litter_id, week) %>%
    summarise(mean = mean(val, na.rm = TRUE), .groups = "drop")
}
line_summ_simple <- function(df, val) {
  df$val <- df[[val]]
  df %>% group_by(line, week) %>%
    summarise(mean = mean(val, na.rm = TRUE), .groups = "drop")
}

s4_bc_litter <- litter_summ(curves, "bc_to_adult")
s4_bc_line   <- line_summ_simple(curves, "bc_to_adult")
s4_mz_litter <- litter_summ(maz_in, "maz")
s4_mz_line   <- line_summ_simple(maz_in, "maz")

make_litter_panel <- function(litter_df, line_df, y_label, ref_zero = FALSE) {
  p <- ggplot()
  if (ref_zero) {
    p <- p + geom_hline(yintercept = 0, lty = 2,
                        color = "grey78", linewidth = 0.4)
  }
  p +
    geom_line(data = litter_df,
              aes(x = week, y = mean, group = litter_id, color = line),
              linewidth = 0.55, alpha = 0.55, lineend = "round") +
    geom_line(data = line_df,
              aes(x = week, y = mean, color = line),
              linewidth = 1.7, alpha = 0.95, lineend = "round") +
    geom_point(data = line_df,
               aes(x = week, y = mean, color = line),
               shape = 21, fill = "white", stroke = 1.2, size = 2.6) +
    scale_color_manual(values = line_colors, name = NULL) +
    scale_x_continuous(breaks = c(3, 6, 9, 12, 15, 17, 21, 24)) +
    facet_wrap(~ line, ncol = 1, strip.position = "right") +
    labs(x = "Week", y = y_label) +
    theme_paper2() +
    theme(strip.text = element_text(face = "bold"),
          panel.spacing.y = unit(8, "pt"))
}

pS4_A <- make_litter_panel(s4_bc_litter, s4_bc_line,
                           "Bray-Curtis to adult") +
  theme(legend.position = "none")
pS4_B <- make_litter_panel(s4_mz_litter, s4_mz_line,
                           "MAZ", ref_zero = TRUE) +
  theme(legend.position = "none")

figS4 <- (pS4_A | pS4_B) +
  plot_annotation(
    title = "Per-litter trajectories within each line",
    subtitle = "Thin lines: per-litter means. Bold lines + points: line population mean.",
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))) &
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey25"))
ggsave("paper/figS4_per_litter_trajectories.png", figS4,
       width = 12, height = 6.5, dpi = 300)

# S5 — Phenotypic vs r_EBV correlation per trait pair
# Scatter of r_pheno against r_EBV (sire-product proxy for r_g) for the
# trait pairs reported in maturation_genetic_correlations.csv. Pairs that
# include a dropped metric (source_* == "drop") are excluded. None of the
# bivariate REML fits converged, so r_g is shown as the r_EBV proxy.
cat("Supplementary S5 (phenotypic vs genetic correlation)...\n")
rgc <- read.csv("results/metrics/maturation_genetic_correlations.csv",
                stringsAsFactors = FALSE)
rgc <- rgc[rgc$source_1 != "drop" & rgc$source_2 != "drop", ]
rgc$sig_pheno <- rgc$r_pheno_p < 0.05
# Display names, not raw column identifiers (cf. the sire_lab map above).
figS5_lab <- c(maz_mean="MAZ mean", maz_slope="MAZ slope",
               bc_to_adult_mean="BC mean", bc_to_adult_slope="BC slope",
               aitchison_slope="Aitchison slope", rel_maturity_w24="Rel. mat. W24",
               volatility_mean="Volatility", intra_bc_var="Intra-BC var",
               shannon_cv="Shannon CV", shannon_slope="Shannon slope",
               BC_FPC1="BC FPC1", BC_FPC2="BC FPC2", t70="t70", t75="t75")
rgc$pair_label <- vapply(strsplit(rgc$pair, " x ", fixed = TRUE), function(z) {
  nm <- ifelse(z %in% names(figS5_lab), figS5_lab[z], gsub("_", " ", z))
  paste(nm, collapse = "  \u00d7  ")
}, character(1))

figS5 <- ggplot(rgc, aes(x = r_pheno, y = r_EBV)) +
  geom_hline(yintercept = 0, lty = 3, color = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, lty = 3, color = "grey75", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, lty = 2,
              color = "grey70", linewidth = 0.4) +
  geom_point(aes(color = sig_pheno),
             shape = 21, fill = "white", stroke = 1.4, size = 3.2) +
  ggrepel::geom_text_repel(aes(label = pair_label),
                           size = 3, color = "grey20",
                           min.segment.length = 0.2,
                           max.overlaps = 20, seed = 7) +
  scale_color_manual(values = c(`TRUE` = pal$SEN, `FALSE` = "grey55"),
                     labels = c(`TRUE` = "p(r_pheno) < 0.05",
                                `FALSE` = "p(r_pheno) ≥ 0.05"),
                     name = NULL) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  coord_equal() +
  labs(x = "Phenotypic correlation r_pheno",
       y = "Sire-product correlation r_EBV (proxy for r_g)",
       title = "Phenotypic vs sire-product correlation across trait pairs",
       subtitle = sprintf(
         "n = %d pairs (n=77 animals). All bivariate REML fits failed to converge, so r_g is shown as the sire-product proxy. Dashed = 1:1.",
         nrow(rgc))) +
  theme_paper2() +
  theme(legend.position = "bottom")
ggsave("paper/figS5_pheno_vs_genetic_correlation.png", figS5,
       width = 8, height = 8.5, dpi = 300)

# S6 — Cross-estimator concordance of h² (parallel-coordinates view)
# Each metric's h² across four estimators (obs full, PACE-imputed, within-line,
# W6–W24 balanced refit). Restricted to the metrics with values in all four
# panels. The W3-robust metrics (heritable on W6+) are highlighted; the
# remaining metric is shown in grey.
cat("Supplementary S6 (cross-estimator concordance)...\n")
conc_h2  <- read.csv("results/metrics/maturation_h2_concordance.csv",
                     stringsAsFactors = FALSE)
conc_sum <- read.csv("results/metrics/maturation_h2_concordance_summary.csv",
                     stringsAsFactors = FALSE)

# Keep metrics with all four estimators (= the all-panel subset that
# the Kendall's W summary is computed on).
conc_full <- conc_h2[complete.cases(conc_h2[, c("h2_obs","h2_imp",
                                                "h2_within_line","h2_w6")]), ]
w3_robust_set <- c("maz_mean", "intra_bc_var", "volatility_mean",
                   "bc_to_adult_slope", "aitchison_slope", "shannon_slope")
conc_full$w3_robust <- conc_full$metric %in% w3_robust_set
conc_full$label <- h2_labels[conc_full$metric]
conc_full$label[is.na(conc_full$label)] <- conc_full$metric[is.na(conc_full$label)]

# Long format for parallel coordinates
conc_long <- reshape2::melt(conc_full,
                            id.vars = c("metric", "label", "w3_robust"),
                            measure.vars = c("h2_obs", "h2_imp",
                                             "h2_within_line", "h2_w6"),
                            variable.name = "estimator",
                            value.name = "h2")
estimator_levels <- c("h2_obs", "h2_imp", "h2_within_line", "h2_w6")
estimator_labels <- c("Observed\n(full window)",
                      "PACE-imputed",
                      "Within-line\n(line fixed)",
                      "W6–W24\nbalanced")
conc_long$estimator <- factor(conc_long$estimator,
                              levels = estimator_levels,
                              labels = estimator_labels)

W_all <- conc_sum$Kendall_W[conc_sum$panel == "all_metrics"]
p_all <- conc_sum$p_perm_9999[conc_sum$panel == "all_metrics"]
W_w3r <- conc_sum$Kendall_W[conc_sum$panel == "w3_robust"]
p_w3r <- conc_sum$p_perm_9999[conc_sum$panel == "w3_robust"]
n_w3r <- sum(conc_full$w3_robust); n_all <- nrow(conc_full)

# End-of-line labels at the W6–W24 endpoint, with light repulsion to avoid
# overlap when multiple metrics converge.
endpoint_lab <- conc_long[conc_long$estimator == "W6–W24\nbalanced", ]

robust_color <- "#234a72"   # navy — W3-robust
fragile_color <- "#9a9a9a"  # grey — attenuating

figS6 <- ggplot(conc_long, aes(x = estimator, y = h2,
                               group = metric, color = w3_robust)) +
  geom_line(linewidth = 1.0, alpha = 0.9) +
  geom_point(size = 3.2, alpha = 0.95) +
  ggrepel::geom_text_repel(data = endpoint_lab,
                           aes(label = label),
                           hjust = 0, nudge_x = 0.20,
                           direction = "y", segment.size = 0.25,
                           segment.color = "grey55",
                           size = 3.1, show.legend = FALSE,
                           max.overlaps = Inf, seed = 20260518) +
  scale_color_manual(values = c(`TRUE` = robust_color, `FALSE` = fragile_color),
                     labels = c(`TRUE` = sprintf("W3-robust (%d metrics)", n_w3r),
                                `FALSE` = sprintf("Attenuates on W6+ (%d metric%s)", n_all - n_w3r,
                                                  ifelse(n_all - n_w3r == 1, "", "s"))),
                     name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.2, 1.1))) +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.1)) +
  labs(x = NULL,
       y = expression(h^2 ~ "(animal model, REML)"),
       title = "Cross-estimator concordance of per-animal h²",
       subtitle = sprintf(
         "Kendall's W = %.2f on the W3-robust %d-metric subset (perm p = %.3f); W = %.2f on the full %d-metric panel (perm p = %.3f).",
         W_w3r, n_w3r, p_w3r, W_all, n_all, p_all)) +
  theme_paper2() +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "grey15"))

ggsave("paper/figS6_cross_estimator_concordance.png", figS6,
       width = 9, height = 6.5, dpi = 300)

cat("\nAll figures saved to paper/\n")
cat("  fig1_panels_BC.png       (design panels B + C; Fig 1A external)\n")
cat("  fig2_maturation.png      (Obj 1)\n")
cat("  fig3_h2_per_trait.png    (Obj 2A — per-trait h²)\n")
cat("  fig4_window_h2_and_correlations.png  (h² by age window + Calo r_EBV)\n")
cat("  fig6_cross_level_h2.png  (Obj 3 — per-MAG/KO/module distributions)\n")
cat("Supplementary figures:\n")
cat("  figS1_w3_attenuation_line_tests.png (line tests under full vs W6+)\n")
cat("  figS2_animal_model_perm.png    (permutation null distributions)\n")
cat("  figS3_line_context.png         (SEN vs TOL across 4 metrics)\n")
cat("  figS4_per_litter_trajectories.png (per-litter BC + MAZ)\n")
cat("  figS5_pheno_vs_genetic_correlation.png (r_pheno vs r_EBV)\n")
cat("  figS6_cross_estimator_concordance.png  (parallel-coordinates of h² across estimators)\n")
cat("\nFig 5 (sire BVs + Pareto + selection response) → run 17_pareto_selection_figure.R\n")

# PART 2: Fig 5 (sire EBV heatmap + selection response)
setwd(here::here())
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
source("scripts/palette.R")

blups <- read.csv("results/metrics/maturation_sire_blups.csv",
                  stringsAsFactors = FALSE)
sel   <- read.csv("results/metrics/maturation_selection_response.csv",
                  stringsAsFactors = FALSE)
ped   <- read.csv("data/metadata/pedigree.csv", stringsAsFactors = FALSE)

blups$sire_id <- as.character(blups$sire_id)
ped$sire_id   <- as.character(ped$sire_id)

sire_line <- ped %>%
  distinct(sire_id, line) %>%
  group_by(sire_id) %>% slice(1) %>% ungroup()
sire_n <- ped %>% count(sire_id, name = "n_progeny")
sire_meta <- sire_line %>% left_join(sire_n, by = "sire_id")
# Anonymized sire codes S1-S7, keyed on the same arrange(line, sire_id) order
# used for the Fig 1A pedigree panel so the labels agree across figures/tables.
sire_meta <- sire_meta %>%
  arrange(line, sire_id) %>%
  mutate(sire_code = sprintf("S%d", dplyr::row_number()))

# Trait set: heritable traits matching Table S12. Use display labels and a fixed
# trait order so Panel A and Panel B read in the same direction.
trait_keys <- c("maz_mean", "BC_FPC2", "aitchison_slope",
                "bc_to_adult_slope", "shannon_slope", "intra_bc_var",
                "t75", "t70")
trait_labels <- c(maz_mean          = "MAZ mean",
                  BC_FPC2           = "BC FPC2",
                  aitchison_slope   = "Aitchison slope",
                  bc_to_adult_slope = "BC slope",
                  shannon_slope     = "Shannon slope",
                  intra_bc_var      = "Intra-BC var",
                  t75               = "t75",
                  t70               = "t70")

# Panel A — Sire EBV heatmap, in phenotypic-SD units (matches Table 3A)
# Per-trait phenotypic SD from the selection-response table; rescales raw
# BLUPs to the σ_P units used in Table 3A and the §Sire BVs narrative.
sigma_P_lookup <- setNames(sel$sigma_P, sel$metric)

blups_long <- blups %>%
  select(sire_id, all_of(trait_keys)) %>%
  pivot_longer(-sire_id, names_to = "metric", values_to = "blup_raw") %>%
  mutate(blup = blup_raw / sigma_P_lookup[metric]) %>%
  left_join(sire_meta, by = "sire_id") %>%
  mutate(trait = factor(trait_labels[metric], levels = unname(trait_labels)))

# Sire row order: MAZ mean BLUP descending (matches Table 3A).
sire_order <- blups_long %>%
  filter(metric == "maz_mean") %>%
  arrange(desc(blup)) %>%
  pull(sire_id)
blups_long$sire_label <- paste0(sire_meta$sire_code[match(blups_long$sire_id,
                                                          sire_meta$sire_id)],
                                " (", sire_meta$line[match(blups_long$sire_id,
                                                           sire_meta$sire_id)],
                                ", n=",
                                sire_meta$n_progeny[match(blups_long$sire_id,
                                                          sire_meta$sire_id)],
                                ")")
sire_label_order <- unique(blups_long$sire_label[match(sire_order,
                                                       blups_long$sire_id)])
blups_long$sire_label <- factor(blups_long$sire_label, levels = rev(sire_label_order))

# Symmetric color limit so red/blue diverge around 0.
blup_lim <- max(abs(blups_long$blup), na.rm = TRUE)

pA <- ggplot(blups_long, aes(x = trait, y = sire_label, fill = blup)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", blup)),
            size = 3.0, color = "grey15") +
  scale_fill_gradient2(low = pal$SEN, mid = "white", high = pal$TOL,
                       midpoint = 0, limits = c(-blup_lim, blup_lim),
                       name = "Sire BLUP",
                       guide = guide_colorbar(barwidth = unit(0.35, "cm"),
                                              barheight = unit(3.2, "cm"))) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_paper2() +
  theme(axis.text.x = element_text(angle = 35, hjust = 0, face = "bold"),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "right")

# Panel B — Expected per-generation response (ΔG / σ_P) at i = 1.40
sel <- sel %>%
  filter(metric %in% trait_keys) %>%
  mutate(trait = factor(trait_labels[metric], levels = unname(trait_labels)))

resp <- bind_rows(
  sel %>% transmute(trait,
                    regime  = "Mass selection",
                    est     = dG_mass_i14,
                    lo      = dG_mass_i14_lo,
                    hi      = dG_mass_i14_hi),
  sel %>% transmute(trait,
                    regime  = "Sire-BLUP selection",
                    est     = dG_sire_i14,
                    lo      = dG_sire_i14_lo,
                    hi      = dG_sire_i14_hi)
)
resp$regime <- factor(resp$regime,
                      levels = c("Mass selection", "Sire-BLUP selection"))
# Trait order on y-axis: matches Panel A column order (MAZ mean at top).
# ggplot stacks y from bottom up, so reverse the panel-A order.
resp$trait <- factor(resp$trait, levels = rev(unname(trait_labels)))

regime_colors <- c("Mass selection"      = pal$neutral,
                   "Sire-BLUP selection" = pal$h2)

pB <- ggplot(resp, aes(x = est, y = trait, color = regime)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey75", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0, linewidth = 0.7,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.7, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = regime_colors, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = expression("Expected per-generation response  " * Delta * "G /" ~ sigma[P]),
       y = NULL) +
  theme_paper2() +
  theme(legend.position = c(0.72, 0.12),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.height = unit(0.9, "lines"),
        legend.text = element_text(size = 9.5),
        panel.grid.major.y = element_blank())

# Assemble
fig5 <- (pA | pB) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A",
                  caption = paste0(
                    "Selection intensity i = 1.40 (top 20% of candidates). ",
                    "Sire-BLUP responses use the trait-specific PEV-based ",
                    "accuracy r averaged across the seven sires."),
                  theme = theme(plot.tag = element_text(face = "bold", size = 13),
                                plot.caption = element_text(size = 8.5,
                                                            color = "grey30",
                                                            hjust = 0)))

ggsave("paper/fig5_sire_bvs_selection.png", fig5,
       width = 13.5, height = 6.4, dpi = 300)
cat("\nWrote paper/fig5_sire_bvs_selection.png\n")
