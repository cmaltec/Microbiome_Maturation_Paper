# Microbiome maturation heritability — analysis pipeline

Analysis code reproducing:

> Maltecca C, Gini C, Musa J, Byrd MH, Brito LF, Johnson JS, Jiang J, Tiezzi F.
> **Host genetics contributes to microbiome maturation trajectories and pathway-level function in swine.**
> *Microbiome* (submitted). DOI: `[to be added on acceptance]`

Twenty R scripts plus a shared plotting theme. Each script is
self-contained and resolves paths relative to the repository root via
`here::here()`; the `.here` file at the root is the anchor.

---

## Data are not in this repository

No data files are stored here. To run the pipeline you need:

| What | Where | Notes |
|---|---|---|
| Feature matrices, pedigree, sample metadata, annotation | Zenodo, DOI [10.5281/zenodo.21761856](https://doi.org/10.5281/zenodo.21761856) (version 1.1.0 or later) | MAG and KO abundance tables, three-generation pedigree, GTDB taxonomy, CheckM2 quality, KO→module map, per-sample read counts |
| Pipeline intermediates | same Zenodo record, `data/derived/` | MAZ scores, per-sample diversity curves, FPCA scores and the h² panels that scripts 00c, 02, 05, 09, 11, 12 and 13 read as inputs |
| Raw whole-genome shotgun reads | NCBI BioProject PRJNA1504886 | Only needed to regenerate the MAG catalog from scratch |

Unpack the deposit so that its `data/` folder sits at `data/`, then run
`00a_prepare_inputs.R`, which derives `data/metadata/pedigree.csv`, pivots
the long-format abundance tables into the wide matrices the pipeline reads, and
copies the derived inputs in `data/derived/` to `results/metrics/`.

The outputs of scripts 06–15 behind the revised manuscript (random regression,
permutation nulls, age-window heritabilities, genetic groups, sequencing depth;
Tables S5 and S15) are tracked under `results/revision/`.

Sire, litter and grandparent identifiers in the deposit are anonymised
(`S1`–`S7`, `L1`–`L13`, `GS*`/`GD*`). The relabelling preserves the pedigree
structure exactly, so every estimate reproduces.

---

## Scripts, in run order

| # | Script | Produces |
|---:|---|---|
| 00a | `scripts/00a_prepare_inputs.R` | `data/metadata/pedigree.csv` and the wide MAG and KO matrices, from the deposited long-format tables; copies `data/derived/` to `results/metrics/` |
| 00b | `scripts/00b_build_a_matrix.R` | Two-generation (sire + litter-as-dam) **A** and **A⁻¹**, used by the pedigree-depth sensitivity in 05 and by 10 |
| 00c | `scripts/00c_revision_inputs.R` | Three- and two-generation **A⁻¹**, per-sample long table and per-animal trait table read by 06, 06b, 13 and 14 (`results/revision/inputs/`) |
| 01 | `scripts/01_pace_imputation.R` | PACE (fdapace FPCA) imputation of per-sample MAG and Shannon trajectories |
| 02 | `scripts/02_per_animal_phenotypes.R` | The 24 per-animal scalars across six metric classes (Table 2) |
| 03 | `scripts/03_growth_curve.R` | Cohort-mean Gompertz-with-Gaussian-dip fit (Fig 2A, Table S3) |
| 04 | `scripts/04_pedigree.R` | Three-generation pedigree (114 individuals) and ASReml-form **A** and **A⁻¹**, from `data/pedigree_3generation.csv` |
| 05 | `scripts/05_heritability_panel.R` | Univariate animal-model h² panel, permutation tests, W6–W24 balanced refit (Table S2, Fig 3) |
| 06 | `scripts/06_rrm_h2t.R` | Random-regression models for community-level variance by week: additive, permanent-environment and joint fits, residual structure, Legendre order (Table S5A) |
| 06b | `scripts/06b_rrm_permutation.R` | Permutation test of pedigree structure in the random regression (B = 999); Table S5A and S5B |
| 07 | `scripts/07_multivariate_test.R` | Hotelling's *T*² and PERMANOVA on the six-class metric vector (Table S7) |
| 08 | `scripts/08_genetic_correlations.R` | Sire BLUPs and accuracies, Calo *r*<sub>EBV</sub>, bootstrap CIs, Bayesian bivariate fits (Table 3, Table S11) |
| 09 | `scripts/09_selection_response.R` | Expected per-generation response at *i* = 1.40 and 2.06 (Table S12, Fig 5B) |
| 10 | `scripts/10_feature_level_h2.R` | Per-MAG, per-KO and per-KEGG-module h² panels (Tables S9, S13, S14; Fig 6) |
| 11 | `scripts/11_sensitivities.R` | Cross-estimator concordance, age and pre-weaning-growth covariate refits, pre-perturbation refit (Table S10, Methods S.9, S.10, S.11) |
| 12 | `scripts/12_figures.R` | All main figures (1–6) and supplementary figures (S1–S6). Sources `palette.R` |
| 13 | `scripts/13_age_window_h2.R` | Heritability of the per-animal descriptors in early (W3–W12) and late (W15–W24) windows, and late without W21, with permutation tests (Fig 4A, Table S15) |
| 14 | `scripts/14_genetic_groups.R` | Animal model with genetic groups for unknown parents compared with selection line as fixed (Table S2b footnote) |
| 15 | `scripts/15_sequencing_depth.R` | Per-sample read pairs after quality control and reads mapped to the 633 core MAGs (Methods) |
| — | `scripts/palette.R` | Shared colour palette and the `theme_paper2()` ggplot theme |

Run order: 00a, 00b, 04, 00c, then 01–03, 05–11, 13–15, and 12 last. Script 02 reads
the pedigree written by 04. `12_figures.R` writes its output to `paper/`; scripts 00c,
06, 06b and 13–15 write to `results/revision/`.

`paper/paper2_metric_source.csv` is configuration, not output: it records the
per-metric observed-vs-imputed source rule (Additional file 1, §S2) and is read
by scripts 05 and 08.

---

## Requirements

Results in the paper were produced with **R 4.4** and the following package
versions. Later versions will generally work, but these are the versions of
record:

| Package | Version | Used for |
|---|---|---|
| ASReml-R | 4.2 | all variance-component fits |
| MCMCglmm | 2.36 | Bayesian bivariate animal models |
| fdapace | 0.6.0 | PACE imputation and FPCA |
| vegan | 2.7.2 | PERMANOVA |
| picante | 1.8 | phylogenetic diversity |
| KEGGREST | 1.42 | KO → KEGG module mapping |

Also required: `ape`, `arrow`, `coda`, `dplyr`, `ggh4x`, `ggnewscale`,
`ggplot2`, `ggrepel`, `here`, `MASS`, `Matrix`, `minpack.lm`, `patchwork`,
`png`, `reshape2`, `scales`, `tidyr`.

> **ASReml-R is commercial software** and requires a licence from VSN
> International. Scripts **00b, 00c, 02, 04, 05, 06, 06b, 08, 10, 11, 13 and 14** will
> not run without it. Scripts **00a, 01, 03, 07, 09, 12 and 15** do not use it: 01 (PACE)
> and 03 (growth curves) run from the MAG matrix, 15 from the deposited read counts, and
> 07, 09 and 12 from the outputs of the earlier scripts.

---

## Licence

Code is released under the MIT Licence — see `LICENSE`. The associated dataset
on Zenodo is released separately under CC BY 4.0. Please cite the article when
reusing either.

## Contact

Christian Maltecca — cmaltec@ncsu.edu
Department of Animal Science, North Carolina State University, Raleigh, NC, USA
