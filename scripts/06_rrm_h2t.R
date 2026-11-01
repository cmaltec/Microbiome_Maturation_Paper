#!/usr/bin/env Rscript
#________________________________________
# - File Name : 06_rrm_h2t.R
# - Purpose : Random-regression animal models for community-level variance by week (Table S5A, h2(t) upper bound)
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
OUT <- "results/revision/rrm"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

Ainv3 <- readRDS(file.path(IN, "Ainv_3gen.rds"))
Ainv2 <- readRDS(file.path(IN, "Ainv_2gen.rds"))
long  <- read.csv(file.path(IN, "per_sample_long.csv"), stringsAsFactors = FALSE)

WEEKS <- c(3, 6, 9, 12, 15, 17, 21, 24)
WMEAN <- mean(WEEKS)                       # week centring of release 1.0.0
xstd  <- function(w) 2 * (w - min(WEEKS)) / (max(WEEKS) - min(WEEKS)) - 1
basis <- function(type, w) {
  x <- xstd(w)
  switch(type,
    raw  = cbind(1, w - WMEAN),
    leg1 = cbind(sqrt(1/2), sqrt(3/2) * x),
    leg2 = cbind(sqrt(1/2), sqrt(3/2) * x, sqrt(5/2) * (3 * x^2 - 1) / 2),
    int  = matrix(1, length(w), 1))
}

long$week_f    <- factor(long$week)
long$line      <- factor(long$line)
long$litter_id <- factor(long$litter_id)
long$wc <- long$week - WMEAN
Bl <- basis("leg2", long$week)
long$L0 <- Bl[, 1]; long$L1 <- Bl[, 2]; long$L2 <- Bl[, 3]

# Random blocks. A = additive (pedigree), P = permanent environment (iid animal)
# The variance model of a str() block is its second formula (ASReml-R manual,
# complex_vs: "the variance structure defined is the structure required across
# all terms in the set"). The additive block must therefore carry
# us(k):vm(animal_id, Ainv); us(k):id(N) imposes an identity covariance among
# animals. Release 1.0.0 of this script used us(2):id(78); idcov = TRUE
# reproduces it.
blk <- function(eff, type, idcov = FALSE) {
  if (type == "int") return(eff)
  cv <- switch(type, raw = c("1", "wc"), leg1 = c("L0", "L1"), leg2 = c("L0", "L1", "L2"))
  terms <- ifelse(cv == "1", eff, paste0(cv, ":", eff))
  k <- length(cv)
  vmod <- if (grepl("^vm", eff) && !idcov) sprintf("us(%d):%s", k, eff) else
    sprintf("us(%d):id(%s)", k, if (grepl("^vm", eff)) "NA_A" else "NA_P")
  sprintf("str(~ %s, ~ %s)", paste(terms, collapse = " + "), vmod)
}
EA <- "vm(animal_id, Ainv_cur)"; EP <- "pe"

models <- list(
  M0_submitted   = list(ped = "2gen", fixed = "lin",   A = "raw",  P = NA,     res = "hom", idcov = TRUE),
  M0b_pedcov     = list(ped = "2gen", fixed = "lin",   A = "raw",  P = NA,     res = "hom"),
  M0c_ped3cov    = list(ped = "3gen", fixed = "lin",   A = "raw",  P = NA,     res = "hom"),
  M1_Arr_hom     = list(ped = "3gen", fixed = "class", A = "raw",  P = NA,     res = "hom"),
  M2_Arr_het     = list(ped = "3gen", fixed = "class", A = "raw",  P = NA,     res = "het"),
  M2L1_Aleg1_het = list(ped = "3gen", fixed = "class", A = "leg1", P = NA,     res = "het"),
  M2L2_Aleg2_het = list(ped = "3gen", fixed = "class", A = "leg2", P = NA,     res = "het"),
  M3_Prr_het     = list(ped = "3gen", fixed = "class", A = NA,     P = "raw",  res = "het"),
  M4_ArrPrr_het  = list(ped = "3gen", fixed = "class", A = "raw",  P = "raw",  res = "het"),
  M5_AintPrr_het = list(ped = "3gen", fixed = "class", A = "int",  P = "raw",  res = "het"),
  M3h_Prr_hom    = list(ped = "3gen", fixed = "class", A = NA,     P = "raw",  res = "hom"),
  M4h_ArrPrr_hom = list(ped = "3gen", fixed = "class", A = "raw",  P = "raw",  res = "hom"),
  M6_rep_hom     = list(ped = "3gen", fixed = "class", A = "int",  P = "int",  res = "hom"),
  M7_rep_het     = list(ped = "3gen", fixed = "class", A = "int",  P = "int",  res = "het")
)

fit_model <- function(trait, sp) {
  Ainv <- if (sp$ped == "2gen") Ainv2 else Ainv3
  ids  <- attr(Ainv, "rowNames")
  df <- long[!is.na(long[[trait]]), ]
  df$y <- as.numeric(scale(df[[trait]]))
  df$animal_id <- factor(as.character(df$animal_id), levels = ids)
  df$pe <- factor(as.character(df$animal_id))
  df <- df[order(df$week_f, df$animal_id), ]
  assign("Ainv_cur", Ainv, envir = .GlobalEnv)
  idc <- isTRUE(sp$idcov)
  rt <- c(if (!is.na(sp$A)) blk(EA, sp$A, idc), if (!is.na(sp$P)) blk(EP, sp$P), "litter_id")
  rt <- gsub("NA_A", length(ids), gsub("NA_P", nlevels(df$pe), rt))
  fixed <- if (sp$fixed == "lin") y ~ 1 + wc else y ~ 1 + week_f + line
  res   <- if (sp$res == "hom") ~ idv(units) else ~ dsum(~ units | week_f)
  msg <- NA_character_
  fit <- tryCatch(suppressWarnings(suppressMessages(
    asreml(fixed = fixed, random = as.formula(paste("~", paste(rt, collapse = " + "))),
           residual = res, data = df, trace = FALSE, maxit = 300,
           na.action = na.method(y = "include")))),
    error = function(e) { msg <<- gsub("\\s+", " ", conditionMessage(e)); NULL })
  list(fit = fit, msg = msg)
}

# V-index matrix of a random block in varcomp
K_idx <- function(vc, eff_regex, k) {
  rn <- rownames(vc)
  if (k == 1) {
    i <- grep(paste0("^", eff_regex, "$"), rn)
    stopifnot(length(i) == 1); return(matrix(i))
  }
  M <- matrix(NA_integer_, k, k)
  for (i in 1:k) for (j in 1:i) {
    h <- grep(sprintf("%s.*!us\\(%d\\)_%d:%d$", eff_regex, k, i, j), rn)
    stopifnot(length(h) == 1); M[i, j] <- M[j, i] <- h
  }
  M
}
qf <- function(M, phi) {
  k <- length(phi); tt <- c()
  for (i in 1:k) for (j in 1:k) tt <- c(tt, sprintf("(%.12g)*V%d", phi[i] * phi[j], M[i, j]))
  paste0("(", paste(tt, collapse = " + "), ")")
}
ev <- function(expr, comp) eval(parse(text = gsub("V(\\d+)", "comp[\\1]", expr)))

per_week <- function(fit, sp, trait, mname) {
  vc <- summary(fit)$varcomp; comp <- vc$component
  KA <- if (!is.na(sp$A)) K_idx(vc, "vm\\(animal_id, Ainv_cur\\)", ncol(basis(sp$A, 3))) else NULL
  KP <- if (!is.na(sp$P)) K_idx(vc, "pe", ncol(basis(sp$P, 3))) else NULL
  ic <- grep("^litter_id$", rownames(vc))
  ie <- if (sp$res == "hom") rep(grep("^units!units$", rownames(vc)), length(WEEKS)) else
    sapply(WEEKS, function(w) grep(sprintf("^week_f_%d!R$", w), rownames(vc)))
  # index sanity: vpredict SE of a single component equals summary SE
  i1 <- if (!is.null(KA)) KA[1, 1] else KP[1, 1]
  chk <- vpredict(fit, as.formula(sprintf("x ~ V%d", i1)))
  if (!is.na(vc$std.error[i1]))
    stopifnot(isTRUE(all.equal(chk$SE, vc$std.error[i1], tolerance = 1e-6)))
  do.call(rbind, lapply(seq_along(WEEKS), function(t) {
    eA <- if (!is.null(KA)) qf(KA, basis(sp$A, WEEKS[t])[1, ]) else "0"
    eP <- if (!is.null(KP)) qf(KP, basis(sp$P, WEEKS[t])[1, ]) else "0"
    tot <- sprintf("(%s + %s + V%d + V%d)", eA, eP, ic, ie[t])
    vp <- function(num) {
      if (num == "0") return(c(0, NA))
      r <- tryCatch(vpredict(fit, as.formula(paste("x ~", num, "/", tot))),
                    error = function(e) data.frame(Estimate = NA, SE = NA))
      c(r$Estimate, r$SE)
    }
    a <- vp(eA); p <- vp(eP)
    an <- vp(sprintf("(%s + %s)", eA, eP))
    data.frame(trait = trait, model = mname, week = WEEKS[t],
               s2_A = ev(eA, comp), s2_PE = ev(eP, comp),
               s2_c = comp[ic], s2_e = comp[ie[t]],
               h2 = a[1], h2_se = a[2], pe2 = p[1], pe2_se = p[2],
               animal2 = an[1], animal2_se = an[2])
  }))
}

traits <- c(bc_to_adult = "BC-to-adult", shannon = "Shannon", maz = "MAZ")
H <- list(); Fm <- list(); V <- list()
for (tr in names(traits)) for (mn in names(models)) {
  sp <- models[[mn]]
  r <- fit_model(tr, sp); fit <- r$fit
  cat(sprintf("\n=== %s | %s ===\n", traits[tr], mn))
  if (is.null(fit)) {
    cat("  FAILED:", r$msg, "\n")
    Fm[[length(Fm) + 1]] <- data.frame(trait = traits[tr], model = mn, fitted = FALSE,
      converge = NA, logLik = NA, npar = NA, AIC = NA, n_boundary = NA, note = r$msg)
    next
  }
  vc <- summary(fit)$varcomp
  npar <- sum(vc$bound != "F"); nb <- sum(vc$bound %in% c("B", "S"))
  cat(sprintf("  converge=%s logLik=%.4f npar=%d boundary=%d\n", fit$converge, fit$loglik, npar, nb))
  print(vc[, c("component", "std.error", "bound")], digits = 4)
  V[[length(V) + 1]] <- data.frame(trait = traits[tr], model = mn, term = rownames(vc),
                                   vc[, c("component", "std.error", "bound")], row.names = NULL)
  Fm[[length(Fm) + 1]] <- data.frame(trait = traits[tr], model = mn, fitted = TRUE,
    converge = fit$converge, logLik = fit$loglik, npar = npar,
    AIC = -2 * fit$loglik + 2 * npar, n_boundary = nb, note = NA)
  h <- tryCatch(per_week(fit, sp, traits[tr], mn),
                error = function(e) { cat("  per-week error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(h)) {
    print(h[, c("week", "h2", "h2_se", "pe2", "animal2", "animal2_se")], digits = 3, row.names = FALSE)
    H[[length(H) + 1]] <- h
  }
}

h2t <- do.call(rbind, H); fits <- do.call(rbind, Fm); vcs <- do.call(rbind, V)
write.csv(h2t,  file.path(OUT, "rrm_by_week.csv"),  row.names = FALSE)
write.csv(fits, file.path(OUT, "rrm_model_fit.csv"), row.names = FALSE)
write.csv(vcs,  file.path(OUT, "rrm_varcomp.csv"),   row.names = FALSE)

# Identifiability: pedigree-structured (M2) vs iid-animal (M3) random regression,
# identical fixed effects, residual and parameter count.
# Residual structure per metric by AIC: week-specific for BC and Shannon,
# homogeneous for MAZ.
RES <- c("BC-to-adult" = "het", "Shannon" = "het", "MAZ" = "hom")
mA <- c(het = "M2_Arr_het", hom = "M1_Arr_hom"); mP <- c(het = "M3_Prr_het", hom = "M3h_Prr_hom")
idf <- do.call(rbind, lapply(names(RES), function(t) {
  ff <- fits[fits$trait == t, ]
  lA <- ff$logLik[ff$model == mA[RES[[t]]]]; lP <- ff$logLik[ff$model == mP[RES[[t]]]]
  data.frame(trait = t, residual = RES[[t]], logLik_A = lA, logLik_PE = lP, dLL_A_minus_PE = lA - lP)
}))
write.csv(idf, file.path(OUT, "rrm_A_vs_PE_identifiability.csv"), row.names = FALSE)
# Residual structure and polynomial order (same fixed effects -> REML comparable)
cmp <- fits %>% filter(model %in% c("M1_Arr_hom", "M2_Arr_het", "M2L1_Aleg1_het", "M2L2_Aleg2_het")) %>%
  select(trait, model, logLik, npar, AIC)
write.csv(cmp, file.path(OUT, "rrm_residual_and_order.csv"), row.names = FALSE)
cat("\n=== fits ===\n"); print(fits[, -9], row.names = FALSE)
cat("\n=== A vs PE ===\n"); print(as.data.frame(idf))
