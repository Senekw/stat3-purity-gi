#!/usr/bin/env Rscript
# 09_null.R -- Part B: null-signature benchmark (B.m) and CMS orthogonality (B.o).
#
# THIS SCRIPT FOLLOWS REGISTERED B.m, NOT THE RUN ORDER OF 2026-08-02.
#
# The run order misstated three B.m parameters from memory. The study author
# confirmed on 2026-08-02 that the registered text governs and the run order was
# not to be followed. Recorded here and in NOTES_FOR_REVIEW so the deviation is
# visible rather than absorbed:
#
#   parameter        run order said        B.m registers (USED HERE)
#   ---------------- --------------------- -----------------------------------
#   N                ~1,000                10,000 per cohort
#   set size         140                   152, "or the cohort's realised panel
#                                          size where genes are missing, matched
#                                          per cohort"
#   matching vars    mean expression x     decile of mean log2 expression x
#                    stromal correlation   decile of expression variance,
#                                          "within that cohort"
#
# SET SIZE, both readings run. B.m gives two: the literal 152 (the locked panel)
# and "the cohort's realised panel size". The signature actually tested in 08 is
# the 140-gene list (Amendment 15), so a 140-gene null is the size-matched
# comparison while 152 is the literal registered number. Neither is silently
# preferred: BOTH are drawn at full N and reported side by side.
#
# STATISTIC. B.m: "each null set is scored by the identical B.h pipeline and
# fitted through M2". M1, M2 and M4 are fitted -- M1 and M2 for the crude log-HR
# and M4 because B.m also requires the same construction for attenuation_total
# (= beta2 - beta4), which is undefined without it.
#
# EXPLORATORY, POST-HOC, NOT PART OF B.m: a third draw of 1,000 sets additionally
# matched on the decile of correlation with the stromal score. Added after the
# primary result was seen, at the author's request, because a set matched only on
# expression and variance need not match the panel's stromal loading. Labelled
# EXPLORATORY_POSTHOC in every output row and never mixed with the registered
# result.
#
# COMPUTATIONAL DECISION: no per-signature bootstrap. Point estimates only -- the
# null distribution IS the reference distribution, so a per-signature interval
# would not be consumed by any p-value computed here.
#
# PERFORMANCE. The per-gene z-scoring in B.h does not depend on set membership, so
# the z-scored matrix is computed ONCE per cohort and each null set is a row
# subset of it. Asserted numerically identical to per-set score_cohort() on a
# sample of sets before it is relied on (see the identity gate below).
#
# The null sets are scored and fitted by the same functions as the real panel,
# sourced from the committed scripts and asserted against their definitions there
# rather than reimplemented.
#
# WHICH SCRIPTS: the run order said "reuse the functions from 06 and 07". The
# scoring functions are indeed in 07 (score_cohort, expression_log2tpm,
# zero_variance_genes), but the MODEL functions are in 08, not 06 --
# complete_case_set, fit_cohort, model_formula, epv_band, meta_one. 06 computes
# purity and defines no function this path needs. Sourcing 07 + 08 is therefore
# what "the identical code path as the real panel" requires; 06 + 07 would not
# reach the models at all. Recorded because the run order named different files.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(survival)
  library(metafor)
  library(withr)
})

halt <- function(section, ...)
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)

assert_n <- function(observed, expected, section, what) {
  if (!identical(as.integer(observed), as.integer(expected)))
    halt(section, what, ": expected ", expected, ", observed ", observed)
  message(sprintf("  ok  %-10s %-56s = %d", section, what, observed))
  invisible(TRUE)
}

COHORTS <- c("TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-ESCA",
             "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL")
CRC     <- c("TCGA-COAD", "TCGA-READ")
N_NULL      <- 10000L              # B.m: N = 10,000 per cohort
N_EXPLOR    <- 1000L               # exploratory post-hoc draw only
BASE_SEED   <- 20260731L
NULL_SEED   <- BASE_SEED + 1000L   # B.m: withr::with_seed(20260731 + 1000 + cohort_index)
N_DECILES   <- 10L
# Matching tolerances. NEITHER IS REGISTERED -- B.m states the matching scheme
# but no acceptance threshold, so both are chosen here and labelled. MATCH_TOL is
# the set-level difference of means; PAIRED_TOL the median per-gene paired
# difference. Measured on COAD: a matched draw gives paired 0.32 SD, a completely
# unmatched draw 1.84 SD, so 0.75 separates them decisively.
MATCH_TOL   <- 0.25                # set-level, in SD units of the variable
PAIRED_TOL  <- 0.75                # paired per-gene, in SD units
N_IDENTITY_CHECK <- 25L            # sets used to prove the z-precompute identity

# Observed values from 08, asserted against the committed outputs rather than
# hardcoded as truth -- if 08 is re-run and moves, this halts.
OBS_M1_POOLED <- 0.121174
OBS_ATT_TOTAL <- -0.024766

# ============================================================ PURE FUNCTIONS

#' Decile index (1..k) of x, by quantile.
decile_of <- function(x, k = N_DECILES) {
  br <- stats::quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(x, breaks = unique(br), labels = FALSE, include.lowest = TRUE))
}

#' Per-gene statistics and decile cells for ONE cohort (B.m: "within that cohort").
cohort_gene_stats <- function(X, stromal = NULL) {
  mu <- rowMeans(X)
  v  <- matrixStats_var(X)
  d <- data.frame(gene = rownames(X), mean_expr = mu, var_expr = v,
                  stringsAsFactors = FALSE)
  d$d_expr <- decile_of(d$mean_expr)
  d$d_var  <- decile_of(d$var_expr)
  if (!is.null(stromal)) {
    cs <- as.numeric(suppressWarnings(stats::cor(t(X), stromal)))
    cs[!is.finite(cs)] <- 0
    d$cor_strom <- cs
    d$d_strom   <- decile_of(cs)
  }
  d
}

#' Row variances without a matrixStats dependency.
matrixStats_var <- function(X) {
  m <- rowMeans(X)
  rowSums((X - m)^2) / (ncol(X) - 1L)
}

#' Draw ONE null gene set matched cell-by-cell to `target_genes`.
#'
#' B.m: "each panel gene is matched to a random gene drawn from the same decile of
#' mean log2 expression and the same decile of expression variance within that
#' cohort, sampled without replacement within a set."
#'
#' The panel's own genes are excluded from the candidate pool, so a null set is an
#' ALTERNATIVE signature rather than a partial copy of the panel; including them
#' would pull the null toward the observed statistic and understate the contrast.
#' Cell exhaustion widens to the expression decile, then to the whole pool, and
#' each widening is counted so the achieved match is audited, not assumed.
draw_null_set <- function(tg_cells, cell_index, pool_genes, exclude_idx) {
  n <- nrow(tg_cells)
  out <- integer(n); used <- logical(length(pool_genes))
  used[exclude_idx] <- TRUE
  widen1 <- 0L; widen2 <- 0L
  keys <- apply(tg_cells, 1, paste, collapse = ":")
  for (i in seq_len(n)) {
    cand <- cell_index$by_cell[[keys[i]]]
    cand <- cand[!used[cand]]
    if (!length(cand)) {
      widen1 <- widen1 + 1L
      cand <- cell_index$by_expr[[as.character(tg_cells[i, 1])]]
      cand <- cand[!used[cand]]
    }
    if (!length(cand)) {
      widen2 <- widen2 + 1L
      cand <- which(!used)
    }
    if (!length(cand)) halt("B.m", "the null gene pool is exhausted")
    pick <- if (length(cand) == 1L) cand else cand[sample.int(length(cand), 1L)]
    out[i] <- pick; used[pick] <- TRUE
  }
  list(idx = out, genes = pool_genes[out], widen_cell = widen1, widen_all = widen2)
}

#' Pre-index the pool by decile cell, so a draw is a lookup rather than a scan.
#' Built once per cohort; the sampling it feeds is identical to scanning the pool.
build_cell_index <- function(gs, vv) {
  k <- do.call(paste, c(as.list(gs[vv]), sep = ":"))
  list(by_cell = split(seq_len(nrow(gs)), k),
       by_expr = split(seq_len(nrow(gs)), as.character(gs[[vv[1]]])))
}

#' Score from a PRE-Z-SCORED matrix: mean across the set, then unit SD.
#'
#' Identical by construction to score_cohort(X, genes) because the z-scoring in
#' B.h is per gene within cohort and does not depend on which genes are in the
#' set. The identity is asserted numerically before this is used.
score_from_z <- function(Z, genes) {
  s <- colMeans(Z[genes, , drop = FALSE])
  s / stats::sd(s)
}
#' Empirical p-value, B.m's add-one construction.
#' B.m: p_emp = (1 + #{|beta_null| >= |beta_observed|}) / (1 + N), two-sided,
#' add-one so p is never exactly zero.
#'
#' D7: the denominator is N -- the number of draws REQUESTED -- not the number
#' that happened to pool. Using the latter would drop non-poolable draws from the
#' denominator non-randomly. `side`: "two" is the registered primary; "upper" and
#' "lower" are directional companions, and for attenuation the LOWER tail is the
#' one matching B.m's question ("whether this panel attenuates more than an
#' arbitrary panel would" -- more attenuation is a larger positive value, so a
#' negative observed value sits in the lower tail).
p_emp <- function(null_vals, observed, side = c("two", "upper", "lower"), N = NULL) {
  side <- match.arg(side)
  v <- null_vals[is.finite(null_vals)]
  n_pooled <- length(v)
  denom <- if (is.null(N)) n_pooled else N
  hits <- switch(side,
                 two   = sum(abs(v) >= abs(observed)),
                 upper = sum(v >= observed),
                 lower = sum(v <= observed))
  list(p = (1 + hits) / (1 + denom), hits = hits, n = n_pooled, denom = denom)
}

# ==================================================================== driver
if (sys.nframe() == 0L) {

  OUTDIR <- "output"; dir.create(OUTDIR, showWarnings = FALSE)
  message("\n== B.m/B.o  null-signature benchmark and CMS orthogonality ==")
  message("  ..  seed for the null draw = ", NULL_SEED,
          " (B.m: 20260731 + 1000, a separate stream from the bootstrap)")

  # ---- reuse the real panel's code path, asserted ------------------------
  e07 <- new.env(); sys.source("07_score.R", envir = e07)
  e08 <- new.env(); sys.source("08_survival.R", envir = e08)
  fns <- list(score_cohort = e07$score_cohort,
              expression_log2tpm = e07$expression_log2tpm,
              zero_variance_genes = e07$zero_variance_genes,
              complete_case_set = e08$complete_case_set,
              fit_cohort = e08$fit_cohort,
              model_formula = e08$model_formula,
              epv_band = e08$epv_band,
              meta_one = e08$meta_one,
              META_ELIGIBLE = e08$META_ELIGIBLE)
  for (n in names(fns))
    if (is.null(fns[[n]])) halt("B.m", "function '", n, "' was not found in 07/08")
  # The earlier form of this guard compared fns$score_cohort to e07$score_cohort
  # -- the object it had just been copied from -- so it could never be FALSE. That
  # is the same class of defect the audits have been catching elsewhere, and it
  # was here in the check meant to prove this script reuses the real panel's code.
  #
  # The real question is whether the functions bound here are the ones the COMMITTED
  # scripts define. identical() on closures compares environments too, so two
  # sourcings of the same file are NOT identical; the comparison must be on the
  # function BODY. A body mismatch means 07/08 changed under this script's feet.
  chk_env <- new.env(); sys.source("07_score.R", envir = chk_env)
  chk8    <- new.env(); sys.source("08_survival.R", envir = chk8)
  body_same <- function(a, b) identical(deparse(body(a)), deparse(body(b)))
  for (nm in c("score_cohort", "expression_log2tpm", "zero_variance_genes"))
    if (!body_same(fns[[nm]], chk_env[[nm]]))
      halt("B.m", "'", nm, "' does not match the definition in 07_score.R")
  for (nm in c("complete_case_set", "fit_cohort", "model_formula", "epv_band", "meta_one"))
    if (!body_same(fns[[nm]], chk8[[nm]]))
      halt("B.m", "'", nm, "' does not match the definition in 08_survival.R")
  # And prove the guard can fire, rather than trusting that it would.
  if (body_same(fns$score_cohort, fns$meta_one))
    halt("B.m", "the body comparison cannot discriminate two different functions")
  message("  ok  B.m        all 8 reused functions match their committed definitions ",
          "in 07_score.R / 08_survival.R (body comparison, negative control passed)")

  # ---- inputs -------------------------------------------------------------
  clin <- read.csv(file.path(OUTDIR, "clinical_analysis_set.csv"), stringsAsFactors = FALSE)
  pur  <- read.csv(file.path(OUTDIR, "purity_per_patient.csv"), stringsAsFactors = FALSE)
  scr  <- read.csv(file.path(OUTDIR, "scores_per_patient.csv"), stringsAsFactors = FALSE)
  real_m1  <- read.csv(file.path(OUTDIR, "meta_analysis.csv"), stringsAsFactors = FALSE)
  obs_m1   <- real_m1$est[real_m1$analysis == "M1"]
  obs_att  <- real_m1$est[real_m1$analysis == "attenuation_total"]
  if (abs(obs_m1 - OBS_M1_POOLED) > 1e-6 || abs(obs_att - OBS_ATT_TOTAL) > 1e-6)
    halt("B.m", "08's committed pooled estimates moved (M1 ", obs_m1, ", att ",
         obs_att, "); the null benchmark must be re-targeted deliberately")
  message("  ok  B.m        observed targets read from 08: pooled M1 = ", obs_m1,
          ", attenuation_total = ", obs_att)

  # B1: LIKE-FOR-LIKE COMPARISON TARGET. 08's -0.024766 was pooled with a PAIRED
  # BOOTSTRAP SE; every null here is pooled with sqrt(se2^2 + se4^2), which is
  # ~1.7x larger on this data because it ignores the positive covariance between
  # b2 and b4. Different weights give a different pooled point estimate, not just
  # a different interval, so comparing the bootstrap-pooled observed value against
  # a naive-pooled null distribution is not like-for-like. p_atten is therefore
  # computed against the observed value RE-POOLED WITH THE NULLS' OWN ESTIMATOR.
  # The bootstrap-pooled value remains 08's reported estimate and is carried in
  # the output beside it.
  pcx <- read.csv(file.path(OUTDIR, "survival_per_cohort.csv"), stringsAsFactors = FALSE)
  atx <- read.csv(file.path(OUTDIR, "attenuation_per_cohort.csv"), stringsAsFactors = FALSE)
  ecx <- atx$cohort[atx$meta_eligible]
  gb <- function(m, w) pcx[[w]][pcx$model == m & pcx$cohort %in% ecx]
  m_naive <- fns$meta_one(gb("M2","beta") - gb("M4","beta"),
                          gb("M2","se")^2 + gb("M4","se")^2, ecx)
  if (is.null(m_naive)) halt("B.m", "could not re-pool the observed attenuation")
  obs_att_naive <- m_naive$est
  message(sprintf("  ..  B.m        observed attenuation re-pooled with the nulls' estimator: %.6f (08 reported %.6f with the paired bootstrap SE)",
                  obs_att_naive, obs_att))
  obs_m2 <- real_m1$est[real_m1$analysis == "M2"]
  if (!length(obs_m2)) halt("B.m", "08's pooled M2 estimate was not found")

  d <- merge(clin, pur[, c("barcode", "purity", "purity_estimate", "StromalScore")], by = "barcode")
  d <- merge(d, scr[, c("barcode", "score", "stromal_score")], by = "barcode")
  d$stage_group <- factor(d$stage_group, levels = c("I/II", "III/IV", "missing"))
  d$sex <- factor(d$sex)

  genes140 <- read.csv("data/panel/final_gene_list_140.csv", stringsAsFactors = FALSE)$gene
  assert_n(length(genes140), 140L, "B.m", "genes in the tested signature")

  # ---- expression, and the eligible null pool -----------------------------
  expr <- lapply(COHORTS, function(cc) {
    se <- readRDS(sprintf("data/tcga/%s_se.rds", cc))
    x <- fns$expression_log2tpm(se, clin$barcode[clin$cohort == cc])
    rm(se); gc(verbose = FALSE); x
  })
  names(expr) <- COHORTS

  # R15: the eligible pool is built PER COHORT below (B.m: "within that
  # cohort"); an aggregate cross-cohort pool would be a number never used.

  # ---- frames: the complete-case set, fixed once per cohort ---------------
  # Built from the REAL panel's score column, which is non-NA everywhere, so this
  # is exactly the patient set 08 used. Every null signature is then evaluated on
  # that same set: a like-for-like comparison in which only the gene set varies.
  frames <- list(); use_sex <- list()
  for (cc in COHORTS) {
    dc <- d[d$cohort == cc, , drop = FALSE]
    us <- isTRUE(unique(dc$use_sex))
    fr <- fns$complete_case_set(dc, us)
    fr$stage_group <- droplevels(fr$stage_group)
    frames[[cc]] <- fr; use_sex[[cc]] <- us
  }

  # ---- per-cohort z-scored matrices, gene stats and decile cells -----------
  # B.m matches "within that cohort", so deciles are computed per cohort and the
  # null sets are drawn per cohort. The z-scored matrix is built once here.
  message("  ..  B.m        building per-cohort z-matrices and decile cells")
  ZC <- list(); GS <- list(); POOL <- list()
  for (cc in COHORTS) {
    X <- expr[[cc]][, frames[[cc]]$barcode, drop = FALSE]
    keep <- setdiff(rownames(X), fns$zero_variance_genes(X, rownames(X)))
    X <- X[keep, , drop = FALSE]
    mu <- rowMeans(X); sdv <- sqrt(matrixStats_var(X))
    ZC[[cc]] <- (X - mu) / sdv
    strom <- d$stromal_score[match(colnames(X), d$barcode)]
    GS[[cc]] <- cohort_gene_stats(X, strom)
    POOL[[cc]] <- rownames(X)
  }
  message(sprintf("  ..  B.m        eligible pool per cohort: %s",
                  paste(sprintf("%s %d", sub("^TCGA-", "", COHORTS),
                                vapply(POOL, length, integer(1))), collapse = " | ")))

  # ---- IDENTITY GATE: the z-precompute must equal per-set score_cohort ------
  # R14: exercise EVERY cohort and BOTH set sizes, not one cohort at one size.
  set.seed(BASE_SEED)
  chk <- unlist(lapply(COHORTS, function(cc) vapply(c(152L, 140L), function(k)
    max(vapply(seq_len(N_IDENTITY_CHECK), function(i) {
      g <- sample(POOL[[cc]], k)
      max(abs(fns$score_cohort(expr[[cc]][, frames[[cc]]$barcode, drop = FALSE], g) -
              score_from_z(ZC[[cc]], g)))
    }, numeric(1))), numeric(1))))
  message(sprintf("  ..  B.m        z-precompute identity over %d cohorts x 2 sizes x %d sets: max |diff| = %.3e",
                  length(COHORTS), N_IDENTITY_CHECK, max(chk)))
  if (max(chk) > 1e-10)
    halt("B.m", "the z-precompute is NOT identical to per-set score_cohort (max ",
         format(max(chk)), "); it must not be used")
  message("  ok  B.m        z-precompute proven identical to the registered path")

  # ---- the signature under test, per cohort --------------------------------
  # B.m: "152 genes, matching the panel exactly (or the cohort's realised panel
  # size where genes are missing, matched per cohort)". Both readings are run.
  panel152 <- read.csv("data/panel/panel_locked.csv", stringsAsFactors = FALSE)$gene
  assert_n(length(panel152), 152L, "B.m", "genes in the locked panel")

  # B.m allows "the cohort's realised panel size where genes are missing, matched
  # per cohort". VERIFIED INERT: all 152 panel genes survive the annotation and
  # zero-variance filters in all seven cohorts (COAD/READ/STAD/ESCA/PAAD/LIHC/CHOL
  # each 152/152), so the realised size IS 152 everywhere and the clause never
  # bites. That realised-size clause is about genes MISSING FROM A COHORT, which
  # is a different thing from the 140-gene post-exclusion list; the two are not
  # conflated here. The 140 configuration is included because it is the signature
  # 08 actually tested, and is labelled by size in every output row.
  SIZES <- list(
    registered_152 = panel152,   # B.m's literal number; realised size in all 7
    tested_140     = genes140)   # the signature 08 tested (Amendment 15)
  EXCLUDE_ALWAYS <- union(panel152, genes140)

  SCHEMES <- list(
    registered      = c("d_expr", "d_var"),     # B.m: expression x variance
    # D3: ADDITIONALLY matched -- expression x variance x stromal correlation.
    # An earlier revision dropped d_var here, which would have reproduced the
    # withdrawn run order's matching under a label claiming the opposite.
    EXPLORATORY_POSTHOC = c("d_expr", "d_var", "d_strom"))

  run_cohort_draws <- function(target_genes, vv, n_draw, seed_off) {
    # returns, per cohort, a matrix of null scores (n_draw x n_patients) plus the
    # achieved-match record
    lapply(seq_along(COHORTS), function(ci) {
      cc <- COHORTS[ci]
      gs <- GS[[cc]]
      tg <- gs[match(intersect(target_genes, gs$gene), gs$gene), , drop = FALSE]
      cell_index <- build_cell_index(gs, vv)
      # D4: exclude BOTH gene lists from EVERY configuration's candidate pool.
      # Excluding only the configuration's own list would let the 12 genes that
      # Amendment 15 and rule 2 removed be drawn into null sets in the 140
      # configuration, pulling the null toward the observed statistic.
      excl_idx <- which(gs$gene %in% EXCLUDE_ALWAYS)
      tgm <- as.matrix(tg[, vv])
      sets <- with_seed(NULL_SEED + ci + seed_off, lapply(seq_len(n_draw), function(i)
        draw_null_set(tgm, cell_index, gs$gene, excl_idx)))
      ach <- do.call(rbind, lapply(seq_along(sets), function(i) {
        m <- gs[match(sets[[i]]$genes, gs$gene), ]
        data.frame(i = i,
                   # set-level difference of means
                   d_expr  = mean(m$mean_expr) - mean(tg$mean_expr),
                   d_var   = mean(m$var_expr)  - mean(tg$var_expr),
                   d_strom = if ("cor_strom" %in% names(gs))
                     mean(m$cor_strom) - mean(tg$cor_strom) else NA_real_,
                   # R9: PAIRED per-gene discrepancy -- each null gene against the
                   # panel gene it replaced. Averaging over 152 genes shrinks the
                   # set-level statistic by ~1/sqrt(152), so the paired version is
                   # the more discriminating test of the match (measured: ~6x
                   # tighter for matched than unmatched draws).
                   paired_expr = median(abs(m$mean_expr - tg$mean_expr)),
                   paired_var  = median(abs(m$var_expr  - tg$var_expr)),
                   paired_strom = if ("cor_strom" %in% names(gs))
                     median(abs(m$cor_strom - tg$cor_strom)) else NA_real_,
                   widen_cell = sets[[i]]$widen_cell,
                   widen_all  = sets[[i]]$widen_all, stringsAsFactors = FALSE)
      }))
      scores <- vapply(sets, function(z) score_from_z(ZC[[cc]], z$genes),
                       numeric(ncol(ZC[[cc]])))
      list(cohort = cc, n_realised = nrow(tg), achieved = ach, scores = scores,
           sd_expr = stats::sd(gs$mean_expr), sd_var = stats::sd(gs$var_expr),
           sd_strom = if ("cor_strom" %in% names(gs)) stats::sd(gs$cor_strom) else NA_real_)
    })
  }

  #' Fit M1, M2, M4 for every null draw in one cohort, from its score matrix.
  fit_null_cohort <- function(cc, scores) {
    fr <- frames[[cc]]; us <- use_sex[[cc]]
    ev <- sum(fr$event == 1L)
    bands <- vapply(c("M1","M2","M4"),
                    function(m) fns$epv_band(ev, fns$model_formula(m, us)$params),
                    character(1))
    out <- lapply(seq_len(ncol(scores)), function(j) {
      fr$score <- scores[, j]
      f <- fns$fit_cohort(fr, us)
      vapply(c("M1","M2","M4"), function(m) {
        if (!isTRUE(f[[m]]$ok) || bands[[m]] == "do_not_fit") c(NA_real_, NA_real_)
        else c(f[[m]]$beta, f[[m]]$se)
      }, numeric(2))
    })
    list(beta = t(vapply(out, function(z) z[1, ], numeric(3))),
         se   = t(vapply(out, function(z) z[2, ], numeric(3))),
         bands = bands, events = ev)
  }

  # ---- run every configuration --------------------------------------------
  CONFIGS <- list()
  for (sz in names(SIZES)) for (sch in names(SCHEMES)) {
    n_draw <- if (sch == "registered") N_NULL else N_EXPLOR
    if (sch != "registered" && sz != "tested_140") next   # exploratory: one size only
    # B2: an explicit ROLE, declared before the run. registered_152 is THE
    # registered primary (B.m's literal 152, and the realised size in all seven
    # cohorts). tested_140 is a size sensitivity. The stromal-matched draw is
    # post-hoc. Two full-N p-values must not arrive with equal standing.
    role <- if (sz == "registered_152" && sch == "registered") "primary_registered"
            else if (sch == "registered") "size_sensitivity"
            else "exploratory_posthoc"
    off <- switch(role, primary_registered = 0L, size_sensitivity = 200L, 500L)
    CONFIGS[[paste(sz, sch, sep = "|")]] <- list(size = sz, scheme = sch, role = role,
                                                 n = n_draw, seed_off = off)
  }
  message("  ..  B.m        configurations: ",
          paste(names(CONFIGS), collapse = "  ;  "))

  t_start <- Sys.time(); first_done <- FALSE
  RES <- list()
  for (key in names(CONFIGS)) {
    cfg <- CONFIGS[[key]]
    message("\n  ..  B.m        === ", key, "  (N = ", cfg$n, ") ===")
    # Seed offsets are distinct PER CONFIGURATION. The registered configuration
    # (152 genes, expression x variance) keeps offset 0, so its stream is exactly
    # B.m's withr::with_seed(20260731 + 1000 + cohort_index). Without distinct
    # offsets, registered_152 and tested_140 would share a stream per cohort and
    # their null sets would be coupled rather than independent draws.
    draws <- run_cohort_draws(SIZES[[cfg$size]], SCHEMES[[cfg$scheme]], cfg$n,
                              seed_off = cfg$seed_off)
    names(draws) <- COHORTS
    fits <- list()
    for (cc in COHORTS) {
      tc <- Sys.time()
      fits[[cc]] <- fit_null_cohort(cc, draws[[cc]]$scores)
      el <- as.numeric(difftime(Sys.time(), tc, units = "secs"))
      message(sprintf("      %-6s %d draws fitted in %.1f s", sub("^TCGA-","",cc),
                      cfg$n, el))
      if (!first_done) {
        first_done <- TRUE
        message(sprintf("      ESTIMATED TOTAL: ~%.1f min (%.1f s x %d cohorts x %d configs, N-scaled)",
                        el * length(COHORTS) * length(CONFIGS) / 60,
                        el, length(COHORTS), length(CONFIGS)))
      }
    }
    RES[[key]] <- list(cfg = cfg, draws = draws, fits = fits)
    # R17: checkpoint each configuration as it lands, so an error in a later one
    # does not discard hours of completed fits.
    saveRDS(RES[[key]], file.path(OUTDIR, paste0("null_ckpt_",
            gsub("[^A-Za-z0-9]", "_", key), ".rds")))
  }

  # ---- pooled statistics per null draw ------------------------------------
  pool_nulls <- function(R) {
    elig <- vapply(COHORTS, function(cc) isTRUE(fns$META_ELIGIBLE[[cc]]), logical(1))
    ec <- COHORTS[elig]
    n <- R$cfg$n
    m1 <- rep(NA_real_, n); m1lo <- rep(NA_real_, n); att <- rep(NA_real_, n)
    m2 <- rep(NA_real_, n)   # D8: B.m fits through M2; the pooled M2 null is reported
    for (j in seq_len(n)) {
      b1 <- vapply(ec, function(cc) R$fits[[cc]]$beta[j, 1], numeric(1))
      s1 <- vapply(ec, function(cc) R$fits[[cc]]$se[j, 1],   numeric(1))
      # D7: a null must pool over the FULL eligible set, or it is not comparable
      # with an observed value pooled over six cohorts. Partial pools are counted
      # and excluded from the numerator, never silently pooled over three.
      k1 <- is.finite(b1) & is.finite(s1) & s1 > 0
      if (all(k1)) {
        mm <- fns$meta_one(b1[k1], s1[k1]^2, ec[k1])
        if (!is.null(mm)) { m1[j] <- mm$est; m1lo[j] <- mm$ci_lo_hksj }
      }
      b2 <- vapply(ec, function(cc) R$fits[[cc]]$beta[j, 2], numeric(1))
      s2m <- vapply(ec, function(cc) R$fits[[cc]]$se[j, 2], numeric(1))
      k2 <- is.finite(b2) & is.finite(s2m) & s2m > 0
      if (all(k2)) {
        mm2 <- fns$meta_one(b2[k2], s2m[k2]^2, ec[k2])
        if (!is.null(mm2)) m2[j] <- mm2$est
      }
      b4 <- vapply(ec, function(cc) R$fits[[cc]]$beta[j, 3], numeric(1))
      s2 <- vapply(ec, function(cc) R$fits[[cc]]$se[j, 2], numeric(1))
      s4 <- vapply(ec, function(cc) R$fits[[cc]]$se[j, 3], numeric(1))
      ka <- is.finite(b2) & is.finite(b4) & is.finite(s2) & is.finite(s4)
      if (all(ka)) {
        # No per-signature bootstrap (declared decision), so the paired SE is
        # unavailable and the components are combined as sqrt(se2^2 + se4^2).
        # This OVERSTATES each null's variance because b2 and b4 are positively
        # correlated. It is applied identically to every null set and affects the
        # WEIGHTS, not the per-cohort differences being pooled.
        ma <- fns$meta_one((b2 - b4)[ka], (s2^2 + s4^2)[ka], ec[ka])
        if (!is.null(ma)) att[j] <- ma$est
      }
    }
    list(m1 = m1, m1_lo = m1lo, att = att, m2 = m2, n_draws = n)
  }
  POOLED <- lapply(RES, pool_nulls)

  # ---- matching quality and the tolerance gate ----------------------------
  tol_tab <- do.call(rbind, lapply(names(RES), function(key) {
    R <- RES[[key]]
    do.call(rbind, lapply(COHORTS, function(cc) {
      dr <- R$draws[[cc]]; a <- dr$achieved
      data.frame(config = key, cohort = cc,
                 n_realised_set = dr$n_realised,
                 matched_on = paste(SCHEMES[[R$cfg$scheme]], collapse = " x "),
                 med_abs_d_expr_sd = round(median(abs(a$d_expr)) / dr$sd_expr, 4),
                 med_abs_d_var_sd  = round(median(abs(a$d_var))  / dr$sd_var, 4),
                 med_abs_d_strom_sd = round(median(abs(a$d_strom)) / dr$sd_strom, 4),
                 p95_abs_d_expr_sd = round(stats::quantile(abs(a$d_expr), .95) / dr$sd_expr, 4),
                 paired_expr_sd = round(median(a$paired_expr) / dr$sd_expr, 4),
                 paired_var_sd  = round(median(a$paired_var)  / dr$sd_var, 4),
                 paired_strom_sd = round(median(a$paired_strom) / dr$sd_strom, 4),
                 widen_cell_median = median(a$widen_cell),
                 widen_all_median  = median(a$widen_all),
                 stringsAsFactors = FALSE)
    }))
  }))
  write.csv(tol_tab, file.path(OUTDIR, "null_matching_quality.csv"), row.names = FALSE)
  for (i in seq_len(nrow(tol_tab))) {
    vv <- SCHEMES[[RES[[tol_tab$config[i]]]$cfg$scheme]]
    chk <- c(d_expr = tol_tab$med_abs_d_expr_sd[i],
             d_var  = tol_tab$med_abs_d_var_sd[i],
             d_strom = tol_tab$med_abs_d_strom_sd[i])[vv]
    # gate on the PAIRED statistic as well; PAIRED_TOL is looser because a paired
    # per-gene difference is not shrunk by averaging. Both are unregistered
    # thresholds chosen here and labelled as such.
    chkp <- c(d_expr = tol_tab$paired_expr_sd[i],
              d_var  = tol_tab$paired_var_sd[i],
              d_strom = tol_tab$paired_strom_sd[i])[vv]
    if (any(chkp > PAIRED_TOL, na.rm = TRUE))
      halt("B.m", "STOP CONDITION: ", tol_tab$config[i], " / ", tol_tab$cohort[i],
           " missed the PAIRED matching tolerance of ", PAIRED_TOL, " SD: ",
           paste(names(chkp), round(chkp, 3), collapse = ", "))
    if (any(chk > MATCH_TOL, na.rm = TRUE))
      halt("B.m", "STOP CONDITION: ", tol_tab$config[i], " / ", tol_tab$cohort[i],
           " missed the ", MATCH_TOL, " SD tolerance on a MATCHED variable: ",
           paste(names(chk), round(chk, 3), collapse = ", "))
  }
  message("  ok  B.m        every configuration met the ", MATCH_TOL,
          " SD matching tolerance on its matched variables")

  # ---- distributions, p-values, per-cohort percentiles --------------------
  q5 <- function(v) stats::quantile(v[is.finite(v)], c(0, .25, .5, .75, 1))
  pctile <- function(v, o) round(100 * mean(v[is.finite(v)] <= o), 3)

  dist_tab <- do.call(rbind, lapply(names(RES), function(key) {
    P <- POOLED[[key]]
    rbind(
      data.frame(config = key, role = RES[[key]]$cfg$role,
                 statistic = "pooled_M1_logHR", n_null = sum(is.finite(P$m1)),
                 min = q5(P$m1)[1], q25 = q5(P$m1)[2], median = q5(P$m1)[3],
                 mean = mean(P$m1, na.rm = TRUE), q75 = q5(P$m1)[4], max = q5(P$m1)[5],
                 observed_reported_08 = obs_m1, observed_percentile = pctile(P$m1, obs_m1),
                 stringsAsFactors = FALSE),
      data.frame(config = key, role = RES[[key]]$cfg$role,
                 statistic = "pooled_M2_logHR", n_null = sum(is.finite(P$m2)),
                 min = q5(P$m2)[1], q25 = q5(P$m2)[2], median = q5(P$m2)[3],
                 mean = mean(P$m2, na.rm = TRUE), q75 = q5(P$m2)[4], max = q5(P$m2)[5],
                 observed_reported_08 = obs_m2, observed_percentile = pctile(P$m2, obs_m2),
                 stringsAsFactors = FALSE),
      data.frame(config = key, role = RES[[key]]$cfg$role,
                 statistic = "pooled_attenuation_total", n_null = sum(is.finite(P$att)),
                 min = q5(P$att)[1], q25 = q5(P$att)[2], median = q5(P$att)[3],
                 mean = mean(P$att, na.rm = TRUE), q75 = q5(P$att)[4], max = q5(P$att)[5],
                 # the DISTRIBUTION is summarised against 08's reported value;
                 # the p-value uses the re-pooled target. Named apart so the two
                 # cannot be read as the same quantity.
                 observed_reported_08 = obs_att,
                 observed_percentile = pctile(P$att, obs_att),
                 stringsAsFactors = FALSE))
  }))

  p_tab <- do.call(rbind, lapply(names(RES), function(key) {
    P <- POOLED[[key]]
    N <- RES[[key]]$cfg$n
    p1t <- p_emp(P$m1,  obs_m1,        "two",   N); p1u <- p_emp(P$m1,  obs_m1,        "upper", N)
    p2t <- p_emp(P$m2,  obs_m2,        "two",   N); p2u <- p_emp(P$m2,  obs_m2,        "upper", N)
    pat <- p_emp(P$att, obs_att_naive, "two",   N); pal <- p_emp(P$att, obs_att_naive, "lower", N)
    data.frame(config = key, role = RES[[key]]$cfg$role,
               statistic = c("p_crude_M1", "p_adjusted_M2", "p_atten"),
               # DISAMBIGUATED: this is the value the null distribution is
               # compared AGAINST, which for p_atten is the observed attenuation
               # RE-POOLED with the nulls' own estimator (B1). It is deliberately
               # NOT 08's reported -0.024766, which is carried beside it. The two
               # files previously both used a bare `observed` for the two
               # different quantities.
               observed_comparison_target = c(obs_m1, obs_m2, obs_att_naive),
               observed_reported_08 = c(obs_m1, obs_m2, obs_att),
               comparison_estimator = c("same", "same",
                                        "re-pooled with sqrt(se2^2+se4^2), per B1"),
               N_requested = N,
               n_pooled = c(p1t$n, p2t$n, pat$n),
               p_two_sided_PRIMARY = c(p1t$p, p2t$p, pat$p),
               hits_two = c(p1t$hits, p2t$hits, pat$hits),
               p_one_sided_directional = c(p1u$p, p2u$p, pal$p),
               directional_tail = c("upper", "upper", "lower"),
               hits_one = c(p1u$hits, p2u$hits, pal$hits),
               p_granularity = 1 / (1 + N),
               stringsAsFactors = FALSE)
  }))

  ci_tab <- do.call(rbind, lapply(names(RES), function(key) {
    lo <- POOLED[[key]]$m1_lo
    data.frame(config = key, role = RES[[key]]$cfg$role,
               n = sum(is.finite(lo)),
               prop_null_CI_excludes_1 = round(mean(lo[is.finite(lo)] > 0), 4),
               stringsAsFactors = FALSE)
  }))

  real_pc <- read.csv(file.path(OUTDIR, "survival_per_cohort.csv"), stringsAsFactors = FALSE)
  pc_tab <- do.call(rbind, lapply(names(RES), function(key) {
    R <- RES[[key]]
    do.call(rbind, lapply(COHORTS, function(cc) {
      do.call(rbind, lapply(seq_along(c("M1","M2","M4")), function(mi) {
        m <- c("M1","M2","M4")[mi]
        v <- R$fits[[cc]]$beta[, mi]
        o <- real_pc$beta[real_pc$cohort == cc & real_pc$model == m]
        if (!length(o) || is.na(o) || !any(is.finite(v))) return(NULL)
        data.frame(config = key, role = R$cfg$role,
                   cohort = cc, model = m, observed_beta = o,
                   null_median = round(median(v, na.rm = TRUE), 6),
                   null_q025 = round(stats::quantile(v[is.finite(v)], .025), 6),
                   null_q975 = round(stats::quantile(v[is.finite(v)], .975), 6),
                   percentile = pctile(v, o), n_null = sum(is.finite(v)),
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

  write.csv(dist_tab, file.path(OUTDIR, "null_distributions.csv"), row.names = FALSE)
  write.csv(p_tab,    file.path(OUTDIR, "null_pvalues.csv"), row.names = FALSE)
  write.csv(ci_tab,   file.path(OUTDIR, "null_ci_exclusion.csv"), row.names = FALSE)
  write.csv(pc_tab,   file.path(OUTDIR, "null_per_cohort_percentile.csv"), row.names = FALSE)
  saveRDS(list(pooled = POOLED, seed = NULL_SEED, n_registered = N_NULL,
               n_exploratory = N_EXPLOR, configs = names(CONFIGS)),
          file.path(OUTDIR, "null_replicates.rds"))

  message("\n-- null distributions --");  print(dist_tab, row.names = FALSE)
  message("\n-- empirical p-values --");  print(p_tab, row.names = FALSE)
  message("\n-- null pooled-M1 CI excluding 1 --"); print(ci_tab, row.names = FALSE)
  message(sprintf("\n  ..  B.m        total elapsed: %.1f min",
                  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  # ========================================================== B.o  CMS
  message("\n== B.o  CMS orthogonality (COAD, READ) ==")
  .libPaths(c("renv/library-local", .libPaths()))
  if (!requireNamespace("CMScaller", quietly = TRUE))
    halt("B.o", "CMScaller is not available; the registered classifier cannot be run")
  cms_ver <- as.character(utils::packageVersion("CMScaller"))
  message("  ok  B.o        classifier: CMScaller ", cms_ver,
          " (Eide et al. 2017), templates.CMS, single-sample NTP")

  utils::data("templates.CMS", package = "CMScaller", envir = environment())
  tmpl <- get("templates.CMS", envir = environment())

  cms_calls <- do.call(rbind, lapply(CRC, function(cc) {
    x <- expr[[cc]]                        # log2(TPM+1), symbols x samples
    # CMScaller expects Entrez rownames. The shipped templates carry both, so the
    # map comes from the classifier's own object rather than an external
    # annotation that might disagree with it.
    map <- unique(tmpl[, c("probe", "symbol")])
    map <- map[map$symbol %in% rownames(x), , drop = FALSE]
    if (!nrow(map)) halt("B.o", cc, ": no template gene matched the expression matrix")
    # R11: an ambiguous symbol->Entrez map would let one expression vector enter
    # twice under two IDs, or silently take an arbitrary duplicate row.
    if (anyDuplicated(rownames(x)))
      halt("B.o", cc, ": duplicate gene symbols in the expression matrix")
    if (anyDuplicated(map$symbol))
      halt("B.o", cc, ": a symbol maps to more than one template probe")
    cov_by_class <- table(tmpl$class[tmpl$symbol %in% rownames(x)])
    tot_by_class <- table(tmpl$class)
    frac <- cov_by_class[names(tot_by_class)] / as.numeric(tot_by_class)
    if (any(frac < 0.70, na.rm = TRUE))
      halt("B.o", cc, ": template coverage below 70% for ",
           paste(names(frac)[frac < 0.70], collapse = ", "))
    em <- x[map$symbol, , drop = FALSE]
    rownames(em) <- as.character(map$probe)
    # RNAseq = TRUE applies the package's own count-scale adjustment; the input
    # here is already log2(TPM+1), so it is passed as continuous data with
    # RNAseq = FALSE and the package's internal row-centring left to do its work.
    r <- CMScaller::CMScaller(emat = em, rowNames = "entrez", RNAseq = FALSE,
                              nPerm = 1000, seed = NULL_SEED, doPlot = FALSE,
                              verbose = FALSE)
    data.frame(cohort = cc, barcode = rownames(r),
               cms = as.character(r$prediction),
               cms_p = r$p.value, cms_fdr = r$FDR,
               n_template_genes = nrow(em), stringsAsFactors = FALSE)
  }))
  cms_calls$cms[is.na(cms_calls$cms)] <- "unclassified"   # B.o: retained, not dropped

  cms_dist <- as.data.frame(table(cms_calls$cohort, cms_calls$cms))
  names(cms_dist) <- c("cohort", "cms", "n")
  cms_dist$pct <- round(100 * cms_dist$n /
                        ave(cms_dist$n, cms_dist$cohort, FUN = sum), 2)
  print(cms_dist, row.names = FALSE)

  # STOP CONDITION: > 10% unclassified in either cohort
  for (cc in CRC) {
    # R12: the stop condition must test the samples that ENTER THE MODELS, not
    # every barcode in the cohort. Both are reported.
    fr_cl <- cms_calls$cms[match(frames[[cc]]$barcode, cms_calls$barcode)]
    u <- 100 * mean(fr_cl == "unclassified")
    u_all <- cms_dist$pct[cms_dist$cohort == cc & cms_dist$cms == "unclassified"]
    message(sprintf("  ..  B.o        %s unclassified: %.2f%% of the analysis frame, %.2f%% of the full cohort",
                    cc, u, if (length(u_all)) u_all else 0))
    message(sprintf("  ..  B.o        %s unclassified: %.2f%%", cc, u))
    if (u > 10)
      halt("B.o", "STOP CONDITION: CMS classification failed for ", round(u, 2),
           "% of ", cc, ", above the 10% ceiling")
  }
  message("  ok  B.o        unclassified below the 10% ceiling in both cohorts")

  # ---- refit M2 and M4 with CMS added -------------------------------------
  cms_fit <- do.call(rbind, lapply(CRC, function(cc) {
    fr <- frames[[cc]]
    fr$cms <- factor(cms_calls$cms[match(fr$barcode, cms_calls$barcode)])
    if (anyNA(fr$cms)) halt("B.o", cc, ": a patient has no CMS call")
    fr$cms <- droplevels(fr$cms)
    us <- use_sex[[cc]]
    do.call(rbind, lapply(c("M2", "M4"), function(m) {
      base <- fns$model_formula(m, us)
      rhs  <- paste(deparse(base$formula[[3]]), collapse = " ")
      f0 <- survival::coxph(base$formula, data = fr, ties = "efron")
      # D5: Amendment 8's EPV rule applies PER MODEL. `+ cms` adds
      # nlevels(cms) - 1 parameters, so the augmented model has its own band.
      p_aug <- base$params + nlevels(fr$cms) - 1L
      band_aug <- fns$epv_band(sum(fr$event == 1L), p_aug)
      f1 <- if (band_aug == "do_not_fit") NULL else
        survival::coxph(stats::as.formula(paste("Surv(time, event) ~", rhs, "+ cms")),
                        data = fr, ties = "efron")
      b0 <- stats::coef(f0)["score"]; s0 <- sqrt(diag(stats::vcov(f0)))["score"]
      b1 <- if (is.null(f1)) NA_real_ else stats::coef(f1)["score"]
      s1 <- if (is.null(f1)) NA_real_ else sqrt(diag(stats::vcov(f1)))["score"]
      data.frame(cohort = cc, model = m, n = nrow(fr), events = sum(fr$event == 1L),
                 beta_without_cms = round(unname(b0), 6), se_without = round(unname(s0), 6),
                 beta_with_cms    = round(unname(b1), 6), se_with = round(unname(s1), 6),
                 params_with_cms = p_aug, epv_with_cms = round(sum(fr$event == 1L) / p_aug, 2),
                 band_with_cms = band_aug, fitted_with_cms = !is.null(f1),
                 HR_without = round(exp(unname(b0)), 4), HR_with = round(exp(unname(b1)), 4),
                 delta_beta = round(unname(b1 - b0), 6),
                 cms_levels = paste(levels(fr$cms), collapse = "/"),
                 stringsAsFactors = FALSE)
    }))
  }))
  # attenuation with and without CMS in the model
  cms_att <- do.call(rbind, lapply(CRC, function(cc) {
    z <- cms_fit[cms_fit$cohort == cc, ]
    data.frame(cohort = cc,
               attenuation_without_cms = round(z$beta_without_cms[z$model == "M2"] -
                                               z$beta_without_cms[z$model == "M4"], 6),
               attenuation_with_cms    = round(z$beta_with_cms[z$model == "M2"] -
                                               z$beta_with_cms[z$model == "M4"], 6),
               stringsAsFactors = FALSE)
  }))

  # ---- score vs CMS4 membership -------------------------------------------
  cms4 <- do.call(rbind, lapply(CRC, function(cc) {
    fr <- frames[[cc]]
    cl <- cms_calls$cms[match(fr$barcode, cms_calls$barcode)]
    is4 <- as.integer(cl == "CMS4")
    data.frame(cohort = cc, n = nrow(fr), n_CMS4 = sum(is4),
               r_point_biserial = round(stats::cor(fr$score, is4), 4),
               mean_score_CMS4 = round(mean(fr$score[is4 == 1]), 4),
               mean_score_other = round(mean(fr$score[is4 == 0]), 4),
               r2_score_on_cms = round(summary(stats::lm(score ~ factor(cl),
                                                         data = fr))$r.squared, 4),
               stringsAsFactors = FALSE)
  }))

  # ---- B.o.1: score tertiles x CMS, with chi-square -----------------------
  cms_xtab <- do.call(rbind, lapply(CRC, function(cc) {
    fr <- frames[[cc]]
    cl <- cms_calls$cms[match(fr$barcode, cms_calls$barcode)]
    ter <- cut(fr$score, breaks = stats::quantile(fr$score, c(0, 1/3, 2/3, 1)),
               labels = c("T1", "T2", "T3"), include.lowest = TRUE)
    # R13: B.o registers "cross-tabulation of score tertiles against CMS1-4".
    # `unclassified` is retained in the MODELS (registered) but is not a CMS
    # class, so the registered test is on CMS1-4; the full table is reported too.
    tb_all <- table(ter, cl)
    tb <- tb_all[, colnames(tb_all) != "unclassified", drop = FALSE]
    ct <- suppressWarnings(stats::chisq.test(tb))
    # simulate.p.value guards the small expected counts that make the asymptotic
    # chi-square unreliable in READ; both p-values are reported.
    cs <- suppressWarnings(stats::chisq.test(tb, simulate.p.value = TRUE, B = 10000))
    z <- as.data.frame(tb); names(z) <- c("tertile", "cms", "n")
    cbind(cohort = cc, z,
          chisq = round(unname(ct$statistic), 4), df = unname(ct$parameter),
          p_asymptotic = signif(ct$p.value, 4),
          p_simulated = signif(cs$p.value, 4),
          min_expected = round(min(ct$expected), 3))
  }))

  # ---- B.o.3b: stratified by CMS4 vs non-CMS4 -----------------------------
  cms_strat <- do.call(rbind, lapply(CRC, function(cc) {
    fr <- frames[[cc]]
    cl <- cms_calls$cms[match(fr$barcode, cms_calls$barcode)]
    fr$is4 <- ifelse(cl == "CMS4", "CMS4", "non_CMS4")
    us <- use_sex[[cc]]
    do.call(rbind, lapply(c("CMS4", "non_CMS4"), function(g) {
      sub <- fr[fr$is4 == g, , drop = FALSE]
      ev <- sum(sub$event == 1L)
      do.call(rbind, lapply(c("M1", "M2", "M4"), function(m) {
        mf <- fns$model_formula(m, us)
        band <- fns$epv_band(ev, mf$params)
        if (nrow(sub) < 10L || ev < 5L || band == "do_not_fit")
          return(data.frame(cohort = cc, stratum = g, model = m, n = nrow(sub),
                            events = ev, epv = round(ev / mf$params, 2), band = band,
                            fitted = FALSE, beta = NA_real_, se = NA_real_,
                            HR = NA_real_, stringsAsFactors = FALSE))
        sub2 <- sub; sub2$stage_group <- droplevels(sub2$stage_group)
        f <- try(survival::coxph(mf$formula, data = sub2, ties = "efron"), silent = TRUE)
        if (inherits(f, "try-error"))
          return(data.frame(cohort = cc, stratum = g, model = m, n = nrow(sub),
                            events = ev, epv = round(ev / mf$params, 2), band = band,
                            fitted = FALSE, beta = NA_real_, se = NA_real_,
                            HR = NA_real_, stringsAsFactors = FALSE))
        b <- stats::coef(f)["score"]; se <- sqrt(diag(stats::vcov(f)))["score"]
        data.frame(cohort = cc, stratum = g, model = m, n = nrow(sub), events = ev,
                   epv = round(ev / mf$params, 2), band = band, fitted = TRUE,
                   beta = round(unname(b), 6), se = round(unname(se), 6),
                   HR = round(exp(unname(b)), 4), stringsAsFactors = FALSE)
      }))
    }))
  }))

  write.csv(cms_xtab,  file.path(OUTDIR, "cms_tertile_crosstab.csv"), row.names = FALSE)
  write.csv(cms_strat, file.path(OUTDIR, "cms_stratified_models.csv"), row.names = FALSE)
  message("\n-- score tertiles x CMS --");  print(cms_xtab, row.names = FALSE)
  message("\n-- stratified by CMS4 --");    print(cms_strat, row.names = FALSE)

  write.csv(cms_calls, file.path(OUTDIR, "cms_calls.csv"), row.names = FALSE)
  write.csv(cms_dist,  file.path(OUTDIR, "cms_distribution.csv"), row.names = FALSE)
  write.csv(cms_fit,   file.path(OUTDIR, "cms_adjusted_models.csv"), row.names = FALSE)
  write.csv(cms_att,   file.path(OUTDIR, "cms_attenuation.csv"), row.names = FALSE)
  write.csv(cms4,      file.path(OUTDIR, "cms4_association.csv"), row.names = FALSE)
  writeLines(c(paste0("classifier: CMScaller ", cms_ver),
               "templates: templates.CMS (Eide et al. 2017), single-sample NTP",
               paste0("nPerm: 1000, seed: ", NULL_SEED),
               "input: log2(TPM+1), Entrez rownames mapped from the templates' own symbol column",
               "unclassified samples retained as a level, never dropped (B.o)"),
             file.path(OUTDIR, "cms_provenance.txt"))

  message("\n-- CMS-adjusted models --"); print(cms_fit, row.names = FALSE)
  message("\n-- attenuation with/without CMS --"); print(cms_att, row.names = FALSE)
  message("\n-- score vs CMS4 --"); print(cms4, row.names = FALSE)

  message("\nB.m/B.o complete.")
}
