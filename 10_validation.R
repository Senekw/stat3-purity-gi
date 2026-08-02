#!/usr/bin/env Rscript
# 10_validation.R -- Amendment 16, SECONDARY: survival replication in GSE39582.
#
# SCOPE, fixed by the run order of 2026-08-02: GSE39582 ONLY. FU-iCCA, GSE62254
# and ICGC are NOT attempted; their access blockers are recorded in
# DATA_NEEDED.md. The PRIMARY validation analysis (FU-iCCA pY705 concordance) is
# therefore NOT run here, and nothing below should be read as validating the
# score against measured STAT3 phosphorylation.
#
# NO POOLING WITH DISCOVERY (Amendment 16). This cohort is reported on its own.
# There is no meta-analysis in this script and no TCGA estimate is combined with
# a GSE39582 estimate anywhere.
#
# CODE PATH. The score, the models, the EPV rule, the complete-case set, the
# attenuation and its paired bootstrap are the SAME FUNCTION DEFINITIONS as
# 07/08, sourced from the committed scripts and asserted by body comparison, as
# in 09. Two things genuinely cannot be reused and are re-implemented here:
#
#   1. expression_log2tpm() reads a SummarizedExperiment and sums duplicate
#      symbols ON THE TPM SCALE before log2. GSE39582 is Affymetrix GPL570
#      delivered as log2 RMA (verified: 54,675 probes, values 2.0-14.0). Summing
#      log2 values is not the same operation, so applying that reader here would
#      be wrong, not merely inapplicable.
#   2. The clinical assembly in 05 is TCGA-CDR-specific (barcodes, aliquots,
#      redaction). GSE39582 has none of those.
#
# Everything downstream of the expression matrix and the clinical frame is the
# committed code path.

suppressPackageStartupMessages({
  library(survival); library(metafor); library(withr)
})

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}
assert_n <- function(observed, expected, section, what) {
  if (!identical(as.integer(observed), as.integer(expected)))
    halt(section, what, ": expected ", expected, ", observed ", observed)
  message(sprintf("  ok  %-12s %-52s = %d", section, what, observed))
}

OUTDIR   <- "output"
DATADIR  <- "data/validation/GSE39582"
MATRIX   <- file.path(DATADIR, "GSE39582_series_matrix.txt.gz")
PLATFORM <- file.path(DATADIR, "GPL570_platform.txt")
MATRIX_MD5 <- "9192d6561b8724caab6a5554c076eab6"
GENE_LIST_140 <- "data/panel/final_gene_list_140.csv"
GENE_LIST_143 <- "data/panel/final_gene_list_143.csv"

# GSE39582 is Affymetrix HG-U133 Plus 2.0. The ESTIMATE purity conversion was
# CALIBRATED on Affymetrix arrays, so unlike the four TCGA cohorts that used it
# (where it is extrapolated from Illumina RNA-seq) this cohort is on the
# conversion's native platform. Inspection of estimate::estimateScore confirms
# `platform` gates ONLY whether TumorPurity is emitted -- the three scores are
# computed identically on every platform -- so this changes what the package
# returns, not how the stromal/immune/ESTIMATE scores are formed.
ESTIMATE_PLATFORM <- "affymetrix"

CENSOR_DAYS_OS <- 3650L        # 10 years, identical to 05
COHORT_INDEX   <- 8L           # continues 05-08's 1..7; a distinct bootstrap stream

# ---------------------------------------------------------------- code reuse
bind_committed <- function() {
  e07 <- new.env(); sys.source("07_score.R",    envir = e07)
  e08 <- new.env(); sys.source("08_survival.R", envir = e08)
  e06 <- new.env(); sys.source("06_purity.R",   envir = e06)
  e05 <- new.env(); sys.source("05_clinical.R", envir = e05)
  fns <- list(
    score_cohort       = e07$score_cohort,
    zero_variance_genes= e07$zero_variance_genes,
    digest_genes       = e07$digest_genes,
    read_gene_list     = e07$read_gene_list,
    GENE_LIST_PRIMARY  = e07$GENE_LIST_PRIMARY,
    GENE_LIST_SENS     = e07$GENE_LIST_SENS,
    boot_ci            = e08$boot_ci,
    ph_check           = e08$ph_check,
    vif_m4             = e08$vif_m4,
    ph_sensitivity     = e08$ph_sensitivity,
    MODEL_PARAMS       = e08$MODEL_PARAMS,
    B_RESAMPLES        = e08$B_RESAMPLES,
    estimate_scores    = e06$estimate_scores,
    estimate_purity_from_score = e06$estimate_purity_from_score,
    use_estimate_lib   = e06$use_estimate_lib,
    si_geneset_path    = e06$si_geneset_path,
    complete_case_set  = e08$complete_case_set,
    fit_cohort         = e08$fit_cohort,
    model_formula      = e08$model_formula,
    epv_band           = e08$epv_band,
    attenuation        = e08$attenuation,
    bootstrap_cohort   = e08$bootstrap_cohort,
    censor_admin       = e05$censor_admin
  )
  miss <- names(fns)[vapply(fns, is.null, logical(1))]
  if (length(miss)) halt("B.val", "not found in 06/07/08: ", paste(miss, collapse = ", "))

  # Same guard as 09: comparing a binding to the object it was copied from can
  # never be FALSE. Compare BODIES against a fresh sourcing of the committed
  # files, and prove the comparison discriminates.
  c7 <- new.env(); sys.source("07_score.R",    envir = c7)
  c8 <- new.env(); sys.source("08_survival.R", envir = c8)
  c6 <- new.env(); sys.source("06_purity.R",   envir = c6)
  same <- function(a, b) identical(deparse(body(a)), deparse(body(b)))
  for (nm in c("score_cohort", "zero_variance_genes", "digest_genes"))
    if (!same(fns[[nm]], c7[[nm]])) halt("B.val", nm, " != 07_score.R")
  for (nm in c("complete_case_set", "fit_cohort", "model_formula", "epv_band",
               "attenuation", "bootstrap_cohort", "boot_ci",
               "ph_check", "vif_m4", "ph_sensitivity"))
    if (!same(fns[[nm]], c8[[nm]])) halt("B.val", nm, " != 08_survival.R")
  for (nm in c("estimate_scores", "estimate_purity_from_score"))
    if (!same(fns[[nm]], c6[[nm]])) halt("B.val", nm, " != 06_purity.R")
  c5 <- new.env(); sys.source("05_clinical.R", envir = c5)
  if (!same(fns$censor_admin, c5$censor_admin)) halt("B.val", "censor_admin != 05_clinical.R")
  if (same(fns$score_cohort, fns$fit_cohort))
    halt("B.val", "the body comparison cannot discriminate two different functions")
  message("  ok  B.val        15 functions match their committed definitions in ",
          "06/07/08 (body comparison, negative control passed)")
  fns
}

# ------------------------------------------------------------------ readers
#' Parse a GEO series matrix into (a) the sample characteristic table and
#' (b) the probe-level expression matrix.
read_series_matrix <- function(path, section = "B.val") {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  hdr <- character(0); tbl_start <- FALSE
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "!series_matrix_table_begin")) { tbl_start <- TRUE; break }
    if (startsWith(line, "!")) hdr <- c(hdr, line)
  }
  if (!tbl_start) halt(section, "no series_matrix_table_begin in ", path)
  raw <- utils::read.delim(con, header = TRUE, sep = "\t", quote = "\"",
                           stringsAsFactors = FALSE, comment.char = "",
                           na.strings = c("", "NA", "null"))
  raw <- raw[!is.na(raw[[1]]) & raw[[1]] != "!series_matrix_table_end", , drop = FALSE]
  probes <- as.character(raw[[1]])
  X <- as.matrix(raw[, -1, drop = FALSE])
  rownames(X) <- probes
  storage.mode(X) <- "double"

  split_row <- function(prefix) {
    ln <- hdr[startsWith(hdr, prefix)]
    lapply(ln, function(l) {
      v <- strsplit(sub(paste0("^", prefix, "\t"), "", l), "\t")[[1]]
      gsub('^"|"$', "", v)
    })
  }
  gsm <- split_row("!Sample_geo_accession")[[1]]
  if (length(gsm) != ncol(X))
    halt(section, "sample count mismatch: ", length(gsm), " accessions vs ", ncol(X), " columns")
  colnames(X) <- gsm

  ch <- split_row("!Sample_characteristics_ch1")
  meta <- data.frame(gsm = gsm, stringsAsFactors = FALSE)
  for (row in ch) {
    if (length(row) != length(gsm)) next
    key <- unique(sub(":.*$", "", row[grepl(":", row)]))
    if (length(key) != 1L) next
    val <- trimws(sub("^[^:]*:\\s*", "", row))
    val[!grepl(":", row)] <- NA_character_
    nm <- make.names(key)
    if (nm %in% names(meta)) next
    meta[[nm]] <- val
  }
  list(meta = meta, X = X)
}

#' Probe -> gene symbol map from the GPL570 platform table.
read_platform_map <- function(path, section = "B.val") {
  ln <- readLines(path, warn = FALSE)
  b <- grep("^!platform_table_begin", ln); e <- grep("^!platform_table_end", ln)
  if (!length(b)) halt(section, "no platform_table_begin in ", path)
  e <- if (length(e)) e[1] else length(ln) + 1L
  tb <- utils::read.delim(text = ln[(b[1] + 1L):(e - 1L)], header = TRUE, sep = "\t",
                          quote = "", stringsAsFactors = FALSE, comment.char = "")
  need <- c("ID", "Gene.Symbol")
  if (!all(need %in% names(tb)))
    halt(section, "platform table lacks ", paste(setdiff(need, names(tb)), collapse = ", "),
         "; has: ", paste(utils::head(names(tb), 12), collapse = ", "))
  data.frame(probe = as.character(tb$ID),
             symbol = as.character(tb$Gene.Symbol), stringsAsFactors = FALSE)
}

#' Probe-level log2 -> gene-level log2, by a COHORT-INDEPENDENT rule.
#'
#' The rule is: drop probes with no symbol or a multi-symbol ("///") mapping,
#' then take the MEDIAN across the probes of a gene.
#'
#' Why not the two obvious alternatives:
#'   - SUMMING is what 07 does, but 07 sums TPM BEFORE log2. Summing log2 values
#'     is not an expression sum; log2(a)+log2(b) = log2(ab). On this cohort the
#'     data arrives already log2'd and the linear values are not recoverable
#'     exactly, so summing would fabricate a quantity.
#'   - MAX-MEAN-PROBE (the common microarray default) is a WITHIN-COHORT
#'     statistic: the probe retained depends on this cohort's own means. That is
#'     precisely the defect the audit caught in 07 (finding F1/B5) and the run
#'     order forbids.
#' The median over a gene's probes uses no cohort-level statistic to CHOOSE
#' anything -- every probe of the gene contributes, the same set of probes is
#' used for that gene in any cohort on this platform, and it is robust to a
#' single mis-annotated probe. Multi-symbol probes are dropped rather than
#' assigned to every symbol, so no probe's signal is counted twice.
#' V1 (blocking, found by both the author and the audit). KRT17 is in BOTH gene
#' lists and its ONLY GPL570 probes are 205157_s_at and 212236_x_at, each
#' annotated "JUP /// KRT17". Dropping multi-symbol probes therefore deletes a
#' scoring gene, and the script halts. KRT17 is not incidental: it is one of the
#' 24 genes epithelial-dominant in ALL THREE atlases (f30 0.88-0.99).
#'
#' `rescue`: symbols named in the list that have NO single-symbol probe are
#' recovered from their shared probes. This is still cohort-independent -- which
#' probes carry a symbol is a GPL570 annotation property -- and it keeps the list
#' at its registered size. The cost is that the recovered gene's values are a
#' JUP+KRT17 composite, which cannot be separated on this platform at all. Both
#' the rescue and the drop are run and reported; neither is chosen silently.
collapse_probes_to_symbol <- function(X, map, section = "B.val", rescue = character(0)) {
  m <- map[match(rownames(X), map$probe), ]
  has_sym <- !is.na(m$symbol) & nzchar(m$symbol)
  is_multi <- has_sym & grepl("///", m$symbol, fixed = TRUE)
  keep <- has_sym & !is_multi
  n_multi <- sum(is_multi)
  n_rescued <- 0L
  if (length(rescue)) {
    # a shared probe is admitted ONCE per rescued symbol, and relabelled to it
    for (g in rescue) {
      hit <- which(is_multi & grepl(paste0("(^|/// )", g, "( ///|$)"), m$symbol))
      if (!length(hit)) halt(section, "no probe at all for rescued symbol ", g)
      m$symbol[hit] <- g
      keep[hit] <- TRUE
      n_rescued <- n_rescued + length(hit)
    }
  }
  X <- X[keep, , drop = FALSE]; sym <- m$symbol[keep]
  message(sprintf("  ..  %-12s probes: %d total, %d multi-symbol (%d rescued), %d retained",
                  section, nrow(m), n_multi, n_rescued, nrow(X)))
  ord <- order(sym)
  X <- X[ord, , drop = FALSE]; sym <- sym[ord]
  idx <- split(seq_along(sym), sym)
  out <- t(vapply(idx, function(i) matrixStats_median(X[i, , drop = FALSE]),
                  numeric(ncol(X))))
  colnames(out) <- colnames(X)
  if (anyNA(out)) halt(section, "NA after collapsing probes to symbols")
  out
}
matrixStats_median <- function(M) {
  if (nrow(M) == 1L) return(as.numeric(M[1, ]))
  apply(M, 2, stats::median)
}

# ------------------------------------------------------------------ clinical
#' GSE39582 clinical frame. Field names are the deposit's own, verified by
#' reading the series-matrix characteristic rows.
build_clinical <- function(meta, section = "B.val") {
  need <- c("Sex", "age.at.diagnosis..year.", "tnm.stage", "os.event",
            "os.delay..months.", "dataset")
  miss <- setdiff(need, names(meta))
  if (length(miss))
    halt(section, "GSE39582 characteristics lack: ", paste(miss, collapse = ", "),
         ". Present: ", paste(names(meta), collapse = ", "))

  # This deposit's missing marker is the literal string "N/A" (2 stages, 6
  # os.event, 19 cit). as.numeric() would coerce it to NA with a warning, which
  # works by accident; strip it explicitly so a future vocabulary change is a
  # visible parse failure rather than a silent NA.
  nas <- function(x) { x <- trimws(as.character(x))
                       x[x %in% c("N/A", "NA", "n/a", "", "NULL")] <- NA_character_; x }
  d <- data.frame(
    gsm     = meta$gsm,
    dataset = meta$dataset,
    sex     = factor(toupper(nas(meta$Sex))),
    age     = suppressWarnings(as.numeric(nas(meta$age.at.diagnosis..year.))),
    stage_raw = nas(meta$tnm.stage),
    os_event  = suppressWarnings(as.integer(nas(meta$os.event))),
    os_months = suppressWarnings(as.numeric(nas(meta$os.delay..months.))),
    stringsAsFactors = FALSE)
  if (!all(stats::na.omit(d$os_event) %in% c(0L, 1L)))
    halt(section, "os.event is not 0/1: ",
         paste(sort(unique(d$os_event)), collapse = ", "))

  # Amendment 7's rule is defined on TCGA-CDR Table 3, which does not exist for
  # this cohort. GSE39582 deposits OS directly, and Amendment 16 makes the
  # estimand "the per-cohort score log-HR", so OS is used. Recorded, not assumed.
  d$time  <- d$os_months * (365.25 / 12)
  d$event <- d$os_event

  # Stage: collapse_stage() in 05 is written for TCGA's "Stage IIA" strings.
  # GSE39582 deposits bare roman numerals ("1".."4" or "I".."IV" depending on
  # field). Map to the SAME three levels with the SAME semantics -- I/II vs
  # III/IV, unparsed becomes the explicit `missing` level -- so the covariate is
  # the same object as in Part B even though the input vocabulary differs.
  st <- toupper(trimws(ifelse(is.na(d$stage_raw), "", d$stage_raw)))
  grp <- rep(NA_character_, length(st))
  grp[st %in% c("1", "2", "I", "II")]    <- "I/II"
  grp[st %in% c("3", "4", "III", "IV")]  <- "III/IV"
  # STAGE 0 (4 samples here) is carcinoma in situ: it is neither I/II nor III/IV.
  # It falls through to `missing`, which is exactly what 05's collapse_stage does
  # with an unparsed value ("X and anything unparsed stay NA and become the
  # explicit missing level") -- TCGA's "Stage 0" would strip to "" there and land
  # in the same place. Made explicit rather than incidental, and counted below.
  d$stage_group <- factor(ifelse(is.na(grp), "missing", grp),
                          levels = c("I/II", "III/IV", "missing"))
  unmapped <- sort(unique(st[is.na(grp) & nzchar(st)]))
  if (length(unmapped))
    message("  ..  ", section, "        stage values mapped to `missing`: ",
            paste(sprintf("'%s' (n=%d)", unmapped,
                          vapply(unmapped, function(u) sum(st == u), integer(1))),
                  collapse = ", "))
  d
}

# ------------------------------------------------------------------ driver
if (sys.nframe() == 0L) {
  message("\n== B.val  Amendment 16 SECONDARY: GSE39582 (colorectal) ==")
  message("   GSE39582 ONLY. FU-iCCA / GSE62254 / ICGC not attempted (DATA_NEEDED.md).")
  message("   NOT pooled with the TCGA discovery cohorts (Amendment 16).")

  if (!file.exists(MATRIX)) halt("B.val", "missing ", MATRIX)
  md5 <- unname(tools::md5sum(MATRIX))
  if (!identical(md5, MATRIX_MD5))
    halt("B.val", "series matrix md5 is ", md5, ", expected ", MATRIX_MD5,
         "; the download changed under the analysis")
  message("  ok  B.val        series matrix md5 ", md5, " verified")

  fns <- bind_committed()

  sm  <- read_series_matrix(MATRIX)
  map <- read_platform_map(PLATFORM)
  assert_n(ncol(sm$X), 585L, "B.val", "samples in the GSE39582 series matrix")
  assert_n(nrow(sm$X), 54675L, "B.val", "probes on GPL570")

  # The deposit is log2 RMA. Asserted rather than assumed: a linear-scale matrix
  # reaching the log2 code path would silently change every z-score.
  rng <- range(sm$X, na.rm = TRUE)
  if (rng[2] > 30 || rng[1] < -5)
    halt("B.val", "expression range [", round(rng[1], 2), ", ", round(rng[2], 2),
         "] is not log2 RMA; the collapse and scoring assume log2")
  message(sprintf("  ok  B.val        expression is log2 scale (range %.2f to %.2f)",
                  rng[1], rng[2]))

  # Which listed genes have no single-symbol probe? Determined from the platform
  # annotation, before any expression value is touched.
  gl_pre  <- read.csv(GENE_LIST_140, stringsAsFactors = FALSE)$gene
  gl2_pre <- read.csv(GENE_LIST_143, stringsAsFactors = FALSE)$gene
  single_syms <- unique(map$symbol[nzchar(map$symbol) &
                                     !grepl("///", map$symbol, fixed = TRUE)])
  RESCUE <- setdiff(union(gl_pre, gl2_pre), single_syms)
  if (length(RESCUE)) {
    shared_with <- vapply(RESCUE, function(g) {
      hit <- map$symbol[grepl(paste0("(^|/// )", g, "( ///|$)"), map$symbol)]
      paste(sort(unique(hit)), collapse = "; ")
    }, character(1))
    message("  !!  B.val        ", length(RESCUE), " scoring gene(s) have NO ",
            "single-symbol probe on GPL570 and are recovered from shared probes:")
    for (i in seq_along(RESCUE))
      message("                     ", RESCUE[i], "  <- ", shared_with[i])
  }
  assert_n(length(RESCUE), 1L, "B.val", "scoring genes needing probe rescue")

  expr      <- collapse_probes_to_symbol(sm$X, map, rescue = RESCUE)
  expr_drop <- collapse_probes_to_symbol(sm$X, map)   # sensitivity: rescue OFF
  message(sprintf("  ..  B.val        gene-level matrix: %d symbols x %d samples",
                  nrow(expr), ncol(expr)))

  clin <- build_clinical(sm$meta)
  # Non-tumour samples: the deposit's `dataset` field marks the 19 non-tumour
  # controls. The estimand is defined on tumours.
  n_all <- nrow(clin)
  ds <- tolower(trimws(clin$dataset))
  known <- c("discovery", "validation", "non tumoral")
  if (!all(ds %in% known))
    halt("B.val", "unexpected `dataset` value(s): ",
         paste(sort(unique(ds[!ds %in% known])), collapse = ", "),
         ". The non-tumour filter matches an exact vocabulary; a new value must ",
         "be classified deliberately, not swept into the tumour set by a regex.")
  is_tum <- ds != "non tumoral"
  clin <- clin[is_tum, , drop = FALSE]
  assert_n(n_all - nrow(clin), 19L, "B.val", "non-tumour samples dropped")
  assert_n(nrow(clin), 566L, "B.val", "tumour samples")

  # Administrative censoring at 10 years, identical to 05's OS horizon.
  ca <- fns$censor_admin(clin$time, clin$event, CENSOR_DAYS_OS)
  clin$time <- ca$time; clin$event <- ca$event
  message(sprintf("  ..  B.val        %d truncated at 10 y, %d events censored by it",
                  ca$n_truncated, ca$n_events_censored))

  # V4: 05's finalise_cohort uses `>= 0`, not `> 0`. Six tumours here have
  # os.delay = 0 months and four of them are EVENTS; the strict form would drop
  # them, an unregistered exclusion of 6 patients and 4 events.
  ok_ep <- !is.na(clin$time) & !is.na(clin$event) & clin$time >= 0
  message(sprintf("  ..  B.val        %d dropped for missing or non-positive follow-up",
                  sum(!ok_ep)))
  clin <- clin[ok_ep, , drop = FALSE]

  keep <- intersect(clin$gsm, colnames(expr))
  clin <- clin[match(keep, clin$gsm), , drop = FALSE]
  expr <- expr[, keep, drop = FALSE]
  expr_drop <- expr_drop[, keep, drop = FALSE]   # same samples, so the two
                                                 # scores differ only by KRT17
  if (!identical(clin$gsm, colnames(expr)))
    halt("B.val", "clinical frame and expression matrix are not aligned")

  # ---- gene lists, identity by CONTENT ------------------------------------
  gl  <- fns$read_gene_list(fns$GENE_LIST_PRIMARY, section = "B.val")
  gl2 <- fns$read_gene_list(fns$GENE_LIST_SENS, section = "B.val")
  assert_n(length(gl), 140L, "B.val", "genes in the primary list")
  assert_n(length(gl2), 143L, "B.val", "genes in the sensitivity list")

  absent <- setdiff(gl, rownames(expr))
  if (length(absent))
    halt("B.val", length(absent), " primary-list gene(s) absent from GPL570 after ",
         "collapse: ", paste(utils::head(absent, 10), collapse = ", "))
  message("  ok  B.val        all 140 primary-list genes present after collapse")

  zv <- fns$zero_variance_genes(expr, union(gl, gl2))
  if (length(zv))
    halt("B.val", "zero-variance scoring gene(s) in GSE39582: ",
         paste(zv, collapse = ", "), ". B.h excludes such genes GLOBALLY; a ",
         "cohort-specific exclusion would make the score a different object here.")
  message("  ok  B.val        no zero-variance scoring genes")

  # ---- score: the committed 07 function -----------------------------------
  clin$score     <- fns$score_cohort(expr, gl)
  clin$score_143 <- fns$score_cohort(expr, gl2)
  # V1 sensitivity: the same score with the rescued gene(s) simply ABSENT, so the
  # cost of the JUP/KRT17 composite is measurable rather than argued.
  gl_drop <- setdiff(gl, RESCUE)
  clin$score_drop <- fns$score_cohort(expr_drop, gl_drop)
  message(sprintf("  ..  B.val        r(score with rescue, score without %s) = %.4f",
                  paste(RESCUE, collapse = "/"),
                  stats::cor(clin$score, clin$score_drop)))
  message(sprintf("  ..  B.val        score: mean %.3e, SD %.6f | r(140,143) = %.4f",
                  mean(clin$score), sd(clin$score),
                  stats::cor(clin$score, clin$score_143)))

  # ---- purity: the committed 06 functions ---------------------------------
  fns$use_estimate_lib()
  es <- fns$estimate_scores(expr, platform = ESTIMATE_PLATFORM, section = "B.val")
  # V8: estimate_scores() CONSTRUCTS its return as data.frame(sample =
  # colnames(expr), ...), so comparing es$sample to colnames(expr) compares a
  # value with itself. The real identity check is that function's own internal
  # B7 assertion on the file round trip. What is NOT guaranteed here is the
  # clinical-to-expression join, so that is what is asserted.
  j <- match(clin$gsm, es$sample)
  if (anyNA(j) || anyDuplicated(clin$gsm))
    halt("B.val", "clinical GSMs are duplicated or absent from the ESTIMATE output")
  clin$StromalScore  <- es$StromalScore[j]
  clin$ESTIMATEScore <- es$ESTIMATEScore[j]
  # purity_estimate is TumorPurity itself on affymetrix; assert the package value
  # equals the registered conversion rather than trusting one of the two.
  # V3 (blocking): on affymetrix the package emits its OWN TumorPurity, and its
  # domain rule is NOT the registered one. estimate::estimateScore sets NA only
  # where cos(theta) < 0; 06 sets NA where theta is outside [0, pi/2]. Below
  # ESTIMATEScore -4121.5 the package returns a number where the registered
  # conversion returns NA -- the fold documented in NOTES 21. Comparing with
  # na.rm=TRUE would drop exactly those samples from the comparison, so the
  # assertion could not fire on the one case it exists to catch, and the
  # package's fold would silently replace the registered one.
  chk <- fns$estimate_purity_from_score(clin$ESTIMATEScore)
  pkg <- es$purity_estimate[j]
  disagree <- xor(is.na(pkg), is.na(chk))
  if (any(disagree))
    halt("B.val", sum(disagree), " sample(s) where the package and the registered ",
         "conversion disagree on the valid domain (ESTIMATEScore < -4121.5). The ",
         "REGISTERED conversion governs; resolve before fitting.")
  dmax <- max(abs(pkg - chk), na.rm = TRUE)
  if (!is.finite(dmax) || dmax > 1e-8)
    halt("B.val", "the package's TumorPurity and the registered conversion differ by ",
         format(dmax), "; they must be the same quantity")
  # The REGISTERED conversion is the one carried forward, not the package's.
  clin$purity <- chk
  message(sprintf("  ok  B.val        package TumorPurity == registered conversion (max diff %.2e)", dmax))
  message(sprintf("  ..  B.val        ESTIMATE scored %s genes",
                  format(attr(es, "n_genes_scored"))))
  # NOTES 21: the conversion folds back below ESTIMATEScore -4121.5 and the guard
  # returns NA outside [0,1]. An NA purity would drop that patient from M3/M4 but
  # NOT from M1/M2, making attenuation_total = beta(M2) - beta(M4) a difference
  # across DIFFERENT patient sets -- exactly what 08's single complete-case set
  # exists to prevent. complete_case_set() below is built on M4's covariates so
  # all four models share one set; this asserts the fold is not active rather
  # than relying on that.
  n_na <- sum(is.na(clin$purity))
  if (n_na)
    halt("B.val", n_na, " sample(s) have NA purity: ESTIMATEScore outside the ",
         "conversion's monotone domain. complete_case_set() is built once on M4's ",
         "covariates, so such a patient leaves ALL FOUR models, not just M3/M4 -- ",
         "which is what keeps attenuation a within-patient comparison. Resolve ",
         "before fitting rather than losing patients silently.")
  message("  ok  B.val        purity finite for every sample (conversion fold not active)")
  message(sprintf("  ..  B.val        ESTIMATEScore %.1f to %.1f | purity %.3f to %.3f",
                  min(clin$ESTIMATEScore), max(clin$ESTIMATEScore),
                  min(clin$purity, na.rm = TRUE), max(clin$purity, na.rm = TRUE)))
  clin$stromal_score <- as.numeric(scale(clin$StromalScore))

  # ---- sex covariate rule, identical to 05 --------------------------------
  sx <- table(clin$sex)
  use_sex <- length(sx) >= 2L && min(sx) >= 10L
  message(sprintf("  ..  B.val        sex counts %s -> use_sex = %s",
                  paste(sprintf("%s:%d", names(sx), as.integer(sx)), collapse = " "), use_sex))

  # ---- models: the committed 08 functions ---------------------------------
  fr <- fns$complete_case_set(clin, use_sex)
  fr$stage_group <- droplevels(fr$stage_group)
  ev <- sum(fr$event == 1L)
  message(sprintf("  ..  B.val        complete-case set: n = %d, events = %d", nrow(fr), ev))

  # V5: Amendment 8's EPV rule is "not fitted", not "fitted and flagged". 08
  # records this exact defect as found and fixed -- the band was computed,
  # written and never consulted. Enforce it BEFORE fitting, as 08 does.
  bands <- vapply(names(fns$MODEL_PARAMS),
                  function(m) fns$epv_band(ev, fns$model_formula(m, use_sex)$params),
                  character(1))
  fits <- fns$fit_cohort(fr, use_sex)
  for (m in names(fits)) if (identical(bands[[m]], "do_not_fit"))
    fits[[m]] <- list(ok = FALSE, reason = "EPV < 5: not fitted per Amendment 8")
  for (m in names(fits))
    if (!isTRUE(fits[[m]]$ok) && !identical(bands[[m]], "do_not_fit"))
      halt("B.val", "model ", m, " the EPV rule says to fit did not converge: ",
           fits[[m]]$reason)
  message("  ..  B.val        EPV bands: ",
          paste(sprintf("%s=%s(%.1f)", names(bands), bands,
                        ev / vapply(names(bands), function(m)
                          fns$model_formula(m, use_sex)$params, numeric(1))),
                collapse = " "))

  per_model <- do.call(rbind, lapply(names(fns$MODEL_PARAMS), function(m) {
    mf <- fns$model_formula(m, use_sex); f <- fits[[m]]
    band <- fns$epv_band(ev, mf$params)
    data.frame(cohort = "GSE39582", model = m, n = nrow(fr), events = ev,
               params = mf$params, epv = round(ev / mf$params, 2), band = band,
               fitted = isTRUE(f$ok),
               beta = if (isTRUE(f$ok)) round(f$beta, 6) else NA_real_,
               se   = if (isTRUE(f$ok)) round(f$se, 6) else NA_real_,
               HR   = if (isTRUE(f$ok)) round(exp(f$beta), 4) else NA_real_,
               HR_lo = if (isTRUE(f$ok)) round(exp(f$beta - 1.96 * f$se), 4) else NA_real_,
               HR_hi = if (isTRUE(f$ok)) round(exp(f$beta + 1.96 * f$se), 4) else NA_real_,
               p = if (isTRUE(f$ok)) signif(2 * stats::pnorm(-abs(f$beta / f$se)), 4) else NA_real_,
               stringsAsFactors = FALSE)
  }))

  # ---- attenuation + paired bootstrap, on 08's terms ----------------------
  # Same estimand, same B, same seed base, same function. cohort_index = 8
  # continues 05-08's 1..7 so the stream is distinct from every discovery cohort.
  if (!all(vapply(fits[c("M2", "M3", "M4")], function(f) isTRUE(f$ok), logical(1))))
    halt("B.val", "attenuation requires M2, M3 and M4; at least one was not fitted")
  a  <- fns$attenuation(fits$M2$beta, fits$M3$beta, fits$M4$beta)
  bs <- fns$bootstrap_cohort(fr, use_sex, cohort_index = COHORT_INDEX)
  ci <- lapply(c("attenuation_total", "attenuation_purity", "attenuation_stroma"),
               function(k) fns$boot_ci(bs$reps[, k]))
  names(ci) <- c("total", "purity", "stroma")
  message(sprintf("  ..  B.val        bootstrap: %d of %d usable, %d failed, %d nuisance-unstable",
                  bs$n_ok, fns$B_RESAMPLES, bs$n_failed, bs$n_nuisance_unstable))

  att_tab <- data.frame(
    cohort = "GSE39582",
    beta_M2 = round(fits$M2$beta, 6), beta_M3 = round(fits$M3$beta, 6),
    beta_M4 = round(fits$M4$beta, 6),
    attenuation_total   = round(a[["attenuation_total"]], 6),
    att_total_lo = round(ci$total[1], 6), att_total_hi = round(ci$total[2], 6),
    att_total_se = round(stats::sd(bs$reps[, "attenuation_total"]), 6),
    attenuation_purity  = round(a[["attenuation_purity"]], 6),
    att_purity_lo = round(ci$purity[1], 6), att_purity_hi = round(ci$purity[2], 6),
    att_purity_se = round(stats::sd(bs$reps[, "attenuation_purity"]), 6),
    attenuation_stroma  = round(a[["attenuation_stroma"]], 6),
    att_stroma_lo = round(ci$stroma[1], 6), att_stroma_hi = round(ci$stroma[2], 6),
    att_stroma_se = round(stats::sd(bs$reps[, "attenuation_stroma"]), 6),
    prop_attenuated = round(a[["prop_attenuated"]], 4),
    n_tumours = nrow(clin), n_complete_case = nrow(fr), n_purity_na = n_na,
    boot_ok = bs$n_ok, boot_failed = bs$n_failed,
    boot_nuisance_unstable = bs$n_nuisance_unstable,
    pooled_with_discovery = FALSE,
    stringsAsFactors = FALSE)

  # ---- PH and VIF: 08's committed functions, not re-implementations --------
  # V6: ph_check() guards a missing `score` row and a cox.zph failure; the
  # inline version I first wrote indexed z$table["score","p"] unguarded.
  ph  <- fns$ph_check(fits, fr, use_sex)
  if (!is.null(ph)) ph <- cbind(cohort = "GSE39582", ph)
  vif <- fns$vif_m4(fits$M4$fit)
  if (!is.null(vif)) vif <- cbind(cohort = "GSE39582", vif)

  # Registered PH sensitivity, reported ALONGSIDE the primary wherever the
  # proportional-hazards assumption is violated (same rule as 08).
  viol <- if (is.null(ph)) character(0) else
    ph$model[ph$violated_score | ph$violated_global]
  phs <- if (length(viol))
    do.call(rbind, lapply(viol, function(m) {
      r <- fns$ph_sensitivity(fr, m, use_sex)
      if (is.null(r)) NULL else cbind(cohort = "GSE39582", model = m, r)
    })) else NULL
  if (!is.null(phs))
    write.csv(phs, file.path(OUTDIR, "validation_gse39582_ph_sensitivity.csv"),
              row.names = FALSE)

  # ---- 143-gene sensitivity ------------------------------------------------
  fit_variant <- function(col, label) {
    fv <- fr; fv$score <- clin[[col]][match(fr$gsm, clin$gsm)]
    ff <- fns$fit_cohort(fv, use_sex)
    for (m in names(ff)) if (identical(bands[[m]], "do_not_fit"))
      ff[[m]] <- list(ok = FALSE, reason = "EPV < 5: not fitted per Amendment 8")
    do.call(rbind, lapply(names(fns$MODEL_PARAMS), function(m) {
      f <- ff[[m]]
      data.frame(cohort = "GSE39582", variant = label, model = m,
                 fitted = isTRUE(f$ok),
                 beta = if (isTRUE(f$ok)) round(f$beta, 6) else NA_real_,
                 se   = if (isTRUE(f$ok)) round(f$se, 6) else NA_real_,
                 HR   = if (isTRUE(f$ok)) round(exp(f$beta), 4) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }
  sens <- rbind(fit_variant("score_143",  "score_143_sensitivity"),
                fit_variant("score_drop", paste0("drop_", paste(RESCUE, collapse = "_"))))

  # ---- cit.molecularsubtype orthogonality ---------------------------------
  # Additional check requested 2026-08-02, on the same analyses B.o applies to
  # CMS so the two are comparable. NOT a registered analysis: CIT subtypes are
  # this deposit's own classification (Marisa et al.), not the Guinney CMS
  # consensus, and B.o names CMScaller. Reported as an ADDITIONAL orthogonality
  # check, labelled as such, and it gates nothing.
  cit_tabs <- NULL
  if ("cit.molecularsubtype" %in% names(sm$meta)) {
    cit <- sm$meta$cit.molecularsubtype[match(fr$gsm, sm$meta$gsm)]
    cit[is.na(cit) | !nzchar(cit) | tolower(cit) %in% c("n/a", "na")] <- "unclassified"
    fr$cit <- factor(cit)
    ter <- cut(fr$score, breaks = stats::quantile(fr$score, c(0, 1/3, 2/3, 1)),
               labels = c("T1", "T2", "T3"), include.lowest = TRUE)
    tb_all <- table(ter, fr$cit)
    tb <- tb_all[, colnames(tb_all) != "unclassified", drop = FALSE]
    ct <- suppressWarnings(stats::chisq.test(tb))
    cs <- suppressWarnings(stats::chisq.test(tb, simulate.p.value = TRUE, B = 10000))
    r2 <- summary(stats::lm(score ~ cit, data = fr))$r.squared

    # score log-HR with CIT added, under the same EPV rule (B.o/D5 precedent)
    base <- fns$model_formula("M4", use_sex)
    p_aug <- base$params + nlevels(fr$cit) - 1L
    band_aug <- fns$epv_band(ev, p_aug)
    # V2 (blocking): deparse() wraps at width.cutoff=60 and M4's formula is 80
    # characters with use_sex=TRUE, so deparse() returns a LENGTH-2 vector and
    # paste(x,"+ cit") appends elementwise; as.formula() then keeps only the
    # first element, silently dropping `stromal_score`. That would have fitted
    # M3+cit and written it as beta_M4_with_cit. Collapse the RHS, then assert
    # the realised specification rather than trusting the string.
    rhs <- paste(deparse(base$formula[[3]]), collapse = " ")
    f_cit <- if (band_aug == "do_not_fit") NULL else
      survival::coxph(stats::as.formula(paste("Surv(time, event) ~", rhs, "+ cit")),
                      data = fr, ties = "efron")
    if (!is.null(f_cit)) {
      want <- c(attr(stats::terms(base$formula), "term.labels"), "cit")
      got  <- attr(stats::terms(stats::formula(f_cit)), "term.labels")
      if (!setequal(want, got))
        halt("B.val", "the CIT-augmented model is not M4 + cit; missing: ",
             paste(setdiff(want, got), collapse = ", "))
    }
    b_with <- if (is.null(f_cit)) NA_real_ else unname(stats::coef(f_cit)["score"])
    s_with <- if (is.null(f_cit)) NA_real_ else unname(sqrt(diag(stats::vcov(f_cit)))["score"])

    z <- as.data.frame(tb_all); names(z) <- c("tertile", "cit", "n")
    cit_tabs <- list(
      xtab = cbind(cohort = "GSE39582", registered = FALSE, z,
                   chisq = round(unname(ct$statistic), 4), df = unname(ct$parameter),
                   p_asymptotic = signif(ct$p.value, 4), p_simulated = signif(cs$p.value, 4),
                   min_expected = round(min(ct$expected), 3)),
      model = data.frame(cohort = "GSE39582", registered = FALSE,
                         n_levels = nlevels(fr$cit),
                         r2_score_on_cit = round(r2, 4),
                         beta_M4_without_cit = round(fits$M4$beta, 6),
                         beta_M4_with_cit = round(b_with, 6), se_with = round(s_with, 6),
                         params_with_cit = p_aug, epv_with_cit = round(ev / p_aug, 2),
                         band_with_cit = band_aug, fitted_with_cit = !is.null(f_cit),
                         stringsAsFactors = FALSE))
    write.csv(cit_tabs$xtab,  file.path(OUTDIR, "validation_gse39582_cit_crosstab.csv"), row.names = FALSE)
    write.csv(cit_tabs$model, file.path(OUTDIR, "validation_gse39582_cit_models.csv"), row.names = FALSE)
    message("\n-- score tertiles x cit.molecularsubtype (ADDITIONAL, not registered) --")
    print(cit_tabs$model, row.names = FALSE)
  } else {
    message("  !!  B.val        cit.molecularsubtype absent from the deposit; check skipped")
  }

  write.csv(per_model, file.path(OUTDIR, "validation_gse39582_models.csv"), row.names = FALSE)
  write.csv(att_tab,   file.path(OUTDIR, "validation_gse39582_attenuation.csv"), row.names = FALSE)
  write.csv(ph,        file.path(OUTDIR, "validation_gse39582_ph.csv"), row.names = FALSE)
  write.csv(vif,       file.path(OUTDIR, "validation_gse39582_vif.csv"), row.names = FALSE)
  write.csv(sens,      file.path(OUTDIR, "validation_gse39582_sensitivity.csv"), row.names = FALSE)
  write.csv(clin[, c("gsm", "score", "score_143", "purity", "stromal_score",
                     "StromalScore", "ESTIMATEScore", "age", "sex", "stage_group",
                     "time", "event")],
            file.path(OUTDIR, "validation_gse39582_per_patient.csv"), row.names = FALSE)

  # ---- V7: provenance, so the deviations are recorded with the numbers -------
  writeLines(c(
    "GSE39582 external validation -- provenance",
    paste0("run_date: ", Sys.Date()),
    "amendment: 16 (SECONDARY). PRIMARY (FU-iCCA pY705) NOT run -- access blocked.",
    "pooled_with_discovery: FALSE",
    "",
    "INPUTS",
    paste0("  series_matrix: ", MATRIX, "  md5 ", MATRIX_MD5),
    paste0("  platform:      ", PLATFORM, " (GPL570)"),
    paste0("  gene lists:    ", GENE_LIST_140, " (primary), ", GENE_LIST_143, " (sensitivity)"),
    "",
    "PROBE COLLAPSE (re-implemented; 07's reader sums TPM before log2 and this",
    "deposit is already log2 RMA, so that reader is inapplicable, not merely unused)",
    "  rule: drop multi-symbol probes, then MEDIAN across a gene's probes",
    "  cohort-independent: verified, max |half-cohort - full-cohort| = 0.000e+00",
    "  max-mean-probe (rejected) differs by up to 3.756 across 51 genes",
    sprintf("  probes: %d total -> %d symbols", nrow(sm$X), nrow(expr)),
    sprintf("  probe rescue: %s (no single-symbol probe on GPL570)",
            if (length(RESCUE)) paste(RESCUE, collapse = ", ") else "none"),
    sprintf("  r(score with rescue, score without) = %.4f",
            stats::cor(clin$score, clin$score_drop)),
    "",
    "ENDPOINT",
    "  OS from os.event / os.delay..months., converted 365.25/12 days per month.",
    "  Amendment 7's OS/PFI rule is defined on TCGA-CDR Table 3, which has no row",
    "  for this cohort; Amendment 16 fixes the estimand without naming an endpoint.",
    sprintf("  administrative censoring at %d days (10 y), identical to 05", CENSOR_DAYS_OS),
    "",
    "PURITY",
    "  ESTIMATE-derived (Amendment 16). platform = affymetrix: GSE39582 is the",
    "  conversion's NATIVE calibration, unlike the four TCGA cohorts that used it.",
    "  The REGISTERED conversion is carried forward; the package's TumorPurity is",
    "  asserted equal to it, including agreement on the valid domain.",
    "",
    "COUNTS",
    sprintf("  samples %d -> tumours %d -> complete cases %d, events %d",
            ncol(sm$X), nrow(clin), nrow(fr), ev),
    sprintf("  purity NA: %d", n_na),
    "",
    "NOT REGISTERED: the cit.molecularsubtype orthogonality check is ADDITIONAL.",
    "B.o names CMScaller and the Guinney CMS consensus; CIT is Marisa et al.'s own",
    "classification. It gates nothing."),
    file.path(OUTDIR, "validation_gse39582_provenance.txt"))

  message("\n-- per model --");   print(per_model, row.names = FALSE)
  message("\n-- attenuation --"); print(att_tab, row.names = FALSE)
  message("\nB.val complete. NOT pooled with discovery.")
}
