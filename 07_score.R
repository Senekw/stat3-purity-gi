#!/usr/bin/env Rscript
# 07_score.R -- Part B: STAT3 activity score construction. NO model fitted.
#
# Score, exactly as registered (analysis_plan.md B.h):
#   1. primary tumour samples only, one aliquot per patient   [done in 05]
#   2. expression = log2(tpm_unstrand + 1)
#   3. subset to the final gene list; halt if any gene is missing
#   4. z-score EACH GENE WITHIN COHORT (mean 0, sd 1 across that cohort's patients)
#   5. score = mean across the final gene list
#   6. scale the resulting score to unit SD within cohort
#
# No dichotomisation, no tertiles, no data-derived cutpoint anywhere.
#
# STRUCTURE: pure functions taking a MATRIX, so external validation cohorts reuse
# this path unchanged. Driver guarded by `if (sys.nframe() == 0L)`.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

halt <- function(section, ...)
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)

assert_n <- function(observed, expected, section, what) {
  if (!identical(as.integer(observed), as.integer(expected)))
    halt(section, what, ": expected ", expected, ", observed ", observed)
  message(sprintf("  ok  %-12s %-54s = %d", section, what, observed))
  invisible(TRUE)
}

COHORTS <- c("TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-ESCA",
             "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL")

# The gene lists 04 locked. Identity is asserted from the FILE's own columns, not
# from its name -- final_gene_list.csv has already changed meaning once (143
# before Amendment 15, 140 after).
#
# LABEL DISCREPANCY, reported not silently reconciled: the run order asks for
# `list_id == "final_140"`. The value committed by 04_lock_gene_list.R (commit
# 744d84c, audited) is "primary_140", with the sensitivity list "sensitivity_143".
# The LIST is unambiguous and matches on every other check -- 140 rows, n_genes
# 140, and the same genes -- so this is a naming mismatch, not a data mismatch.
# Renaming the locked file's list_id to match the instruction would be editing a
# committed artefact to make an assertion pass, which is the one move the standing
# rules forbid. The assertion therefore uses the committed value and the
# discrepancy is reported.
# Content digest of the sorted symbols: identity that survives reordering and
# column changes but not a substituted gene. Values recorded from the artefacts
# 04_lock_gene_list.R committed at 744d84c.
digest_genes <- function(g) {
  tf <- tempfile(); on.exit(unlink(tf), add = TRUE)
  writeLines(sort(unique(as.character(g))), tf)
  unname(tools::md5sum(tf))
}

GENE_LIST_PRIMARY <- list(path = "data/panel/final_gene_list_140.csv",
                          list_id = "primary_140",     n = 140L,
                          instructed_list_id = "final_140",
                          digest = "eb167e8c7a33b4202bd609a17defa629")
GENE_LIST_SENS    <- list(path = "data/panel/final_gene_list_143.csv",
                          list_id = "sensitivity_143", n = 143L,
                          digest = "8a54835eedcfe3b23dcb56ce74805a87")

# B.h: the stromal subscore is built from panel genes that are NOT
# epithelial-dominant. If it correlates with the main score above this, the
# registered ESTIMATE stromal fallback applies.
STROMAL_COLLINEARITY_THRESHOLD <- 0.9

# Registered tolerance for zero variance (analysis_plan.md:662-663), NOT exact
# equality: a gene with sd = 4.5e-12 would otherwise pass and be divided through,
# producing a column of ~1e12 that dominates the mean across genes.
ZERO_VAR_TOL <- 1e-8

# ============================================================ PURE FUNCTIONS

#' Read a locked gene list and assert its identity from its own contents.
read_gene_list <- function(spec, section = "B.h") {
  if (!file.exists(spec$path)) halt(section, "gene list not found: ", spec$path)
  d <- read.csv(spec$path, stringsAsFactors = FALSE)
  for (col in c("list_id", "n_genes", "gene"))
    if (!col %in% names(d))
      halt(section, spec$path, " lacks the '", col, "' column; it is not a locked list")
  if (!identical(unique(d$list_id), spec$list_id))
    halt(section, spec$path, ": list_id is '", paste(unique(d$list_id), collapse = ","),
         "', expected '", spec$list_id, "'")
  # The Part B run order named a different string for the primary list
  # ("final_140"). That instructed label was erroneous -- written from memory --
  # and the locked artefact is canonical (confirmed 2026-08-02). Recorded as a
  # one-line note, not a warning: there is nothing to decide. Identity is enforced
  # by the sorted-symbol digest below, which is what stop condition 6 rests on.
  if (!is.null(spec$instructed_list_id) &&
      !identical(spec$instructed_list_id, spec$list_id))
    message("  ..  ", section, "          note: an earlier run order said list_id == '",
            spec$instructed_list_id, "'; the canonical artefact value is '",
            spec$list_id, "' (identity enforced by symbol digest)")
  if (!identical(unique(as.integer(d$n_genes)), spec$n))
    halt(section, spec$path, ": n_genes is ", paste(unique(d$n_genes), collapse = ","),
         ", expected ", spec$n)
  if (nrow(d) != spec$n)
    halt(section, spec$path, ": ", nrow(d), " rows but n_genes says ", spec$n)
  if (anyDuplicated(d$gene)) halt(section, spec$path, ": duplicate gene symbols")
  # F6: row count and label cannot detect the WRONG 140 genes -- a stale copy or a
  # truncate-and-append passes all three. Assert the content digest, recorded from
  # the artefact 04 committed at 744d84c.
  dg <- digest_genes(d$gene)
  if (!identical(dg, spec$digest))
    halt(section, spec$path, ": gene-set digest is ", dg, ", expected ", spec$digest,
         ". The file has the right shape but not the right genes.")
  message(sprintf("  ok  %-12s %-54s = %d", section,
                  paste0("gene list '", spec$list_id, "' verified"), nrow(d)))
  d$gene
}

#' log2(TPM + 1) matrix for given barcodes, genes x samples, symbols as rownames.
#'
#' Duplicate symbols are collapsed by MAXIMUM ROW MEAN before subsetting, so gene
#' selection is by name against a matrix in which each symbol appears once. Taking
#' whichever row happened to come first would make the score depend on annotation
#' order.
expression_log2tpm <- function(se, barcodes, section = "B.h") {
  if (!"tpm_unstrand" %in% SummarizedExperiment::assayNames(se))
    halt(section, "assay 'tpm_unstrand' is absent; available: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  j <- match(barcodes, as.character(SummarizedExperiment::colData(se)$barcode))
  if (any(is.na(j)))
    halt(section, "barcode(s) from the analysis set absent from the cohort: ",
         paste(utils::head(barcodes[is.na(j)], 3), collapse = ", "))
  x <- SummarizedExperiment::assay(se, "tpm_unstrand")[, j, drop = FALSE]
  if (any(x < 0, na.rm = TRUE)) halt(section, "negative TPM values")
  sym <- as.character(SummarizedExperiment::rowData(se)$gene_name)
  ok <- !is.na(sym) & nzchar(sym)
  x <- x[ok, , drop = FALSE]; sym <- sym[ok]
  # F1 (blocking): duplicate symbols are SUMMED on the TPM scale, before log2.
  # An earlier revision kept the row with the largest mean IN THIS COHORT, which
  # made the retained transcript cohort-dependent -- gene G could be scored from
  # one Ensembl row in COAD and another in STAD. B.h forbids exactly that: "The
  # score is computed on a single final gene list used identically in every
  # cohort." Summing is additively correct for TPM and cohort-invariant. It must
  # happen BEFORE log2, since log2(a+1)+log2(b+1) != log2(a+b+1).
  x <- rowsum(x, group = sym, reorder = TRUE)
  if (anyNA(x)) halt(section, "NA in the expression matrix after collapsing symbols")
  x <- log2(x + 1)
  colnames(x) <- barcodes
  x
}

#' Genes with zero variance WITHIN a cohort, which cannot be z-scored.
zero_variance_genes <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  # F4: B.h registers the tolerance explicitly -- "Zero variance is evaluated on
  # log2 TPM across each cohort's primary-tumour samples, using the tolerance
  # `sd > 1e-8` rather than exact equality." (analysis_plan.md:662-663)
  s <- apply(expr[g, , drop = FALSE], 1, sd, na.rm = TRUE)
  g[!is.finite(s) | s <= ZERO_VAR_TOL]
}

#' The registered score for one cohort.
#'
#' @param expr log2(TPM+1), genes x samples
#' @param genes the scoring gene set, already excluding globally-dropped genes
#'
#' Steps 3-6 of B.h. Pure: the z-scoring is WITHIN this matrix, so the function is
#' a closed transform of its arguments -- which is exactly why the validation
#' cohorts can reuse it.
score_cohort <- function(expr, genes, section = "B.h") {
  missing <- setdiff(genes, rownames(expr))
  if (length(missing))
    halt(section, length(missing), " scoring gene(s) absent from the expression ",
         "matrix: ", paste(utils::head(missing, 10), collapse = ", "),
         ". The list was built to be present in every cohort.")
  x <- expr[genes, , drop = FALSE]
  mu <- rowMeans(x, na.rm = TRUE)
  sdv <- apply(x, 1, sd, na.rm = TRUE)
  if (any(!is.finite(sdv) | sdv <= ZERO_VAR_TOL))
    halt(section, "zero-variance gene(s) reached score_cohort: ",
         paste(utils::head(genes[!is.finite(sdv) | sdv <= ZERO_VAR_TOL], 5), collapse = ", "),
         ". These must be excluded GLOBALLY before scoring, never per cohort.")
  z <- (x - mu) / sdv                       # 4. z-score each gene within cohort
  s <- colMeans(z, na.rm = TRUE)            # 5. mean across the gene list
  s <- s / sd(s)                            # 6. unit SD within cohort
  if (!isTRUE(all.equal(sd(s), 1))) halt(section, "score SD is not 1 after scaling")
  s
}

#' Stromal subscore (B.h): panel genes that are NOT epithelial-dominant.
#'
#' Dominance is read from Part A's committed dominance matrix. A gene is
#' epithelial-dominant if it was called dominant in at least two atlases -- the
#' same two-of-three rule that defines k, so the subscore's complement is defined
#' identically to the quantity Part A reported.
stromal_gene_set <- function(dominance_csv, panel_genes, section = "B.h") {
  dm <- read.csv(dominance_csv, stringsAsFactors = FALSE)
  p <- dm[which(dm$in_panel), , drop = FALSE]
  # F3: count evaluability and dominance separately. `sum(na.rm=TRUE)` alone
  # collapses "evaluated, not epithelial-dominant" and "never evaluable" into the
  # same bucket, so a gene with no compartment evidence would be asserted as
  # non-dominant and admitted to the stromal set. A gene evaluable in fewer than
  # two atlases cannot satisfy or refute a two-of-three rule and is excluded from
  # both sides. (Realised: 0 such genes -- every panel gene is evaluable in >= 2.)
  nd <- tapply(p$dominant, p$gene, function(v) sum(v, na.rm = TRUE))
  ne <- tapply(p$dominant, p$gene, function(v) sum(!is.na(v)))
  dom        <- names(nd)[nd >= 2]
  undecidable <- names(ne)[ne < 2]
  if (length(undecidable))
    message("  ..  B.h          ", length(undecidable), " gene(s) evaluable in < 2 ",
            "atlases, excluded from both the dominant and stromal sets: ",
            paste(undecidable, collapse = ", "))
  setdiff(setdiff(panel_genes, dom), undecidable)
}

#' ESTIMATE stromal score as the registered fallback when the subscore is
#' collinear with the main score.
estimate_stromal_fallback <- function(purity_csv, cohort, barcodes, section = "B.h") {
  p <- read.csv(purity_csv, stringsAsFactors = FALSE)
  p <- p[p$cohort == cohort, , drop = FALSE]
  v <- p$StromalScore[match(barcodes, p$barcode)]
  if (any(is.na(v)))
    halt(section, cohort, ": ESTIMATE StromalScore missing for ", sum(is.na(v)),
         " sample(s); the fallback cannot be applied")
  as.numeric(scale(v))
}

# ==================================================================== driver
if (sys.nframe() == 0L) {

  OUTDIR <- "output"; dir.create(OUTDIR, showWarnings = FALSE)
  message("\n== B.h  score construction (no model fitted) ==")

  clin <- read.csv(file.path(OUTDIR, "clinical_analysis_set.csv"), stringsAsFactors = FALSE)
  if (!nrow(clin)) halt("B.h", "clinical analysis set is empty; run 05 first")
  genes140 <- read_gene_list(GENE_LIST_PRIMARY)
  genes143 <- read_gene_list(GENE_LIST_SENS)
  panel <- read.csv("data/panel/panel_locked.csv", stringsAsFactors = FALSE)$gene
  if (!all(genes140 %in% genes143))
    halt("B.h", "the 140-gene primary list is not a subset of the 143-gene ",
         "sensitivity list; one of the locked artefacts is wrong")

  # F2 -- REGISTERED TEXT vs AUDIT OBJECTION, both computed, decision deferred.
  #
  # analysis_plan.md:718-719 reads: "The stromal subscore used as a covariate in
  # B.j model 4 is constructed identically from panel genes that are not
  # epithelial-dominant." PANEL is the locked 152, so the registered wording gives
  # the 106-gene set and that is what `stromal_score` carries.
  #
  # The Implementation Auditor objects that this readmits 9 of the 12 genes that
  # prespecification rules 2 and 3.3 removed from the scoring list -- including
  # TIMP1, excluded ONLY for sitting on chrX yet detected in 98.4% of stromal
  # cells, and IL3RA/IL2RG/CSF2RA at 61%/57%/44%. Rule 3.3 exists to keep
  # sex-differential expression out of the models; entering through the covariate
  # re-opens that door.
  #
  # Both are computed. The 97-gene variant is carried alongside as
  # `stromal_score_140` so the choice can be made on the numbers rather than
  # re-derived. NOT resolved here: it changes a covariate in the primary model.
  stromal_genes_all <- stromal_gene_set("output/compartment_dominance_matrix.csv", panel)
  stromal_genes_140 <- stromal_gene_set("output/compartment_dominance_matrix.csv", genes140)

  # --- pass 1: expression per cohort, and the GLOBAL zero-variance set --------
  # Registered decision D: zero-variance genes are excluded GLOBALLY, never per
  # cohort, so every cohort scores an identical gene set. Computed as the UNION
  # over cohorts, since a gene invariant in any one cohort cannot be z-scored
  # there and dropping it only in that cohort would break the identity.
  expr <- lapply(COHORTS, function(cc) {
    cl <- clin[clin$cohort == cc, , drop = FALSE]
    se <- readRDS(sprintf("data/tcga/%s_se.rds", cc))
    x <- expression_log2tpm(se, cl$barcode)
    rm(se); gc(verbose = FALSE)
    x
  })
  names(expr) <- COHORTS

  # F5: keep per-cohort attribution before flattening. B.h requires the exclusion
  # be recorded as zero_variance_in_<cohort>, which the union alone destroys.
  zv_by <- lapply(COHORTS, function(cc)
    zero_variance_genes(expr[[cc]], Reduce(union, list(genes143, stromal_genes_all, stromal_genes_140))))
  names(zv_by) <- COHORTS
  zv <- sort(unique(unlist(zv_by)))
  if (length(zv))
    message("  !!  B.h          zero-variance genes excluded GLOBALLY (", length(zv), "): ",
            paste(zv, collapse = ", "))
  else
    message("  ..  B.h          no zero-variance genes in any cohort")

  g140 <- setdiff(genes140, zv)
  g143 <- setdiff(genes143, zv)
  gstr    <- setdiff(stromal_genes_all, zv)   # registered wording: panel 152 domain
  gstr140 <- setdiff(stromal_genes_140, zv)   # audit F2 variant: final 140 domain
  message(sprintf("  ..  B.h          scoring sets after global exclusion: primary %d, sensitivity %d, stromal %d (panel-152 domain), stromal_140 %d (final-140 domain)",
                  length(g140), length(g143), length(gstr), length(gstr140)))
  if (!length(gstr)) halt("B.h", "the stromal gene set is empty")

  # --- pass 2: scores --------------------------------------------------------
  out <- do.call(rbind, lapply(COHORTS, function(cc) {
    cl <- clin[clin$cohort == cc, , drop = FALSE]
    x <- expr[[cc]]
    s140 <- score_cohort(x, g140)
    s143 <- score_cohort(x, g143)
    sstr <- score_cohort(x, gstr)
    sstr140 <- score_cohort(x, gstr140)
    data.frame(cohort = cc, patient = cl$patient, barcode = cl$barcode,
               score = unname(s140), score_143 = unname(s143),
               stromal_score = unname(sstr), stromal_score_140 = unname(sstr140),
               stringsAsFactors = FALSE)
  }))

  # --- registered collinearity check and fallback ----------------------------
  coll <- do.call(rbind, lapply(COHORTS, function(cc) {
    d <- out[out$cohort == cc, ]
    data.frame(cohort = cc, r = cor(d$score, d$stromal_score),
               r_stromal140 = cor(d$score, d$stromal_score_140),
               r_between_variants = cor(d$stromal_score, d$stromal_score_140),
               stringsAsFactors = FALSE)
  }))
  # F7: `any(c(FALSE, NA))` is NA and `if (NA)` errors at the registered decision
  # point. An NA correlation is a defect, not a "no", so it halts explicitly.
  if (anyNA(coll$r))
    halt("B.h", "the score/stromal correlation is NA in ",
         paste(coll$cohort[is.na(coll$r)], collapse = ", "),
         "; the registered |r| > 0.9 decision cannot be evaluated on an unknown")
  coll$exceeds <- abs(coll$r) > STROMAL_COLLINEARITY_THRESHOLD
  fallback_triggered <- any(coll$exceeds)
  message("\n-- stromal subscore vs main score --")
  print(coll, row.names = FALSE)

  if (fallback_triggered) {
    message("  !!  B.h          |r| > 0.9 in ", sum(coll$exceeds), " cohort(s): ",
            paste(sub("^TCGA-", "", coll$cohort[coll$exceeds]), collapse = ", "),
            " -- applying the REGISTERED ESTIMATE stromal fallback in ALL cohorts, ",
            "so the stromal covariate has one definition throughout.")
    out$stromal_score_subscore <- out$stromal_score
    for (cc in COHORTS) {
      i <- out$cohort == cc
      out$stromal_score[i] <- estimate_stromal_fallback(
        file.path(OUTDIR, "purity_per_patient.csv"), cc, out$barcode[i])
    }
    out$stromal_source <- "ESTIMATE_StromalScore"
  } else {
    out$stromal_score_subscore <- out$stromal_score
    out$stromal_source <- "panel_non_dominant_genes"
  }

  # --- distribution and invariants -------------------------------------------
  dist <- do.call(rbind, lapply(COHORTS, function(cc) {
    d <- out[out$cohort == cc, ]
    data.frame(cohort = cc, n = nrow(d),
               mean = round(mean(d$score), 6), sd = round(sd(d$score), 6),
               min = round(min(d$score), 4), max = round(max(d$score), 4),
               sd_143 = round(sd(d$score_143), 6),
               r_140_143 = round(cor(d$score, d$score_143), 4),
               sd_stromal = round(sd(d$stromal_score), 6),
               stringsAsFactors = FALSE)
  }))
  if (any(abs(dist$sd - 1) > 1e-8)) halt("B.h", "within-cohort score SD is not 1")
  if (any(abs(dist$sd_143 - 1) > 1e-8)) halt("B.h", "143-gene score SD is not 1")
  assert_n(nrow(out), nrow(clin), "B.h", "scored patients == clinical rows")
  if (anyDuplicated(out$patient)) halt("B.h", "duplicate patient in the score table")

  write.csv(out,  file.path(OUTDIR, "scores_per_patient.csv"), row.names = FALSE)
  write.csv(dist, file.path(OUTDIR, "score_distribution.csv"), row.names = FALSE)
  write.csv(coll, file.path(OUTDIR, "stromal_collinearity.csv"), row.names = FALSE)
  writeLines(c(
    paste0("zero_variance_excluded_globally: ", length(zv)),
    if (length(zv)) paste0("  ", zv) else "  (none)",
    "# per-cohort attribution (B.h: excluded_reason = zero_variance_in_<cohort>)",
    unlist(lapply(COHORTS, function(cc)
      paste0("  zero_variance_in_", sub("^TCGA-", "", cc), ": ",
             if (length(zv_by[[cc]])) paste(zv_by[[cc]], collapse = ", ") else "(none)"))),
    paste0("list_id_discrepancy: run order says 'final_140'; the locked artefact ",
           "from 04 (744d84c) says 'primary_140'. The committed value is asserted ",
           "and the artefact is unmodified."),
    paste0("n_genes_primary_140_after_exclusion: ", length(g140)),
    paste0("n_genes_sensitivity_143_after_exclusion: ", length(g143)),
    paste0("n_genes_stromal_subscore_panel152_REGISTERED: ", length(gstr)),
    paste0("n_genes_stromal_subscore_final140_AUDIT_VARIANT: ", length(gstr140)),
    paste0("stromal_genes_readmitted_by_panel_domain: ",
           paste(sort(setdiff(gstr, gstr140)), collapse = ", ")),
    paste0("stromal_fallback_triggered: ", fallback_triggered),
    paste0("stromal_source: ", out$stromal_source[1])),
    file.path(OUTDIR, "score_gene_sets.txt"))

  message("\n-- score distribution per cohort --")
  print(dist, row.names = FALSE)
  message("\nB.h complete. HARD STOP: no model fitted.")
}
