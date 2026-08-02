#!/usr/bin/env Rscript
# 06_purity.R -- Part B: tumour purity per patient. NO score, NO model.
#
# ESTIMATE (Yoshihara et al. 2013) stromal/immune/ESTIMATE scores and derived
# purity, plus the Aran et al. 2015 consensus CPE table where available, and the
# registered switch rule between them.
#
# STRUCTURE: pure functions taking a MATRIX (genes x samples), so the external
# validation cohorts reuse this code path unchanged. The driver is guarded by
# `if (sys.nframe() == 0L)` -- sourcing this file yields the functions and runs
# nothing.
#
# TWO INPUT FACTS ESTABLISHED BEFORE WRITING, both recorded in NOTES_FOR_REVIEW:
#
#  1. data/manual/aran_purity.xlsx IS NOT ON DISK. analysis_plan.md section 3.4
#     anticipates exactly this ("aran_purity.xlsx is not on disk, so the
#     proportion of patients in each cohort with a CPE value is unknown") and
#     prespecifies the resolution: CPE primary where coverage >= 80%, ESTIMATE
#     primary otherwise, reported per cohort. Coverage is therefore 0% and
#     ESTIMATE is primary in every cohort BY THE REGISTERED RULE, not by choice.
#     The script still looks for the file and uses it if it appears.
#
#  2. The ESTIMATE package derives TumorPurity ONLY for platform="affymetrix".
#     Read from the installed source: the purity branch is inside
#     `if (platform != "affymetrix") {...} else {...}`, so on illumina it returns
#     the three scores and no purity. TCGA STAR-Counts is Illumina RNA-seq. The
#     conversion cos(0.6049872018 + 0.0001467884 * ESTIMATEScore) is applied here
#     EXPLICITLY, with its Affymetrix calibration recorded as a limitation --
#     see `estimate_purity_from_score()`.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# The authentic ESTIMATE package (v1.0.13, pure R) is installed in the system R
# library, not the project env -- it is not on CRAN and not in the configured
# conda channels. Pure-R packages port across minor R versions; the signature
# file's identity is asserted by md5 below rather than assumed.
#
# B6: the .libPaths() mutation is NOT executed at source() time. This file's
# contract is that sourcing defines functions and changes nothing, and rewriting
# the library search path of a sourcing process would violate it. The driver
# calls use_estimate_lib() explicitly; a caller that sources this file for its
# functions must do the same.
ESTIMATE_LIB <- "/opt/homebrew/lib/R/4.6/site-library"

#' Make the system-library `estimate` loadable WITHOUT shadowing other packages.
#'
#' Prepending the whole system library to .libPaths() makes every package resolve
#' there first -- including readxl, whose Homebrew build links a second OpenMP
#' runtime and aborts the process ("OMP: Error #15 ... libomp.dylib already
#' initialized") once the conda BLAS has initialised its own. Diagnosed by
#' bisection, not worked around: KMP_DUPLICATE_LIB_OK is documented as unsafe and
#' can "silently produce incorrect results", which is not acceptable for a
#' covariate feeding the primary model.
#'
#' Instead the path is APPENDED, so it is searched only for packages the project
#' environment does not provide -- which is `estimate` and nothing else.
use_estimate_lib <- function(lib = ESTIMATE_LIB) {
  if (dir.exists(lib) && !lib %in% .libPaths()) .libPaths(c(.libPaths(), lib))
  invisible(.libPaths())
}

halt <- function(section, ...)
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)

assert_n <- function(observed, expected, section, what) {
  if (!identical(as.integer(observed), as.integer(expected)))
    halt(section, what, ": expected ", expected, ", observed ", observed)
  message(sprintf("  ok  %-14s %-52s = %d", section, what, observed))
  invisible(TRUE)
}

COHORTS <- c("TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-ESCA",
             "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL")

# Registered switch rule (analysis_plan.md 3.4): CPE primary at coverage >= 80%
# of the cohort's patients, ESTIMATE-derived purity primary below that.
CPE_COVERAGE_THRESHOLD <- 0.80

CPE_PATH <- "data/manual/aran_purity.xlsx"

# Yoshihara et al. 2013, Nat Commun 4:2612 -- the published conversion, taken
# verbatim from the estimate package source rather than retyped from the paper.
PURITY_A <- 0.6049872018
PURITY_B <- 0.0001467884

# Identity of the signature file, so a different ESTIMATE release cannot silently
# change the scores.
SI_GENESET_MD5 <- "9628130acbd541d8a33f21b4ab1e864f"
SI_N_PER_SIG   <- 141L

# ============================================================ PURE FUNCTIONS

#' Locate and verify the ESTIMATE signature file.
si_geneset_path <- function(section = "B.purity") {
  p <- system.file("extdata", "SI_geneset.gmt", package = "estimate")
  if (!nzchar(p) || !file.exists(p))
    halt(section, "the ESTIMATE signature file SI_geneset.gmt was not found. ",
         "The estimate package (v1.0.13) must be available; it is not on CRAN.")
  md5 <- unname(tools::md5sum(p))
  if (!identical(md5, SI_GENESET_MD5))
    halt(section, "SI_geneset.gmt md5 is ", md5, ", expected ", SI_GENESET_MD5,
         ". This is a different ESTIMATE release; scores would not be comparable ",
         "to the registered method.")
  gmt <- strsplit(readLines(p), "\t")
  n <- vapply(gmt, function(x) length(x) - 2L, integer(1))
  if (!identical(sort(n), c(SI_N_PER_SIG, SI_N_PER_SIG)))
    halt(section, "signature sizes are ", paste(n, collapse = "/"),
         ", expected 141/141")
  p
}

#' Convert an ESTIMATE score to tumour purity.
#'
#' cos(a + b * ESTIMATEScore), Yoshihara et al. 2013.
#'
#' CALIBRATION LIMITATION, recorded rather than hidden: this conversion was fitted
#' on AFFYMETRIX arrays, and the estimate package computes it only when
#' platform="affymetrix". TCGA STAR-Counts is Illumina RNA-seq, so applying it
#' here is an extrapolation across platforms. It is applied because the plan
#' specifies "ESTIMATE-derived purity" and this is the published derivation; the
#' alternative -- inventing an RNA-seq recalibration -- would be worse. Values
#' outside [0,1] are set NA exactly as the package does, never clamped, since a
#' clamp would manufacture a boundary value from an out-of-range score.
#' DOMAIN GUARD (audit finding B1). cos(a + bS) is monotone-decreasing in S only
#' while the ANGLE theta = a + bS lies in [0, pi/2]; that is
#' S in [-4121.5, +6579.6]. Below -4121.5 theta goes negative and the cosine folds
#' back: purity DECREASES as the tumour gets purer, and two different scores map
#' to one purity (S = -6000 and S = -2243 both give 0.962223 -- verified). Testing
#' the COSINE for [0,1] cannot catch this, because the folded values are in range.
#' The test is therefore on theta, not on p. Verified against the realised data:
#' observed ESTIMATEScore spans about -2689 to +4891, so no sample is in the fold
#' region -- this guard is latent here, and exists so it cannot become active
#' silently in a validation cohort.
estimate_purity_from_score <- function(estimate_score) {
  theta <- PURITY_A + PURITY_B * estimate_score
  p <- cos(theta)
  p[!is.na(theta) & (theta < 0 | theta > pi/2)] <- NA_real_
  p
}

#' ESTIMATE scores for one expression matrix.
#'
#' @param expr genes (rownames = HGNC symbols) x samples, on a LINEAR scale.
#'   ESTIMATE ranks within sample, so a monotone transform of the columns is
#'   immaterial, but the input is documented as linear TPM for clarity.
#' @return data.frame: sample, StromalScore, ImmuneScore, ESTIMATEScore, purity
#'
#' Pure: no globals beyond the registered constants, no file reads except the
#' package's own signature file, and the temporary files the package's file-based
#' API requires are created and removed inside this call.
estimate_scores <- function(expr, platform = "illumina", section = "B.purity") {
  if (!requireNamespace("estimate", quietly = TRUE))
    halt(section, "the estimate package is not available")
  si_geneset_path(section)
  if (is.null(rownames(expr))) halt(section, "expression matrix has no gene rownames")
  if (anyDuplicated(rownames(expr)))
    halt(section, "expression matrix has duplicate gene symbols; ESTIMATE ranks ",
         "by symbol and would use only the first")
  if (any(!is.finite(expr))) halt(section, "expression matrix has non-finite values")

  # filterCommonGenes() resolves `common_genes` from its own package environment,
  # which is not populated when the function is reached via `::` without the
  # package being attached. Attach it for the duration of this call rather than
  # globally, so the function stays self-contained.
  if (!"package:estimate" %in% search()) {
    suppressWarnings(suppressMessages(attachNamespace(asNamespace("estimate"))))
    on.exit(try(detach("package:estimate", unload = FALSE), silent = TRUE), add = TRUE)
  }
  utils::data("common_genes", package = "estimate", envir = environment())

  td <- tempfile("estimate_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  f_in <- file.path(td, "expr.txt"); f_gct <- file.path(td, "common.gct")
  f_out <- file.path(td, "scores.gct")

  df <- data.frame(GeneSymbol = rownames(expr), expr, check.names = FALSE,
                   stringsAsFactors = FALSE)
  write.table(df, f_in, sep = "\t", quote = FALSE, row.names = FALSE)
  # B8: capture the merge report rather than swallowing it -- the genes actually
  # scored are the intersection with ESTIMATE's 10,412-gene common set, not the
  # rows passed in, and the driver asserts that intersection is identical across
  # cohorts (it must be, or the scores are not comparable).
  utils::capture.output(
    estimate::filterCommonGenes(input.f = f_in, output.f = f_gct, id = "GeneSymbol"))
  # Count the merged GCT's rows rather than parsing the package's printed message:
  # the row count IS the number of genes scored, and a message-format change would
  # otherwise yield NA and make the cross-cohort comparability check pass vacuously.
  n_common <- nrow(utils::read.delim(f_gct, skip = 2, header = TRUE,
                                     check.names = FALSE, stringsAsFactors = FALSE))
  if (!is.finite(n_common) || n_common < 1L)
    halt(section, "could not determine how many genes ESTIMATE scored")
  utils::capture.output(
    estimate::estimateScore(input.ds = f_gct, output.ds = f_out, platform = platform))

  g <- utils::read.delim(f_out, skip = 2, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
  rownames(g) <- g$NAME
  vals <- t(as.matrix(g[, -(1:2), drop = FALSE]))
  need <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
  if (!all(need %in% colnames(vals)))
    halt(section, "ESTIMATE output lacks: ",
         paste(setdiff(need, colnames(vals)), collapse = ", "))
  # B7: sample identity travels through three file round trips and comes back
  # through read.table(header=TRUE), whose check.names mangles TCGA barcodes
  # (TCGA-W5-AA39-... -> TCGA.W5.AA39....). Assert the identity and ORDER rather
  # than trusting position: a reordering would attach every purity value to the
  # wrong patient with no error anywhere.
  if (!identical(rownames(vals), make.names(colnames(expr))))
    halt(section, "ESTIMATE returned samples in a different identity or order ",
         "than they were supplied; positional assignment would mis-attribute purity")

  out <- data.frame(
    sample        = colnames(expr),
    StromalScore  = as.numeric(vals[, "StromalScore"]),
    ImmuneScore   = as.numeric(vals[, "ImmuneScore"]),
    ESTIMATEScore = as.numeric(vals[, "ESTIMATEScore"]),
    stringsAsFactors = FALSE)
  # The package emits TumorPurity only for affymetrix; derive it explicitly so the
  # illumina path yields the same quantity from the same formula.
  out$purity_estimate <- if ("TumorPurity" %in% colnames(vals))
    as.numeric(vals[, "TumorPurity"]) else estimate_purity_from_score(out$ESTIMATEScore)
  if (nrow(out) != ncol(expr))
    halt(section, "ESTIMATE returned ", nrow(out), " rows for ", ncol(expr), " samples")
  # true count of genes ESTIMATE actually scored (B8), read from its own report
  attr(out, "n_genes_scored") <- if (length(n_common) && !is.na(n_common)) n_common else NA_integer_
  out
}

#' Read the Aran CPE table if it exists; NULL if it does not.
#'
#' Absence is a registered contingency (plan 3.4), not an error: the switch rule
#' resolves it. Returning NULL rather than halting is what lets the rule apply.
#' Aran et al. 2015 (Nat Commun 6:8971) Supplementary Data 1, "Tumor purity
#' estimates for TCGA samples". Downloaded 2026-08-02 from the publisher; md5
#' asserted so a different supplementary file cannot be substituted silently.
#'
#' Layout: three header rows before the column names, then
#' Sample ID | Cancer type | ESTIMATE | ABSOLUTE | LUMP | IHC | CPE.
#'
#' The table covers 21 TCGA cancer types from the 2015 freeze. FOUR of this
#' study's seven -- STAD, ESCA, PAAD, CHOL -- are not among them, so their CPE
#' coverage is genuinely 0% and the registered switch rule makes ESTIMATE primary
#' there. That is cohort coverage in the source, not a join failure: verified by
#' listing the table's own cancer types.
ARAN_MD5 <- "c459e6a965789b96860fc77bd346c681"
ARAN_ROWS <- 9364L

read_cpe <- function(path = CPE_PATH, section = "B.purity") {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("readxl", quietly = TRUE))
    halt(section, "readxl is needed to read the CPE workbook")
  md5 <- unname(tools::md5sum(path))
  if (!identical(md5, ARAN_MD5))
    halt(section, "aran_purity.xlsx md5 is ", md5, ", expected ", ARAN_MD5,
         ". This is not the registered Aran 2015 Supplementary Data 1.")
  d <- as.data.frame(suppressWarnings(suppressMessages(
    readxl::read_excel(path, sheet = 1, skip = 3))), stringsAsFactors = FALSE)
  if (ncol(d) < 7L) halt(section, "CPE workbook has ", ncol(d), " columns, expected >= 7")
  names(d)[1:7] <- c("sample", "cancer_type", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
  if (nrow(d) != ARAN_ROWS)
    halt(section, "CPE workbook has ", nrow(d), " rows, expected ", ARAN_ROWS)
  out <- data.frame(patient = substr(as.character(d$sample), 1, 12),
                    cancer_type = as.character(d$cancer_type),
                    CPE = suppressWarnings(as.numeric(d$CPE)),
                    stringsAsFactors = FALSE)
  # One CPE per patient: the table is sample-level and a patient may appear twice.
  out <- out[!is.na(out$CPE), , drop = FALSE]
  out[!duplicated(out$patient), , drop = FALSE]
}

#' Apply the registered switch rule for one cohort.
#'
#' CPE primary at coverage >= 80% of the cohort's eligible patients; ESTIMATE
#' primary below. Coverage is computed over the patients actually in the analysis
#' set, not over the CPE table, since that is the quantity the rule names.
select_purity <- function(df, threshold = CPE_COVERAGE_THRESHOLD) {
  n <- nrow(df)
  cov <- if (!"CPE" %in% names(df)) 0 else sum(!is.na(df$CPE)) / n
  primary <- if (cov >= threshold) "CPE" else "ESTIMATE"
  df$purity_source <- primary
  df$purity <- if (primary == "CPE") df$CPE else df$purity_estimate
  list(data = df, coverage = cov, source = primary)
}

#' Cross-method agreement, the registered sensitivity. NA when CPE is absent.
purity_agreement <- function(df) {
  if (!"CPE" %in% names(df)) return(c(n = 0, pearson = NA_real_, spearman = NA_real_))
  ok <- !is.na(df$CPE) & !is.na(df$purity_estimate)
  if (sum(ok) < 3L) return(c(n = sum(ok), pearson = NA_real_, spearman = NA_real_))
  c(n = sum(ok),
    pearson  = suppressWarnings(cor(df$CPE[ok], df$purity_estimate[ok], method = "pearson")),
    spearman = suppressWarnings(cor(df$CPE[ok], df$purity_estimate[ok], method = "spearman")))
}

#' Linear TPM matrix for the analysis samples of one cohort.
#'
#' Reads the assay -- which 05 was forbidden to do and 06 requires. Restricted to
#' the barcodes 05 retained, so purity is computed on exactly the analysis set.
cohort_expression <- function(se, barcodes, section = "B.purity") {
  if (!"tpm_unstrand" %in% SummarizedExperiment::assayNames(se))
    halt(section, "assay 'tpm_unstrand' is absent; available: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  keep <- match(barcodes, as.character(SummarizedExperiment::colData(se)$barcode))
  if (any(is.na(keep)))
    halt(section, "a barcode from the analysis set is absent from the cohort: ",
         paste(utils::head(barcodes[is.na(keep)], 3), collapse = ", "))
  x <- SummarizedExperiment::assay(se, "tpm_unstrand")[, keep, drop = FALSE]
  sym <- as.character(SummarizedExperiment::rowData(se)$gene_name)
  ok <- !is.na(sym) & nzchar(sym)
  x <- x[ok, , drop = FALSE]; sym <- sym[ok]
  # Duplicate symbols are collapsed by SUMMING TPM within symbol (audit finding
  # B5). An earlier revision kept the row with the largest mean in THIS cohort,
  # which made the retained transcript cohort-dependent: two cohorts could score
  # the same symbol from different Ensembl rows, and the purity covariate is
  # meta-analysed across cohorts. Summing is additively correct for TPM (the
  # units are per-million transcript fractions, so a gene's total is the sum over
  # its transcripts) and gives the identical answer in every cohort and in any
  # validation set.
  if (anyDuplicated(sym)) x <- rowsum(x, group = sym, reorder = TRUE)
  else rownames(x) <- sym
  colnames(x) <- barcodes
  x
}

# ==================================================================== driver
if (sys.nframe() == 0L) {

  OUTDIR <- "output"; dir.create(OUTDIR, showWarnings = FALSE)
  message("\n== B  purity (ESTIMATE + CPE switch rule). No score, no model. ==")

  use_estimate_lib()                       # B6: side effect confined to the driver
  clin <- read.csv(file.path(OUTDIR, "clinical_analysis_set.csv"), stringsAsFactors = FALSE)
  if (!nrow(clin)) halt("B.purity", "clinical analysis set is empty; run 05 first")

  si <- si_geneset_path()
  message("  ok  B.purity       SI_geneset.gmt verified (md5 ", SI_GENESET_MD5, ", 141+141 genes)")

  cpe <- read_cpe()
  if (is.null(cpe)) {
    # The plan (3.4) anticipates absence and the switch rule resolves it -- but a
    # missing DOWNLOAD is not the same fact as low coverage, and letting it fall
    # through the rule would invert B.j's registered CPE-primary ordering on an
    # acquisition gap rather than on the data (audit finding B2). Halt instead.
    halt("B.purity", "data/manual/aran_purity.xlsx is absent. B.j registers CPE as ",
         "the PRIMARY purity covariate; treating a missing file as 0% coverage ",
         "would silently invert that ordering for every cohort. Obtain Aran et al. ",
         "2015 Nat Commun 6:8971 Supplementary Data 1, or record an explicit ",
         "decision to proceed on ESTIMATE alone.")
  }
  message("  ok  B.purity       Aran CPE table loaded: ", nrow(cpe), " patients, ",
          length(unique(cpe$cancer_type)), " cancer types")
  absent <- setdiff(sub("^TCGA-", "", COHORTS), unique(cpe$cancer_type))
  if (length(absent))
    message("  ..  B.purity       NOT in the Aran 2015 freeze (coverage genuinely 0%): ",
            paste(absent, collapse = ", "),
            " -- the registered switch rule makes ESTIMATE primary for these.")

  res <- lapply(COHORTS, function(cc) {
    cl <- clin[clin$cohort == cc, , drop = FALSE]
    se <- readRDS(sprintf("data/tcga/%s_se.rds", cc))
    x  <- cohort_expression(se, cl$barcode)
    rm(se); gc(verbose = FALSE)
    sc <- estimate_scores(x)
    d  <- data.frame(cohort = cc, patient = cl$patient, barcode = cl$barcode,
                     sc[match(cl$barcode, sc$sample), -1, drop = FALSE],
                     stringsAsFactors = FALSE)
    if (!is.null(cpe)) d$CPE <- cpe$CPE[match(d$patient, cpe$patient)]
    s <- select_purity(d)
    a <- purity_agreement(s$data)
    list(data = s$data, coverage = s$coverage, source = s$source, agree = a,
         n_genes = attr(sc, "n_genes_scored"), n_genes_in = nrow(x))
  })
  names(res) <- COHORTS

  purity <- do.call(rbind, lapply(res, `[[`, "data"))
  summ <- do.call(rbind, lapply(COHORTS, function(cc) data.frame(
    cohort        = cc,
    n             = nrow(res[[cc]]$data),
    genes_scored  = res[[cc]]$n_genes,
    cpe_coverage  = round(res[[cc]]$coverage, 4),
    purity_source = res[[cc]]$source,
    # B4: the platform caveat belongs to the ESTIMATE-derived values only. CPE is
    # a consensus of four orthogonal methods and carries no Affymetrix conversion,
    # so labelling those rows "affymetrix_extrapolated" would be wrong.
    purity_calibration = if (res[[cc]]$source == "CPE") "aran_consensus_cpe"
                         else "estimate_affymetrix_extrapolated",
    n_purity_na   = sum(is.na(res[[cc]]$data$purity)),
    purity_mean   = round(mean(res[[cc]]$data$purity, na.rm = TRUE), 4),
    purity_sd     = round(sd(res[[cc]]$data$purity, na.rm = TRUE), 4),
    purity_min    = round(min(res[[cc]]$data$purity, na.rm = TRUE), 4),
    purity_max    = round(max(res[[cc]]$data$purity, na.rm = TRUE), 4),
    agree_n       = unname(res[[cc]]$agree["n"]),
    agree_pearson = round(unname(res[[cc]]$agree["pearson"]), 4),
    agree_spearman= round(unname(res[[cc]]$agree["spearman"]), 4),
    stringsAsFactors = FALSE)))

  # -- invariants -----------------------------------------------------------
  assert_n(nrow(purity), nrow(clin), "B.purity", "purity rows == clinical rows")
  if (anyDuplicated(purity$patient)) halt("B.purity", "duplicate patient in the purity table")
  ok <- !is.na(purity$purity)
  if (any(purity$purity[ok] < 0 | purity$purity[ok] > 1))
    halt("B.purity", "purity outside [0,1]: ",
         paste(round(range(purity$purity[ok]), 4), collapse = " to "))
  message("  ok  B.purity       purity bounded [0,1] for all ", sum(ok), " non-NA values")
  if (!all(summ$purity_source %in% c("CPE", "ESTIMATE")))
    halt("B.purity", "purity_source is neither CPE nor ESTIMATE")
  # B8: the genes ESTIMATE actually scored must be identical across cohorts, or
  # the scores are not on a common footing and the meta-analysis pools quantities
  # computed from different gene sets.
  if (anyNA(summ$genes_scored))
    halt("B.purity", "the scored-gene count is NA for ",
         paste(summ$cohort[is.na(summ$genes_scored)], collapse = ", "),
         "; the comparability check cannot pass on an unknown")
  if (length(unique(summ$genes_scored)) != 1L)
    halt("B.purity", "ESTIMATE scored a different number of genes per cohort (",
         paste(sprintf("%s=%d", sub("^TCGA-", "", summ$cohort), summ$genes_scored),
               collapse = ", "), "); the scores are not comparable")
  message("  ok  B.purity       ESTIMATE scored an identical gene set in all cohorts = ",
          summ$genes_scored[1])
  # B3: with CPE absent, the registered cross-method purity sensitivity cannot be
  # computed at all. Reported as a named deviation rather than emitted as a
  # silent column of NAs.
  if (all(is.na(summ$agree_pearson)))
    message("  !!  B.purity       REGISTERED SENSITIVITY NOT COMPUTABLE: the CPE-vs-ESTIMATE\n",
            "                     agreement analysis requires the Aran table, which is absent.\n",
            "                     agree_n = 0 in every cohort. Recorded in NOTES_FOR_REVIEW.")

  write.csv(purity, file.path(OUTDIR, "purity_per_patient.csv"), row.names = FALSE)
  write.csv(summ,   file.path(OUTDIR, "purity_summary.csv"), row.names = FALSE)

  message("\n-- purity per cohort --")
  print(summ, row.names = FALSE)
  message("\nB.purity complete. HARD STOP: no score constructed, no model fitted.")
}
