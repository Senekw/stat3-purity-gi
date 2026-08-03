#!/usr/bin/env Rscript
# 14_rppa_control.R
#
# ============================ EXPLORATORY, POST-HOC ==========================
# A POSITIVE CONTROL on the RPPA assay, not on the panel.
#
# 13_rppa_exploratory.R found the 140-gene score correlates with RPPA STAT3_pY705
# at a pooled r = 0.0964 (Wald 0.0041-0.1872, I2 59.8%) across the six discovery
# cohorts. That is weak. Two explanations are separable with data already in hand:
#   (a) the panel does not track STAT3 phosphorylation in these tissues, or
#   (b) the RPPA phospho measurement is itself noisy in these cohorts.
#
# The discriminating comparison is a pair of correlations whose expected value is
# NOT in question:
#   1. total STAT3 RPPA protein vs STAT3 mRNA -- if the assay and the transcript
#      agree for the SAME gene, the platform is working.
#   2. the 140-gene score vs STAT3 mRNA -- the panel's internal coherence, with
#      the committed FU-iCCA value (0.5409) as the comparator.
#   3. total STAT3 RPPA vs STAT3_pY705 RPPA -- within-platform, same lysate.
#
# NOTHING here is registered. No amendment. No survival model. No registered
# quantity is altered, and no output file of 13 is overwritten.
#
# The two antibodies are DISTINCT and must not be conflated:
#   STAT3_pY705  AGID00388  catalog 9131  -- the phosphosite
#   Stat3        AGID00185  catalog 4904  -- TOTAL protein
# =============================================================================

suppressPackageStartupMessages({ library(stats); library(metafor) })

OUTDIR  <- "output"
SECT    <- "EXPLORATORY"
TAG     <- "EXPLORATORY_POSTHOC"
BOTH    <- "data/validation/RPPA/tcga_rppa_stat3_both.csv"
BOTH_MD5 <- "1bc5d6cedb9ad2e2fc15ef6af71c90e2"
COHORTS <- c("TCGA-COAD","TCGA-READ","TCGA-STAD","TCGA-ESCA","TCGA-PAAD","TCGA-LIHC")
MIN_N   <- 50L
FUICCA_SCORE_VS_MRNA <- 0.5409     # committed comparator, output/fuicca_correlations.csv

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}

#' Bind 07's committed scoring path and prove it is the committed one -- bodies
#' against a FRESH sourcing, with a negative control, as 09/10/12/13 do.
bind_committed <- function() {
  e07 <- new.env(); sys.source("07_score.R",    envir = e07)
  e08 <- new.env(); sys.source("08_survival.R", envir = e08)
  fns <- list(score_cohort = e07$score_cohort, zero_variance_genes = e07$zero_variance_genes,
              digest_genes = e07$digest_genes, read_gene_list = e07$read_gene_list,
              expression_log2tpm = e07$expression_log2tpm,
              GENE_LIST_PRIMARY = e07$GENE_LIST_PRIMARY, meta_one = e08$meta_one)
  if (any(vapply(fns, is.null, logical(1)))) halt(SECT, "a function is missing from 07/08")
  c7 <- new.env(); sys.source("07_score.R",    envir = c7)
  c8 <- new.env(); sys.source("08_survival.R", envir = c8)
  same <- function(a, b) identical(deparse(body(a)), deparse(body(b)))
  for (nm in c("score_cohort","zero_variance_genes","digest_genes","read_gene_list",
               "expression_log2tpm"))
    if (!same(fns[[nm]], c7[[nm]])) halt(SECT, nm, " != 07_score.R")
  if (!same(fns$meta_one, c8$meta_one)) halt(SECT, "meta_one != 08_survival.R")
  if (same(fns$score_cohort, fns$meta_one))
    halt(SECT, "the body comparison cannot discriminate two different functions")
  message("  ok  ", SECT, "  6 functions match their committed definitions (negative control passed)")
  fns
}

#' Fisher z for Pearson; Bonett-Wright SE for Spearman. Same as 12 and 13.
cor_report <- function(x, y, method, cohort, comparison) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 4L) return(data.frame(cohort = cohort, comparison = comparison, method = method,
                                n = n, r = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                                p = NA_real_, label = TAG, stringsAsFactors = FALSE))
  r <- suppressWarnings(cor(x[ok], y[ok], method = method))
  se <- if (method == "pearson") 1/sqrt(n-3) else sqrt((1 + r^2/2)/(n-3))
  ci <- tanh(atanh(r) + c(-1,1) * qnorm(0.975) * se)
  tt <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.frame(cohort = cohort, comparison = comparison, method = method, n = n, r = r,
             ci_lo = ci[1], ci_hi = ci[2], p = unname(tt$p.value), label = TAG,
             stringsAsFactors = FALSE)
}

#' Pool Pearson correlations on the Fisher z scale, reporting heterogeneity.
pool_z <- function(rows, meta_one, comparison) {
  pe <- rows[rows$method == "pearson" & is.finite(rows$r), ]
  if (nrow(pe) < 2L) return(NULL)
  mz <- meta_one(atanh(pe$r), 1/(pe$n - 3), pe$cohort)
  data.frame(comparison = comparison, k = mz$k, n_total = sum(pe$n),
             r_pooled = tanh(mz$est),
             ci_lo_wald = tanh(mz$ci_lo_wald), ci_hi_wald = tanh(mz$ci_hi_wald),
             ci_lo_hksj = tanh(mz$ci_lo_hksj), ci_hi_hksj = tanh(mz$ci_hi_hksj),
             p_wald = mz$p_wald, tau2 = mz$tau2, I2 = mz$I2, Q = mz$Q, Q_df = mz$Q_df,
             Q_p = mz$Q_p, pi_lo = tanh(mz$pi_lo), pi_hi = tanh(mz$pi_hi),
             label = TAG, stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------- driver
if (sys.nframe() == 0L) {
  message("\n== EXPLORATORY, POST-HOC: RPPA positive control ==")
  dir.create(OUTDIR, showWarnings = FALSE)
  if (!file.exists("14_rppa_control.R")) halt(SECT, "run from the repository root")
  fns <- bind_committed()

  if (!file.exists(BOTH)) halt(SECT, "RPPA two-antibody extraction absent: ", BOTH)
  got <- unname(tools::md5sum(BOTH))
  if (!identical(got, BOTH_MD5))
    halt(SECT, "RPPA extraction md5 is ", got, ", expected ", BOTH_MD5)

  rp <- read.csv(BOTH, stringsAsFactors = FALSE, colClasses = "character")
  for (cc in c("pY705","total_STAT3")) rp[[cc]] <- suppressWarnings(as.numeric(rp[[cc]]))
  if (nrow(rp) != 1282L) halt(SECT, "expected 1282 RPPA rows, got ", nrow(rp))

  # THE SAME sample-type and aliquot rules as 13, restated so they cannot drift.
  # Sample type is re-derived from the barcode: read.csv coerces an extracted
  # "01" to integer 1, which silently emptied the analysis set in 13 before the
  # audit caught it.
  rp$stype <- substr(vapply(strsplit(rp$aliquot, "-"), function(x)
                            if (length(x) >= 4) x[4] else "", character(1)), 1, 2)
  if (!all(nzchar(rp$stype))) halt(SECT, "unparseable sample-type code")
  # F6 (audit): 13's filter is `stype == "01" & is.finite(pY705)`. Matched
  # verbatim rather than relying on the is.finite clause being a no-op here (all
  # 1,282 extracted pY705 values are finite, so both select the same 1,278 rows).
  rp <- rp[rp$stype == "01" & is.finite(rp$pY705), , drop = FALSE]
  if (anyDuplicated(rp$patient))
    halt(SECT, "more than one primary-tumour aliquot for: ",
         paste(unique(rp$patient[duplicated(rp$patient)]), collapse = ", "),
         "; a tie-break must be registered before deduplication")
  rp <- rp[order(rp$patient, rp$aliquot), ]     # deterministic, as in 13

  # the pY705 column must be the SAME numbers 13 used
  old <- read.csv("data/validation/RPPA/tcga_rppa_stat3_py705.csv",
                  stringsAsFactors = FALSE, colClasses = "character")
  old$pY705 <- suppressWarnings(as.numeric(old$pY705))
  j <- match(rp$aliquot, old$aliquot)
  if (anyNA(j) || !isTRUE(all.equal(rp$pY705, old$pY705[j], tolerance = 1e-12)))
    halt(SECT, "the pY705 column does not reproduce the extraction 13 used")
  message(sprintf("  ok  %s  %d primary-tumour aliquots | pY705 reproduces 13's extraction exactly",
                  SECT, nrow(rp)))
  message(sprintf("  ..  %s  total STAT3 present for %d of %d", SECT,
                  sum(is.finite(rp$total_STAT3)), nrow(rp)))

  clin <- read.csv(file.path(OUTDIR, "clinical_analysis_set.csv"), stringsAsFactors = FALSE)
  scr  <- read.csv(file.path(OUTDIR, "scores_per_patient.csv"), stringsAsFactors = FALSE)
  g140 <- fns$read_gene_list(fns$GENE_LIST_PRIMARY, SECT)

  all_rows <- list(); pp <- list()
  for (cc in COHORTS) {
    cl <- clin[clin$cohort == cc, ]
    se_obj <- readRDS(file.path("data/tcga", paste0(cc, "_se.rds")))
    expr <- fns$expression_log2tpm(se_obj, cl$barcode, SECT)
    if (!"STAT3" %in% rownames(expr)) halt(SECT, cc, ": STAT3 absent from the expression matrix")
    if (sum(rownames(expr) == "STAT3") != 1L) halt(SECT, cc, ": STAT3 is not a unique row")
    zv <- fns$zero_variance_genes(expr, g140)
    if (length(zv)) halt(SECT, cc, ": zero-variance scoring gene(s)")
    s <- fns$score_cohort(expr, g140, SECT)
    stat3_mrna <- setNames(expr["STAT3", ], colnames(expr))
    names(s) <- cl$patient[match(names(s), cl$barcode)]
    names(stat3_mrna) <- cl$patient[match(names(stat3_mrna), cl$barcode)]
    if (anyNA(names(s)) || anyNA(names(stat3_mrna)))
      halt(SECT, cc, ": a scored barcode has no patient id")

    # the recomputation must reproduce the committed registered scores
    ref <- scr[scr$cohort == cc, c("patient","score")]
    kk <- intersect(names(s), ref$patient)
    if (!isTRUE(all.equal(unname(s[kk]), ref$score[match(kk, ref$patient)], tolerance = 1e-10)))
      halt(SECT, cc, ": the recomputed score does not reproduce scores_per_patient.csv")

    r <- rp[rp$cohort == cc, ]
    yP <- setNames(r$pY705, r$patient); yT <- setNames(r$total_STAT3, r$patient)
    k  <- intersect(names(s), names(yP))
    if (length(k) < MIN_N) halt(SECT, cc, ": joined n = ", length(k), " below ", MIN_N)

    # COMPARISON 2 IS PARTLY A SELF-CORRELATION: STAT3 is itself one of the 140
    # scoring genes (1/140 = 0.71% of the mean). The committed FU-iCCA comparator
    # 0.5409 has the identical property, so the two are like-for-like -- but the
    # leave-STAT3-out variant is computed alongside so the size of the artefact
    # is visible rather than argued.
    s139 <- fns$score_cohort(expr, setdiff(g140, "STAT3"), SECT)
    names(s139) <- cl$patient[match(names(s139), cl$barcode)]

    for (mth in c("pearson","spearman")) {
      all_rows[[length(all_rows)+1]] <- rbind(
        cor_report(yT[k], stat3_mrna[k], mth, cc, "1. total STAT3 RPPA vs STAT3 mRNA"),
        cor_report(s[k],  stat3_mrna[k], mth, cc, "2. score_140 vs STAT3 mRNA"),
        cor_report(s139[k], stat3_mrna[k], mth, cc,
                   "2b. score_139 (STAT3 dropped) vs STAT3 mRNA"),
        cor_report(yT[k], yP[k],         mth, cc, "3. total STAT3 RPPA vs STAT3_pY705 RPPA"),
        cor_report(s[k],  yP[k],         mth, cc, "4. score_140 vs STAT3_pY705 RPPA (from 13)"))
    }
    pp[[cc]] <- data.frame(cohort = cc, patient = k, score_140 = unname(s[k]),
                           stat3_mrna = unname(stat3_mrna[k]),
                           rppa_total_STAT3 = unname(yT[k]), rppa_pY705 = unname(yP[k]),
                           label = TAG, stringsAsFactors = FALSE)
    message(sprintf("  ..  %s  %-10s n = %3d | scored, reproduces committed values", SECT, cc, length(k)))
  }
  cr <- do.call(rbind, all_rows); rownames(cr) <- NULL
  cr <- cr[order(cr$comparison, cr$method, cr$cohort), ]

  pool <- do.call(rbind, lapply(unique(cr$comparison),
                                function(x) pool_z(cr[cr$comparison == x, ], fns$meta_one, x)))
  # F7 (audit): pool_z returns NULL below two poolable cohorts; rbind of all-NULL
  # is NULL, and the next line would fail with "object of type NULL is not
  # subsettable" -- an uninformative R error where this project's rule is that a
  # fired assertion is a finding.
  if (is.null(pool) || !nrow(pool))
    halt(SECT, "no comparison had two or more poolable cohorts")
  # attach the committed FU-iCCA comparator to comparison 2 only (not 2b, which
  # has no FU-iCCA counterpart)
  pool$fuicca_comparator <- ifelse(pool$comparison == "2. score_140 vs STAT3 mRNA",
                                   FUICCA_SCORE_VS_MRNA, NA_real_)
  if (sum(!is.na(pool$fuicca_comparator)) != 1L)
    halt(SECT, "the FU-iCCA comparator attached to ",
         sum(!is.na(pool$fuicca_comparator)), " rows, expected exactly 1")
  # WHAT EACH COMPARISON CAN AND CANNOT ESTABLISH. Without this the script emits
  # four numbers and no rule for reading them, which invites the reader to treat
  # whichever is largest as validation.
  interp <- c(
    "1." = paste0("Protein-mRNA agreement for ONE gene across two independent assays. ",
                  "Typical TCGA protein-mRNA r is modest (~0.4-0.5 for well-behaved ",
                  "antibodies), so a value in that range is consistent with a working ",
                  "platform and a value near zero would indict the total-STAT3 antibody ",
                  "in that cohort. NO threshold is prespecified; this is descriptive."),
    "2." = paste0("Panel internal coherence. PARTLY A SELF-CORRELATION: STAT3 is 1 of the ",
                  "140 scoring genes (0.71% of the mean). Comparator: FU-iCCA 0.5409, ",
                  "which has the same property. See 2b for the leave-STAT3-out value."),
    "2b." = "As 2, with STAT3 excluded from the score. The difference from 2 is the self-correlation artefact.",
    "3." = paste0("WITHIN-PLATFORM: both antibodies are measured on the SAME lysate, the ",
                  "same array, and are normalised together. It is NOT independent evidence ",
                  "that either antibody is valid -- shared loading and normalisation inflate ",
                  "it. A high value is uninformative about assay validity; a LOW value would ",
                  "be surprising and would suggest one antibody is noise."),
    "4." = "Reproduces 13's committed score-vs-pY705 correlations; the quantity under test.")
  pool$can_establish <- unname(interp[sub("^([0-9]+b?\\.).*$", "\\1", pool$comparison)])
  pool$note <- paste0("EXPLORATORY, POST-HOC. Back-transformed pooled Fisher z (random ",
                      "effects, REML); Wald interval quoted, as in 13. tanh of a pooled z is ",
                      "not an unbiased mean of the per-cohort r. The four comparisons are NOT ",
                      "on equal evidential footing -- see can_establish.")

  # verify comparison 4 reproduces 13 exactly -- same patients, same numbers
  prev <- read.csv(file.path(OUTDIR, "exploratory_rppa_correlations.csv"), stringsAsFactors = FALSE)
  chk <- merge(cr[grepl("^4\\.", cr$comparison), c("cohort","method","n","r")],
               prev[prev$cohort != "POOLED", c("cohort","method","n","r")],
               by = c("cohort","method"), suffixes = c("_new","_13"))
  # F5 (audit): two load-bearing properties, stated so a future edit cannot
  # silently defeat them. (a) `nrow(chk) != 12L` is what forecloses a VACUOUS
  # pass -- a 0-row merge fails here rather than sailing through an empty any().
  # (b) prev$r is stored at exactly 4 dp (verified), so round(r_new, 4) is the
  # right comparison and 1e-9 is the right tolerance for it.
  if (nrow(chk) != 12L || any(chk$n_new != chk$n_13) ||
      any(abs(round(chk$r_new, 4) - chk$r_13) > 1e-9))
    halt(SECT, "comparison 4 does not reproduce 13's committed score-vs-pY705 values")
  message(sprintf("  ok  %s  comparison 4 reproduces all %d of 13's committed correlations",
                  SECT, nrow(chk)))

  cr$r <- round(cr$r, 4); cr$ci_lo <- round(cr$ci_lo, 4); cr$ci_hi <- round(cr$ci_hi, 4)
  write.csv(cr,   file.path(OUTDIR, "exploratory_rppa_control_correlations.csv"), row.names = FALSE)
  write.csv(pool, file.path(OUTDIR, "exploratory_rppa_control_pooled.csv"), row.names = FALSE)
  write.csv(do.call(rbind, pp),
            file.path(OUTDIR, "exploratory_rppa_control_per_patient.csv"), row.names = FALSE)

  writeLines(c(
    "# EXPLORATORY, POST-HOC positive control -- not registered, not an amendment",
    paste0("run_date: ", Sys.Date()),
    "Purpose: separate panel validity from RPPA assay performance, using the same",
    "1,282 files already downloaded. No re-fetch.",
    "  STAT3_pY705  AGID00388  catalog 9131  -- phosphosite",
    "  Stat3        AGID00185  catalog 4904  -- TOTAL protein, a DIFFERENT antibody",
    paste0("  extraction md5: ", BOTH_MD5),
    "Aliquot rule and sample-type derivation identical to 13_rppa_exploratory.R.",
    "The recomputed 140-gene score reproduces scores_per_patient.csv in every cohort,",
    "and comparison 4 reproduces 13's committed score-vs-pY705 correlations exactly.",
    paste0("script_md5: ", unname(tools::md5sum("14_rppa_control.R"))),
    "COMPARISONS 1 and 3 rest on fewer patients than 2, 2b and 4: total STAT3 is",
    "absent from 19 of the 1,282 files, so their n is smaller. Each row reports its own n.",
    "COMPARISON 2 IS PARTLY A SELF-CORRELATION: STAT3 is 1 of the 140 scoring genes",
    "(0.71% of the mean). The FU-iCCA comparator 0.5409 has the same property.",
    "Comparison 2b repeats it with STAT3 dropped from the score.",
    "COMPARISON 3 IS WITHIN-PLATFORM: both antibodies come from the same lysate and",
    "array and are normalised together, so it is not independent evidence of validity.",
    "NO survival model. NO amendment. No registered quantity altered."),
    file.path(OUTDIR, "exploratory_rppa_control_provenance.txt"))

  message("\n-- per cohort (EXPLORATORY, POST-HOC) --")
  print(cr[cr$method == "pearson", c("comparison","cohort","n","r","ci_lo","ci_hi","p")],
        row.names = FALSE)
  message("\n-- pooled (EXPLORATORY, POST-HOC) --")
  print(pool[, c("comparison","k","n_total","r_pooled","ci_lo_wald","ci_hi_wald","I2","Q_p")],
        row.names = FALSE)
  message("\nEXPLORATORY control complete. HARD STOP: no survival model, no amendment.")
}
