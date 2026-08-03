#!/usr/bin/env Rscript
# 13_rppa_exploratory.R
#
# ============================ EXPLORATORY, POST-HOC ==========================
# NEITHER analysis in this script is registered. Both were requested after the
# Part B primary result, the null benchmark and the FU-iCCA concordance were all
# known. Every output file and every printed table carries the label.
#
#   1. Colorectal pooling: inverse-variance meta-analysis of the M1 score log-HR
#      across TCGA-COAD, TCGA-READ and GSE39582. This is NOT the registered
#      discovery meta-analysis: GSE39582 is an EXTERNAL validation cohort, and
#      Amendment 16 states validation cohorts "are not meta-analysed with
#      discovery". Pooling them here is a post-hoc precision exercise, reported
#      as such, and it does not replace or amend the registered pooled estimate.
#
#   2. TCGA RPPA STAT3_pY705 concordance: the 140-gene score against directly
#      measured pY705 within each discovery cohort, closing the transportability
#      gap left by FU-iCCA (intrahepatic cholangiocarcinoma, not one of the six).
#
# PRIOR-KNOWLEDGE DISCLOSURE, required by the run order: the author's prior ESMO
# Asia work used TCGA RPPA STAT3_pY705 in TCGA-CHOL. This data source is
# therefore NOT new to the author. TCGA-CHOL is excluded from analysis 2 on that
# ground and because it is not a meta-eligible discovery cohort (Amendment 8).
#
# NO survival model is fitted here. NO amendment is written. Nothing in this
# script alters a registered quantity.
# =============================================================================

suppressPackageStartupMessages({ library(stats); library(metafor) })

OUTDIR <- "output"
SECT   <- "EXPLORATORY"
TAG    <- "EXPLORATORY_POSTHOC"
RPPA   <- "data/validation/RPPA/tcga_rppa_stat3_py705.csv"
COHORTS <- c("TCGA-COAD","TCGA-READ","TCGA-STAD","TCGA-ESCA","TCGA-PAAD","TCGA-LIHC")
MIN_N  <- 50L      # run order: stop if coverage is below 50 in EVERY cohort

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}

#' Bind the committed functions and prove they are the committed ones -- the
#' same guard 09/10/12 use: bodies compared against a FRESH sourcing, with a
#' negative control, because comparing a binding to its own source cannot fail.
bind_committed <- function() {
  e07 <- new.env(); sys.source("07_score.R",    envir = e07)
  e08 <- new.env(); sys.source("08_survival.R", envir = e08)
  fns <- list(score_cohort = e07$score_cohort, zero_variance_genes = e07$zero_variance_genes,
              digest_genes = e07$digest_genes, read_gene_list = e07$read_gene_list,
              expression_log2tpm = e07$expression_log2tpm,
              GENE_LIST_PRIMARY = e07$GENE_LIST_PRIMARY, GENE_LIST_SENS = e07$GENE_LIST_SENS,
              meta_one = e08$meta_one)
  miss <- names(fns)[vapply(fns, is.null, logical(1))]
  if (length(miss)) halt(SECT, "not found in 07/08: ", paste(miss, collapse = ", "))
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

#' Fisher z for Pearson; Bonett-Wright SE for Spearman. Same function as 12.
cor_report <- function(x, y, method, label, n_lab = "") {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 4L) return(data.frame(comparison = label, method = method, n = n, r = NA_real_,
                                ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_,
                                label = TAG, stringsAsFactors = FALSE))
  r <- suppressWarnings(cor(x[ok], y[ok], method = method))
  se <- if (method == "pearson") 1/sqrt(n-3) else sqrt((1 + r^2/2)/(n-3))
  ci <- tanh(atanh(r) + c(-1,1) * qnorm(0.975) * se)
  tt <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.frame(comparison = label, method = method, n = n, r = r,
             ci_lo = ci[1], ci_hi = ci[2], p = unname(tt$p.value),
             label = TAG, stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------- driver
if (sys.nframe() == 0L) {
  message("\n== EXPLORATORY, POST-HOC: colorectal pooling + TCGA RPPA pY705 concordance ==")
  dir.create(OUTDIR, showWarnings = FALSE)
  fns <- bind_committed()

  # ================= 1. colorectal pooling (EXPLORATORY) ===================
  pc <- read.csv(file.path(OUTDIR, "survival_per_cohort.csv"), stringsAsFactors = FALSE)
  vg <- read.csv(file.path(OUTDIR, "validation_gse39582_models.csv"), stringsAsFactors = FALSE)
  crc <- rbind(pc[pc$cohort %in% c("TCGA-COAD","TCGA-READ") & pc$model == "M1",
                  c("cohort","n","events","beta","se")],
               vg[vg$model == "M1", c("cohort","n","events","beta","se")])
  if (nrow(crc) != 3L) halt(SECT, "expected 3 colorectal cohorts, got ", nrow(crc))
  if (any(!is.finite(crc$beta) | !is.finite(crc$se) | crc$se <= 0))
    halt(SECT, "a colorectal input is missing or has non-positive SE")
  crc$source <- ifelse(crc$cohort == "GSE39582", "external validation", "discovery")

  # meta_one takes vi = VARIANCE (passed straight to metafor::rma). My first
  # draft passed `se`, and my comment then claimed that would give a spuriously
  # TIGHT interval. That is backwards, and the audit caught it: every |se| here
  # is below 1, so se > se^2 -- treating SEs as variances INFLATES them,
  # understates the weights, and widens the interval (Wald width 0.718 against
  # the correct 0.218). Verified by running both.
  m <- fns$meta_one(crc$beta, crc$se^2, crc$cohort)
  if (is.null(m)) halt(SECT, "the colorectal pooling returned NULL")

  crc_out <- data.frame(
    analysis = "M1_colorectal_pooled", label = TAG,
    k = m$k, cohorts = paste(m$cohorts, collapse = "+"),
    n_total = sum(crc$n), events_total = sum(crc$events),
    est = m$est, se_hksj = m$se_hksj,
    ci_lo_hksj = m$ci_lo_hksj, ci_hi_hksj = m$ci_hi_hksj, p_hksj = m$p_hksj,
    ci_lo_wald = m$ci_lo_wald, ci_hi_wald = m$ci_hi_wald, p_wald = m$p_wald,
    fe_est = m$fe_est, fe_ci_lo = m$fe_ci_lo, fe_ci_hi = m$fe_ci_hi,
    HR = exp(m$est), HR_lo_hksj = exp(m$ci_lo_hksj), HR_hi_hksj = exp(m$ci_hi_hksj),
    HR_lo_wald = exp(m$ci_lo_wald), HR_hi_wald = exp(m$ci_hi_wald),
    tau2 = m$tau2, I2 = m$I2, Q = m$Q, Q_df = m$Q_df, Q_p = m$Q_p,
    pi_lo = m$pi_lo, pi_hi = m$pi_hi,
    HR_pi_lo = exp(m$pi_lo), HR_pi_hi = exp(m$pi_hi),
    # F3 (audit): renamed from excludes_HR_above_*. A 95% upper bound does not
    # "exclude" everything above it; it is the largest per-SD HR compatible with
    # these data at two-sided alpha = 0.05 under this model.
    compat_upper_HR_per_SD_wald = exp(m$ci_hi_wald),
    compat_upper_HR_per_SD_hksj = exp(m$ci_hi_hksj),
    note = paste0("EXPLORATORY, POST-HOC. NOT the registered discovery meta-analysis: ",
                  "GSE39582 is an external validation cohort and Amendment 16 states ",
                  "validation cohorts are not meta-analysed with discovery. With k = 3 and ",
                  "Q p = ", signif(m$Q_p, 3), ", the Hartung-Knapp interval is NARROWER than ",
                  "the Wald interval -- the known small-k behaviour when heterogeneity is ",
                  "near zero. The WALD/FE interval is the conservative one and is the bound ",
                  "quoted. Registered discovery pooled M1 for comparison: 0.121174 ",
                  "(HR 1.1288, HKSJ 1.0002-1.2740) over SIX cohorts; this row is two of those ",
                  "six plus one external cohort and is NOT a re-estimate of it. The three ",
                  "inputs are not exchangeable: COAD and READ are RNA-seq scored by ",
                  "expression_log2tpm, GSE39582 is GPL570 microarray scored through a ",
                  "probe-collapse path with a KRT17 rescue; and READ's registered endpoint is ",
                  "PFI (Amendment 7) while COAD and GSE39582 are OS. Pooling across differing ",
                  "endpoints and platforms is the substantive limitation, not the arithmetic."),
    stringsAsFactors = FALSE)
  write.csv(crc_out, file.path(OUTDIR, "exploratory_colorectal_pooled.csv"), row.names = FALSE)
  write.csv(crc, file.path(OUTDIR, "exploratory_colorectal_inputs.csv"), row.names = FALSE)
  message(sprintf("  ok  %s  colorectal pooled: HR %.4f (Wald %.4f-%.4f), I2 %.2f%%, k=%d, n=%d, events=%d",
                  SECT, exp(m$est), exp(m$ci_lo_wald), exp(m$ci_hi_wald), m$I2,
                  m$k, sum(crc$n), sum(crc$events)))

  # ================= 2. RPPA pY705 concordance (EXPLORATORY) ===============
  if (!file.exists(RPPA)) halt(SECT, "RPPA extraction absent: ", RPPA)
  # F5 (audit): this is the only input in the project not bound to a digest.
  # Pin it, so the AGID00388 provenance claim refers to a specific file.
  RPPA_MD5 <- "71c6bbff3ab4207cfd520c1e883dfa19"
  got_md5 <- unname(tools::md5sum(RPPA))
  if (!identical(got_md5, RPPA_MD5))
    halt(SECT, "RPPA extraction md5 is ", got_md5, ", expected ", RPPA_MD5)
  rp <- read.csv(RPPA, stringsAsFactors = FALSE, colClasses = "character")
  rp$pY705 <- suppressWarnings(as.numeric(rp$pY705))
  if (!identical(names(rp)[1:6], c("cohort","file_id","aliquot","patient",
                                   "sample_type_code","pY705")))
    halt(SECT, "RPPA extraction columns changed: ", paste(names(rp), collapse = ", "))
  if (nrow(rp) != 1282L) halt(SECT, "RPPA extraction has ", nrow(rp),
                              " rows, expected 1282")
  # SAMPLE TYPE is re-derived from the barcode, not read from the CSV column.
  # read.csv coerces the extracted "01" to INTEGER 1, so a `== "01"` test matched
  # NOTHING and silently produced an empty analysis set -- caught before the run.
  # colClasses="character" fixes the read; deriving from the barcode as well means
  # the filter cannot depend on how the intermediate was serialised.
  rp$stype <- substr(vapply(strsplit(rp$aliquot, "-"), function(x)
                            if (length(x) >= 4) x[4] else "", character(1)), 1, 2)
  if (!all(nzchar(rp$stype))) halt(SECT, sum(!nzchar(rp$stype)),
                                   " RPPA aliquot(s) have no parseable sample-type code")
  n_all <- nrow(rp)
  rp <- rp[rp$stype == "01" & is.finite(rp$pY705), , drop = FALSE]
  if (!nrow(rp)) halt(SECT, "no primary-tumour RPPA rows survived the filter")
  # F6 (audit): the dedup below has no registered tie-break, so assert that it
  # never has to choose. On this download no patient has two primary-tumour
  # aliquots; if a re-download introduced one, a rule must be registered first.
  if (anyDuplicated(rp$patient))
    halt(SECT, "more than one primary-tumour RPPA aliquot for: ",
         paste(unique(rp$patient[duplicated(rp$patient)]), collapse = ", "),
         "; a tie-break must be registered before deduplication")
  message(sprintf("  ..  %s  RPPA: %d of %d rows are primary tumour (code 01) with a finite pY705",
                  SECT, nrow(rp), n_all))
  # one aliquot per patient: the analysis set's own barcode where possible
  clin <- read.csv(file.path(OUTDIR, "clinical_analysis_set.csv"), stringsAsFactors = FALSE)
  # F7 (audit): this was read and never used -- the trace of a cross-check that
  # was not performed. The loop RE-COMPUTES the score from the SE objects, so it
  # must reproduce the committed registered values exactly; asserted below.
  scr  <- read.csv(file.path(OUTDIR, "scores_per_patient.csv"), stringsAsFactors = FALSE)

  # F9 (audit): CHOL appears in the COVERAGE table for completeness and is
  # excluded from analysis. Both stated grounds are recorded: it is not
  # meta-eligible (Amendment 8, all 35 rows meta_eligible = FALSE), and the
  # author's prior ESMO Asia work used TCGA RPPA pY705 in exactly this cohort.
  cov <- do.call(rbind, lapply(c(COHORTS, "TCGA-CHOL"), function(cc) {
    # ONE ALIQUOT PER PATIENT, deterministically. `!duplicated()` alone keeps
    # whichever row the tar happened to yield first, which is the ordering
    # dependence the clinical script had to register a tie-break for. Sorting the
    # full barcode first makes the pick reproducible regardless of extraction
    # order. (On this download no patient has two primary-tumour aliquots, so the
    # rule is latent -- it exists so a re-download cannot change the result.)
    r <- rp[rp$cohort == cc, ]
    r <- r[order(r$patient, r$aliquot), ]
    n_dup <- sum(duplicated(r$patient))
    r <- r[!duplicated(r$patient), ]
    cl <- clin[clin$cohort == cc, ]
    data.frame(cohort = cc, label = TAG,
               rppa_files = sum(rp$cohort == cc), rppa_patients = nrow(r),
               dup_aliquots_dropped = n_dup,
               analysis_set_patients = nrow(cl),
               overlap = length(intersect(r$patient, cl$patient)),
               analysed = cc %in% COHORTS,
               note = if (cc == "TCGA-CHOL")
                 paste0("descriptive only: not meta-eligible (Amendment 8) AND the ",
                        "author's prior ESMO Asia work used TCGA RPPA STAT3_pY705 in ",
                        "this cohort -- prior knowledge, disclosed") else "",
               stringsAsFactors = FALSE)
  }))
  write.csv(cov, file.path(OUTDIR, "exploratory_rppa_coverage.csv"), row.names = FALSE)
  message("\n  RPPA STAT3_pY705 coverage:")
  print(cov[, c("cohort","rppa_patients","analysis_set_patients","overlap","analysed")],
        row.names = FALSE)

  # run order: REPORT AND STOP if coverage is below 50 in EVERY cohort
  if (all(cov$overlap[cov$analysed] < MIN_N))
    halt(SECT, "RPPA pY705 overlap is below ", MIN_N, " in every cohort (max ",
         max(cov$overlap), "); the run order stops here.")
  qual <- cov$cohort[cov$analysed & cov$overlap >= MIN_N]
  message(sprintf("  ok  %s  %d cohort(s) at or above n = %d: %s",
                  SECT, length(qual), MIN_N, paste(sub("TCGA-","",qual), collapse = ", ")))

  # score by the committed path, per qualifying cohort
  g140 <- fns$read_gene_list(fns$GENE_LIST_PRIMARY, SECT)
  res <- list(); per <- list()
  for (cc in qual) {
    # Signature read from 07, not recalled: expression_log2tpm(se, barcodes,
    # section). It selects the analysis-set barcodes itself and halts if any is
    # absent, so no post-hoc subsetting is needed -- and passing only `se` would
    # have errored on a missing `barcodes`.
    cl <- clin[clin$cohort == cc, ]
    se_obj <- readRDS(file.path("data/tcga", paste0(cc, "_se.rds")))
    expr <- fns$expression_log2tpm(se_obj, cl$barcode, SECT)
    zv <- fns$zero_variance_genes(expr, g140)
    if (length(zv)) halt(SECT, cc, ": zero-variance scoring gene(s): ",
                         paste(zv, collapse = ", "))
    s <- fns$score_cohort(expr, g140, SECT)
    names(s) <- cl$patient[match(names(s), cl$barcode)]
    if (anyNA(names(s))) halt(SECT, cc, ": a scored barcode has no patient id")
    # F7: the recomputation must reproduce the committed registered scores
    # exactly -- same reader, same gene list, same patients. If it does not, the
    # concordance would be computed on a score that is not the study's.
    ref <- scr[scr$cohort == cc, c("patient", "score")]
    kk <- intersect(names(s), ref$patient)
    if (length(kk) < 1L) halt(SECT, cc, ": no overlap with scores_per_patient.csv")
    if (!isTRUE(all.equal(unname(s[kk]), ref$score[match(kk, ref$patient)],
                          tolerance = 1e-10)))
      halt(SECT, cc, ": the recomputed score does not reproduce the committed ",
           "scores_per_patient.csv values")
    message(sprintf("  ..  %s  %-10s recomputed score reproduces the committed values (%d patients)",
                    SECT, cc, length(kk)))
    r <- rp[rp$cohort == cc, ]
    r <- r[order(r$patient, r$aliquot), ]
    r <- r[!duplicated(r$patient), ]
    y <- setNames(r$pY705, r$patient)
    k <- intersect(names(s), names(y))
    if (length(k) < MIN_N) halt(SECT, cc, ": joined n = ", length(k), ", below ", MIN_N)
    res[[cc]] <- rbind(
      cbind(cohort = cc, cor_report(s[k], y[k], "pearson",
                                    paste0("score_140 vs RPPA STAT3_pY705 (", cc, ")"))),
      cbind(cohort = cc, cor_report(s[k], y[k], "spearman",
                                    paste0("score_140 vs RPPA STAT3_pY705 (", cc, ")"))))
    per[[cc]] <- data.frame(cohort = cc, patient = k, score_140 = unname(s[k]),
                            rppa_pY705 = unname(y[k]), label = TAG,
                            stringsAsFactors = FALSE)
    message(sprintf("  ..  %s  %-10s n = %3d | Pearson %+.4f | Spearman %+.4f",
                    SECT, cc, length(k), res[[cc]]$r[1], res[[cc]]$r[2]))
  }
  cr <- do.call(rbind, res); rownames(cr) <- NULL
  pp <- do.call(rbind, per)

  # pooled across qualifying cohorts, on the Fisher z scale
  pe <- cr[cr$method == "pearson", ]
  if (nrow(pe) > 1L) {
    z <- atanh(pe$r); vz <- 1/(pe$n - 3)
    mz <- fns$meta_one(z, vz, pe$cohort)
    pooled <- data.frame(
      cohort = "POOLED", comparison = "score_140 vs RPPA STAT3_pY705 (Fisher z, random effects)",
      method = "pearson_pooled", n = sum(pe$n), r = tanh(mz$est),
      ci_lo = tanh(mz$ci_lo_wald), ci_hi = tanh(mz$ci_hi_wald), p = mz$p_wald,
      label = TAG, stringsAsFactors = FALSE)
    het <- data.frame(statistic = "pooled Fisher z, random effects (REML)",
                      k = mz$k, n_total = sum(pe$n), z_est = mz$est,
                      r_backtransformed = tanh(mz$est),
                      ci_lo_wald = tanh(mz$ci_lo_wald), ci_hi_wald = tanh(mz$ci_hi_wald),
                      ci_lo_hksj = tanh(mz$ci_lo_hksj), ci_hi_hksj = tanh(mz$ci_hi_hksj),
                      tau2 = mz$tau2, I2 = mz$I2, Q = mz$Q, Q_df = mz$Q_df, Q_p = mz$Q_p,
                      pi_lo = tanh(mz$pi_lo), pi_hi = tanh(mz$pi_hi), label = TAG,
                      note = paste0("EXPLORATORY, POST-HOC. tanh of the pooled z is the ",
                                    "back-transformed pooled correlation, not an unbiased ",
                                    "mean of the per-cohort r. Wald interval quoted, as in ",
                                    "the colorectal pooling. Cohorts differ in tissue and in ",
                                    "n, and the RPPA antibody is one phosphosite assay."),
                      stringsAsFactors = FALSE)
    write.csv(het, file.path(OUTDIR, "exploratory_rppa_pooled_heterogeneity.csv"),
              row.names = FALSE)
    cr <- rbind(cr, pooled[, names(cr)])
    message(sprintf("  ok  %s  pooled Pearson (Fisher z, %d cohorts, n = %d): %.4f (%.4f to %.4f), I2 %.1f%%",
                    SECT, mz$k, sum(pe$n), tanh(mz$est), tanh(mz$ci_lo_wald),
                    tanh(mz$ci_hi_wald), mz$I2))
  }
  cr$r <- round(cr$r, 4); cr$ci_lo <- round(cr$ci_lo, 4); cr$ci_hi <- round(cr$ci_hi, 4)
  write.csv(cr, file.path(OUTDIR, "exploratory_rppa_correlations.csv"), row.names = FALSE)
  write.csv(pp, file.path(OUTDIR, "exploratory_rppa_per_patient.csv"), row.names = FALSE)

  writeLines(c(
    "# EXPLORATORY, POST-HOC -- not registered, not an amendment",
    paste0("run_date: ", Sys.Date()),
    "",
    "1. Colorectal pooling of TCGA-COAD + TCGA-READ + GSE39582 (M1 score log-HR).",
    "   NOT the registered discovery meta-analysis. Amendment 16: validation cohorts",
    "   are not meta-analysed with discovery. Reported as a post-hoc precision exercise.",
    "",
    "2. TCGA RPPA STAT3_pY705 concordance within the discovery cohorts.",
    "   source   : GDC API, data_category 'Proteome Profiling', open access",
    "   antibody : AGID00388 / catalog 9131 / peptide_target 'STAT3_pY705'",
    "   NOTE     : 'Stat3' (AGID00185) in the same files is TOTAL STAT3, a different",
    "              antibody, and is NOT used here.",
    readLines("data/validation/RPPA/provenance.txt"),
    "",
    "PRIOR-KNOWLEDGE DISCLOSURE: the author's prior ESMO Asia work used TCGA RPPA",
    "STAT3_pY705 in TCGA-CHOL. This data source is not new to the author. TCGA-CHOL",
    "is excluded here (not meta-eligible, Amendment 8; and the prior-use ground).",
    "",
    paste0("script_md5: ", unname(tools::md5sum("13_rppa_exploratory.R"))),
    paste0("gene list: ", fns$GENE_LIST_PRIMARY$list_id, " (", fns$digest_genes(g140), ")"),
    "NO survival model fitted here. NO amendment written. No registered quantity altered."),
    file.path(OUTDIR, "exploratory_provenance.txt"))

  message("\n-- colorectal pooling (EXPLORATORY, POST-HOC) --")
  print(crc_out[, c("k","n_total","events_total","est","HR","HR_lo_wald","HR_hi_wald",
                    "tau2","I2","Q","Q_p")], row.names = FALSE)
  message("\n-- RPPA concordance (EXPLORATORY, POST-HOC) --")
  print(cr[, c("cohort","method","n","r","ci_lo","ci_hi","p")], row.names = FALSE)
  message("\nEXPLORATORY complete. HARD STOP: no survival model, no amendment.")
}
