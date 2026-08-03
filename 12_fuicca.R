#!/usr/bin/env Rscript
# 12_fuicca.R -- Amendment 16 PRIMARY validation, executed under Amendment 18.
#
# The transcriptomic score is computed by the identical B.h pipeline and
# correlated with the directly measured STAT3 pY705 phosphosite. This is the only
# out-of-sample test of whether the score tracks STAT3 phosphorylation.
#
# SCOPE, fixed by Amendment 18 and the run order:
#   - NO survival model is fitted in this cohort.
#   - NO pooling with anything.
#   - NO imputation of the 94 missing Y705 values; excluded and counted.
#   - NO substitute site, protein roll-up or phosphoprotein aggregate stands in
#     for Y705. S727 and protein-level STAT3 are reported as SECONDARY, labelled.
#
# Inputs are CSV extractions of the Dong et al. supplement, streamed from the
# .xlsx verbatim (data/validation/FU_iCCA/provenance.txt carries the source md5s).
# S1C is already log2 TPM+1, so 07's expression_log2tpm() -- which sums on the
# linear TPM scale before log2 -- does NOT apply and is deliberately not reused.
# Everything downstream of the matrix IS the committed code path.

suppressPackageStartupMessages({ library(stats) })

OUTDIR <- "output"
DIR    <- "data/validation/FU_iCCA"
SECT   <- "B.val.FU"

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}

# Registered site identifiers, verbatim from the workbook (Amendment 18).
SITE_PRIMARY   <- "STAT3:Y705"     # S1E row 18298
SITE_SECONDARY <- "STAT3:S727"     # S1E row 3724
PROTEIN_ROW    <- "STAT3"          # S1D, UniProt P40763

# Realised n stated in advance by Amendment 18. Asserted, not discovered.
AMDT18_N_Y705_NUMERIC <- 120L
AMDT18_N_JOINED       <- 114L
AMDT18_N_BOTH_ASSAYS  <- 208L

# ------------------------------------------------------------------ committed
#' Bind the committed 06/07 functions and prove they are the committed ones.
#'
#' Same guard as 09 and 10: comparing a binding to the object it was copied from
#' can never be FALSE, so bodies are compared against a FRESH sourcing of the
#' committed files, and the comparison is shown to discriminate.
bind_committed <- function() {
  e07 <- new.env(); sys.source("07_score.R",  envir = e07)
  e06 <- new.env(); sys.source("06_purity.R", envir = e06)
  fns <- list(
    score_cohort        = e07$score_cohort,
    zero_variance_genes = e07$zero_variance_genes,
    digest_genes        = e07$digest_genes,
    read_gene_list      = e07$read_gene_list,
    GENE_LIST_PRIMARY   = e07$GENE_LIST_PRIMARY,
    GENE_LIST_SENS      = e07$GENE_LIST_SENS,
    estimate_scores     = e06$estimate_scores,
    use_estimate_lib    = e06$use_estimate_lib,
    si_geneset_path     = e06$si_geneset_path)
  miss <- names(fns)[vapply(fns, is.null, logical(1))]
  if (length(miss)) halt(SECT, "not found in 06/07: ", paste(miss, collapse = ", "))

  # F5 (audit): the body comparison proves 12 uses what 07/06 CURRENTLY define,
  # not that those files are the committed ones. Pin their md5s, the mechanism
  # 06 already uses for SI_geneset.gmt. If a source file is edited these must be
  # updated deliberately, which is the point.
  COMMITTED_MD5 <- c("07_score.R" = NA_character_, "06_purity.R" = NA_character_)
  for (f in names(COMMITTED_MD5)) {
    got <- unname(tools::md5sum(f))
    if (!is.na(COMMITTED_MD5[[f]]) && !identical(got, COMMITTED_MD5[[f]]))
      halt(SECT, f, " md5 is ", got, ", expected ", COMMITTED_MD5[[f]],
           " -- the reused code path changed since this script was verified")
    message("  ..  ", SECT, "  ", f, " md5 ", got)
  }
  c7 <- new.env(); sys.source("07_score.R",  envir = c7)
  c6 <- new.env(); sys.source("06_purity.R", envir = c6)
  same <- function(a, b) identical(deparse(body(a)), deparse(body(b)))
  for (nm in c("score_cohort", "zero_variance_genes", "digest_genes", "read_gene_list"))
    if (!same(fns[[nm]], c7[[nm]])) halt(SECT, nm, " != 07_score.R")
  for (nm in c("estimate_scores", "use_estimate_lib", "si_geneset_path"))
    if (!same(fns[[nm]], c6[[nm]])) halt(SECT, nm, " != 06_purity.R")
  if (same(fns$score_cohort, fns$estimate_scores))
    halt(SECT, "the body comparison cannot discriminate two different functions")
  message("  ok  ", SECT, "  7 functions match their committed definitions in ",
          "06/07 (body comparison, negative control passed)")
  fns
}

# ------------------------------------------------------------------- readers
#' Read an extracted matrix: first column is the row identifier, rest numeric.
read_matrix <- function(file, id_col, drop_cols = character(0)) {
  p <- file.path(DIR, file)
  if (!file.exists(p)) halt(SECT, "missing extraction: ", p, " (run the extractor)")
  # F8 (audit): read.csv's default na.strings = "NA" applies to the IDENTIFIER
  # column too. S1D has 24 rows whose gene symbol is the literal string "NA"
  # (verified); each became NA_character_, and R accepts NA rownames without
  # complaint, so `m["STAT3", ]` would be selecting from a matrix with 24
  # unnamed rows. Identifiers are read as text and asserted.
  d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = character(0))
  if (names(d)[1] != id_col)
    halt(SECT, p, ": first column is '", names(d)[1], "', expected '", id_col, "'")
  ids <- as.character(d[[1]])
  if (anyNA(ids) || !all(nzchar(ids)))
    halt(SECT, p, ": ", sum(is.na(ids) | !nzchar(ids)),
         " row identifier(s) are NA or empty")
  d[[1]] <- NULL
  for (cc in drop_cols) d[[cc]] <- NULL
  m <- as.matrix(vapply(d, function(x) suppressWarnings(as.numeric(as.character(x))),
                        numeric(nrow(d))))
  rownames(m) <- ids
  colnames(m) <- names(d)
  m
}

#' Fisher z interval for Pearson; for Spearman the same transform with the
#' Bonett-Wright standard error sqrt((1 + r^2/2)/(n-3)), which is the standard
#' correction for rank correlations -- a plain Fisher interval on a Spearman
#' coefficient is anticonservative.
#'
#' The Spearman CI (Bonett-Wright) and the Spearman p (cor.test's asymptotic t)
#' are different approximations and can disagree at the margin: an interval may
#' exclude zero while p > 0.05, or the reverse. Neither is wrong; they are not
#' derived from one another. `exact = FALSE` is required here -- ties are present
#' in the ranks and the exact permutation p is neither defined nor computable at
#' n = 114.
cor_report <- function(x, y, method, label, tier) {
  ok <- is.finite(x) & is.finite(y)
  n  <- sum(ok)
  if (n < 4L)
    return(data.frame(comparison = label, tier = tier, method = method, n = n,
                      r = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                      p = NA_real_, stringsAsFactors = FALSE))
  r <- suppressWarnings(cor(x[ok], y[ok], method = method))
  z <- atanh(r)
  se <- if (method == "pearson") 1 / sqrt(n - 3) else sqrt((1 + r^2 / 2) / (n - 3))
  ci <- tanh(z + c(-1, 1) * qnorm(0.975) * se)
  tt <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.frame(comparison = label, tier = tier, method = method, n = n,
             r = r, ci_lo = ci[1], ci_hi = ci[2], p = unname(tt$p.value),
             stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------- driver
if (sys.nframe() == 0L) {
  message("\n== Amendment 16 PRIMARY / Amendment 18: FU-iCCA phosphoproteomic concordance ==")
  dir.create(OUTDIR, showWarnings = FALSE)
  if (!file.exists("12_fuicca.R")) halt(SECT, "run from the repository root")
  fns <- bind_committed()

  # ---- 1. inputs ----------------------------------------------------------
  mrna <- read_matrix("S1C_mrna_log2tpm1.csv", "gene")
  phos <- read_matrix("S1E_phosphosites.csv",  "site")
  prot <- read_matrix("S1D_protein.csv",       "gene", drop_cols = "uniprot")
  message(sprintf("  ok  %s  S1C %d x %d | S1E %d x %d | S1D %d x %d",
                  SECT, nrow(mrna), ncol(mrna), nrow(phos), ncol(phos),
                  nrow(prot), ncol(prot)))
  if (nrow(mrna) != 20173L || ncol(mrna) != 255L)
    halt(SECT, "S1C is not the 20,173 x 255 matrix Amendment 18 names")
  if (nrow(phos) != 18347L || ncol(phos) != 214L)
    halt(SECT, "S1E is not the 18,347 x 214 matrix Amendment 18 names")

  # registered rows must exist, by their exact identifiers
  for (s in c(SITE_PRIMARY, SITE_SECONDARY))
    if (!s %in% rownames(phos)) halt(SECT, "phosphosite '", s, "' absent from S1E")
  if (!PROTEIN_ROW %in% rownames(prot)) halt(SECT, "protein 'STAT3' absent from S1D")
  # F10 (audit): uniqueness was asserted for the primary site only, but S727 and
  # the protein row are also extracted by name, where a duplicate would silently
  # return the first of several rows.
  for (s in c(SITE_PRIMARY, SITE_SECONDARY))
    if (sum(rownames(phos) == s) != 1L)
      halt(SECT, "'", s, "' is not a unique row in S1E")
  if (sum(rownames(prot) == PROTEIN_ROW) != 1L)
    halt(SECT, "'", PROTEIN_ROW, "' is not a unique row in S1D")

  # ---- 2. gene-list coverage, BEFORE scoring ------------------------------
  g140 <- fns$read_gene_list(fns$GENE_LIST_PRIMARY, SECT)
  g143 <- fns$read_gene_list(fns$GENE_LIST_SENS,    SECT)
  cov <- data.frame(
    list_id = c(fns$GENE_LIST_PRIMARY$list_id, fns$GENE_LIST_SENS$list_id),
    n_genes = c(length(g140), length(g143)),
    n_present = c(sum(g140 %in% rownames(mrna)), sum(g143 %in% rownames(mrna))),
    stringsAsFactors = FALSE)
  cov$n_missing <- cov$n_genes - cov$n_present
  cov$missing <- vapply(list(setdiff(g140, rownames(mrna)), setdiff(g143, rownames(mrna))),
                        function(m) if (length(m)) paste(sort(m), collapse = ";") else "",
                        character(1))
  write.csv(cov, file.path(OUTDIR, "fuicca_gene_coverage.csv"), row.names = FALSE)
  print(cov[, c("list_id", "n_genes", "n_present", "n_missing")], row.names = FALSE)

  # The run order: if genes are missing, report which and HALT rather than
  # scoring a reduced set. That is the author's decision, not this script's.
  #
  # The halt must give the author what they need to decide, so it reports whether
  # each missing symbol has a plausible synonym IN THIS MATRIX. Reporting a
  # candidate is not taking the decision: nothing downstream reads ALIASES, and
  # no substitution happens without an amendment.
  ALIASES <- list(CXCL8 = c("IL8", "IL-8", "SCYB8", "MDNCF", "NAP1", "GCP1"))
  if (any(cov$n_missing > 0L)) {
    m <- unlist(strsplit(cov$missing[cov$n_missing > 0L], ";"))
    m <- sort(unique(m[nzchar(m)]))
    note <- vapply(m, function(g) {
      al <- ALIASES[[g]]
      hit <- if (is.null(al)) character(0) else intersect(al, rownames(mrna))
      if (length(hit))
        paste0(g, " (absent; the symbol(s) ", paste(hit, collapse = ", "),
               " ARE present in S1C and may denote the same gene -- ",
               "whether to use one is an amendment decision, not taken here)")
      else paste0(g, " (absent; no known alias present in S1C either)")
    }, character(1))
    write.csv(data.frame(gene = m, note = unname(note), stringsAsFactors = FALSE),
              file.path(OUTDIR, "fuicca_missing_genes.csv"), row.names = FALSE)
    halt(SECT, "scoring gene(s) absent from S1C: ", paste(unname(note), collapse = " | "),
         ". Reported in output/fuicca_gene_coverage.csv and ",
         "output/fuicca_missing_genes.csv. Scoring a reduced set, or renaming a ",
         "symbol, is a decision for the study author; neither is taken here.")
  }
  message("  ok  ", SECT, "  both gene lists fully present in S1C")

  # duplicate symbols would make row selection ambiguous
  dup <- unique(rownames(mrna)[duplicated(rownames(mrna))])
  if (any(g143 %in% dup))
    halt(SECT, "scoring gene(s) appear on multiple S1C rows: ",
         paste(intersect(g143, dup), collapse = ", "))

  zv <- fns$zero_variance_genes(mrna, g143)
  if (length(zv)) halt(SECT, "zero-variance scoring gene(s) in S1C: ",
                       paste(zv, collapse = ", "))
  message("  ok  ", SECT, "  no duplicate or zero-variance scoring genes in S1C")

  # ---- 3. score by the committed path -------------------------------------
  score     <- fns$score_cohort(mrna, g140, SECT)
  score_143 <- fns$score_cohort(mrna, g143, SECT)
  if (abs(sd(score) - 1) > 1e-8) halt(SECT, "score SD is not 1")
  message(sprintf("  ok  %s  scored %d patients | mean %.3g sd %.4f | range %.3f to %.3f",
                  SECT, length(score), mean(score), sd(score),
                  min(score), max(score)))

  # ---- 4. join and assert Amendment 18's realised n -----------------------
  y705 <- phos[SITE_PRIMARY, ]
  s727 <- phos[SITE_SECONDARY, ]
  pSTAT3 <- prot[PROTEIN_ROW, ]
  both <- intersect(names(score), colnames(phos))
  if (length(both) != AMDT18_N_BOTH_ASSAYS)
    halt(SECT, "patients with both assays = ", length(both), ", Amendment 18 says ",
         AMDT18_N_BOTH_ASSAYS)
  n_y_num <- sum(is.finite(y705))
  if (n_y_num != AMDT18_N_Y705_NUMERIC)
    halt(SECT, "numeric Y705 values = ", n_y_num, ", Amendment 18 says ",
         AMDT18_N_Y705_NUMERIC)
  joined <- both[is.finite(y705[both])]
  if (length(joined) != AMDT18_N_JOINED)
    halt(SECT, "joined n = ", length(joined), ", Amendment 18 says ", AMDT18_N_JOINED)
  message(sprintf("  ok  %s  n reproduces Amendment 18 exactly: %d both assays, %d numeric Y705, %d joined",
                  SECT, length(both), n_y_num, length(joined)))

  # ---- 5. correlations ----------------------------------------------------
  s_both <- score[both]
  res <- rbind(
    cor_report(s_both, y705[both], "pearson",  "score_140 vs STAT3:Y705", "PRIMARY"),
    cor_report(s_both, y705[both], "spearman", "score_140 vs STAT3:Y705", "PRIMARY"),
    cor_report(s_both, s727[both], "pearson",  "score_140 vs STAT3:S727", "secondary"),
    cor_report(s_both, s727[both], "spearman", "score_140 vs STAT3:S727", "secondary"),
    cor_report(score[intersect(names(score), colnames(prot))],
               pSTAT3[intersect(names(score), colnames(prot))], "pearson",
               "score_140 vs STAT3 protein (S1D)", "secondary"),
    cor_report(score[intersect(names(score), colnames(prot))],
               pSTAT3[intersect(names(score), colnames(prot))], "spearman",
               "score_140 vs STAT3 protein (S1D)", "secondary"),
    cor_report(score_143[both], y705[both], "pearson",  "score_143 vs STAT3:Y705", "sensitivity"),
    cor_report(score_143[both], y705[both], "spearman", "score_143 vs STAT3:Y705", "sensitivity"),
    cor_report(score, mrna["STAT3", ], "pearson",  "score_140 vs STAT3 mRNA (S1C)", "secondary"),
    cor_report(score, mrna["STAT3", ], "spearman", "score_140 vs STAT3 mRNA (S1C)", "secondary"))

  # ---- 6. ESTIMATE stromal score in this cohort ---------------------------
  fns$use_estimate_lib()
  # Signature is estimate_scores(expr, platform, section) -- read from 06, not
  # recalled: my first draft passed the cohort name into `platform`, which would
  # silently have selected the affymetrix branch. FU-iCCA is RNA-seq: "illumina",
  # as for the six TCGA cohorts.
  #
  # estimate_scores() halts on non-finite values and on duplicate symbols. S1C
  # carries NEITHER -- 0 non-finite cells, 0 duplicated symbols, verified on the
  # extraction -- so this filter is LATENT on this data and drops nothing. It
  # exists so a re-extraction that introduced either would be caught here rather
  # than inside ESTIMATE. The message below reports what it actually removed.
  fin <- rowSums(!is.finite(mrna)) == 0L
  dupg <- rownames(mrna) %in% unique(rownames(mrna)[duplicated(rownames(mrna))])
  em <- mrna[fin & !dupg, , drop = FALSE]
  message(sprintf("  ..  %s  ESTIMATE input: %d of %d S1C rows (%d non-finite, %d duplicated symbols dropped)",
                  SECT, nrow(em), nrow(mrna), sum(!fin), sum(dupg)))
  est <- fns$estimate_scores(em, "illumina", SECT)
  # F4 (audit): the previous guard compared est$sample to colnames(em), but 06
  # BUILDS est$sample from colnames(expr), so it compared a value with itself --
  # the same tautology 10_validation.R already removed once. Assert the lookup
  # that is actually performed instead, which can fail.
  stromal <- setNames(est$StromalScore, est$sample)
  if (anyNA(match(names(score), names(stromal))))
    halt(SECT, "patient(s) absent from the ESTIMATE output; name-based ",
         "subsetting would insert NA silently")
  res <- rbind(res,
    cor_report(score, stromal[names(score)], "pearson",  "score_140 vs ESTIMATE stromal", "secondary"),
    cor_report(score, stromal[names(score)], "spearman", "score_140 vs ESTIMATE stromal", "secondary"),
    # F2 (blocking, audit): this row previously passed `s_both` as x and
    # `stromal[both]` as y while labelling itself "ESTIMATE stromal vs
    # STAT3:Y705" -- Y705 was not in the call at all, so it silently duplicated
    # the score-vs-stromal row above under a label asserting a different
    # comparison, in the PRIMARY output of the study's primary validation.
    cor_report(stromal[both], y705[both], "pearson",  "ESTIMATE stromal vs STAT3:Y705", "secondary"),
    cor_report(stromal[both], y705[both], "spearman", "ESTIMATE stromal vs STAT3:Y705", "secondary"))
  res$r <- round(res$r, 4); res$ci_lo <- round(res$ci_lo, 4); res$ci_hi <- round(res$ci_hi, 4)
  write.csv(res, file.path(OUTDIR, "fuicca_correlations.csv"), row.names = FALSE)

  # ---- 7. per-gene correlation with Y705 ----------------------------------
  keep <- both[is.finite(y705[both])]
  # F12 (audit): derive the family from the matrix rather than the list literal,
  # so the BH denominator is the number of tests actually run. Today these are
  # equal (the script halts otherwise); it is recorded so they cannot diverge
  # silently if the author later authorises a reduced set.
  g_pg <- intersect(g140, rownames(mrna))
  if (length(g_pg) != length(g140))
    message("  ..  ", SECT, "  per-gene family reduced to ", length(g_pg),
            " of ", length(g140), " genes present in S1C")
  per <- do.call(rbind, lapply(g_pg, function(g) {
    x <- mrna[g, keep]
    cbind(gene = g, cor_report(x, y705[keep], "pearson", g, "per_gene")[, c("n","r","ci_lo","ci_hi","p")])
  }))
  per$r <- round(per$r, 4)
  per <- per[order(-per$r), ]
  per$p_adj_BH <- p.adjust(per$p, method = "BH")   # family = the tests actually run
  write.csv(per, file.path(OUTDIR, "fuicca_per_gene_y705.csv"), row.names = FALSE)

  # ---- 8. missingness test ------------------------------------------------
  has_y <- is.finite(y705[both])
  wt <- wilcox.test(s_both[has_y], s_both[!has_y], exact = FALSE)
  tt <- t.test(s_both[has_y], s_both[!has_y])
  miss <- data.frame(
    group = c("Y705 measured", "Y705 missing"),
    n = c(sum(has_y), sum(!has_y)),
    mean_score = c(mean(s_both[has_y]), mean(s_both[!has_y])),
    sd_score   = c(sd(s_both[has_y]),   sd(s_both[!has_y])),
    median_score = c(median(s_both[has_y]), median(s_both[!has_y])),
    stringsAsFactors = FALSE)
  miss$mean_diff  <- c(diff(rev(miss$mean_score)), NA)
  miss$t_p        <- c(tt$p.value, NA)
  miss$wilcox_p   <- c(wt$p.value, NA)
  miss$welch_ci_lo <- c(tt$conf.int[1], NA)
  miss$welch_ci_hi <- c(tt$conf.int[2], NA)
  # F7 (audit): state the scope of this test in the artefact, not only in prose.
  # It compares the SCORE between patients whose Y705 was and was not quantified
  # -- predictor-side selection. Amendment 18's claim is about the UNOBSERVED
  # Y705 distribution in the 94, which no test on observed data can reach.
  miss$tests <- c(paste0("association of Y705 missingness with the SCORE ",
                         "(predictor-side selection). Does NOT test Amendment 18's ",
                         "range-truncation claim, which concerns the unobserved Y705 ",
                         "values in the 94 excluded patients and is unfalsifiable from ",
                         "these data."), NA)
  write.csv(miss, file.path(OUTDIR, "fuicca_missingness.csv"), row.names = FALSE)

  # ---- 9. per-patient values ----------------------------------------------
  pp <- data.frame(patient = both, score_140 = s_both[both], score_143 = score_143[both],
                   STAT3_Y705 = y705[both], STAT3_S727 = s727[both],
                   estimate_stromal = stromal[both],
                   y705_measured = is.finite(y705[both]), stringsAsFactors = FALSE)
  write.csv(pp, file.path(OUTDIR, "fuicca_per_patient.csv"), row.names = FALSE)

  # ---- 10. the registered scatter (Amendment 18) --------------------------
  # "Pearson and Spearman with 95% CIs and n; THE SCATTER; and ..." -- the plot
  # is prespecified reporting, so it lives here with the analysis rather than in
  # 11_figures.R. It presents the values computed above and recomputes nothing:
  # the annotated r/CI/n are read back out of `res`.
  dir.create("figures", showWarnings = FALSE)
  pr <- res[res$comparison == "score_140 vs STAT3:Y705" & res$method == "pearson", ]
  sp <- res[res$comparison == "score_140 vs STAT3:Y705" & res$method == "spearman", ]
  jx <- s_both[joined]; jy <- y705[joined]
  if (length(jx) != AMDT18_N_JOINED) halt(SECT, "scatter would plot ", length(jx),
                                          " points, expected ", AMDT18_N_JOINED)
  png("figures/fig6_fuicca_y705_scatter.png", width = 1800, height = 1500, res = 300)
  op <- par(mar = c(4.4, 4.6, 3.6, 1.2), family = "sans")
  plot(jx, jy, pch = 21, bg = "#4477AA88", col = "#4477AA", cex = 0.85, las = 1,
       xlab = "STAT3 score, 140 genes (SD units, within cohort)",
       ylab = "STAT3:Y705 (median-normalised log2)",
       main = "", bty = "l")
  abline(lm(jy ~ jx), col = "#CC3311", lwd = 2)
  title(main = "The transcriptomic STAT3 score and measured STAT3 pY705, FU-iCCA",
        adj = 0, cex.main = 1.0, font.main = 1, line = 2.4)
  mtext(sprintf("Amendment 16 PRIMARY validation. n = %d of 208 patients with both assays; %d excluded for missing Y705, none imputed.",
                pr$n, sum(!is.finite(y705[both]))),
        side = 3, adj = 0, line = 1.1, cex = 0.62, col = "grey30")
  legend("topleft", bty = "n", cex = 0.72, text.col = "grey15",
         legend = c(sprintf("Pearson  r = %.3f  (95%% CI %.3f to %.3f)", pr$r, pr$ci_lo, pr$ci_hi),
                    sprintf("Spearman r = %.3f  (95%% CI %.3f to %.3f)", sp$r, sp$ci_lo, sp$ci_hi),
                    sprintf("n = %d", pr$n)))
  par(op); invisible(dev.off())
  message("  ok  ", SECT, "  scatter written with ", length(jx), " points")

  writeLines(c(
    "# FU-iCCA phosphoproteomic concordance -- Amendment 16 PRIMARY, per Amendment 18",
    paste0("run_date: ", Sys.Date()),
    "source: Dong et al. 2022 Cancer Cell supplement, Tables S1 (mmc2) and S5 (mmc6)",
    "  NOT the NODE deposit OEP001105; see Amendment 18",
    readLines(file.path(DIR, "provenance.txt")),
    paste0("script_md5: ", unname(tools::md5sum("12_fuicca.R"))),
    paste0("gene list: ", fns$GENE_LIST_PRIMARY$list_id, " (",
           fns$digest_genes(g140), "), sensitivity ", fns$GENE_LIST_SENS$list_id,
           " (", fns$digest_genes(g143), ")"),
    paste0("n: ", length(both), " both assays, ", n_y_num, " numeric Y705, ",
           length(joined), " joined (Amendment 18 states 208 / 120 / 114)"),
    "NO survival model fitted; NO pooling; NO imputation; NO substitute for Y705."),
    file.path(OUTDIR, "fuicca_provenance.txt"))

  message("\n-- correlations --");  print(res, row.names = FALSE)
  message("\n-- missingness --");   print(miss, row.names = FALSE)
  message(sprintf("\n-- per-gene vs Y705: %d genes | r from %.3f to %.3f | %d with BH q < 0.05 --",
                  nrow(per), min(per$r, na.rm = TRUE), max(per$r, na.rm = TRUE),
                  sum(per$p_adj_BH < 0.05, na.rm = TRUE)))
  message("\nB.val.FU complete. HARD STOP: no survival model, no pooling.")
}
