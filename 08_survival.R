#!/usr/bin/env Rscript
# 08_survival.R -- Part B: the PRIMARY RESULT.
#
# Four nested Cox models per cohort on an IDENTICAL complete-case set (B.j), the
# attenuation estimand and its paired bootstrap (B.k), and the meta-analysis
# (B.l).
#
#   M1: Surv(time, event) ~ score                                   1 parameter
#   M2: + age + sex + stage_group                                   5
#   M3: + purity                                                    6
#   M4: + stromal_score                                             7
#
# PRIMARY ESTIMAND: attenuation_total = beta(M2) - beta(M4), meta-analysed.
# Components: attenuation_purity = b2 - b3, attenuation_stroma = b3 - b4.
# prop_attenuated is SECONDARY (unstable near b2 = 0, undefined on sign flip).
#
# The interval is a paired nonparametric bootstrap over PATIENTS WITHIN COHORT,
# B = 2000, seed 20260731, refitting ALL FOUR models inside each resample so that
# b2 and b4 carry their correlation. Two independent SEs would overstate the
# width: the estimates come from nested models on identical data.
#
# STRUCTURE: pure functions taking a data.frame; driver guarded by sys.nframe().

suppressPackageStartupMessages({
  library(survival)
  library(metafor)
  library(withr)      # F12: the primary result's RNG stream depends on it
  library(car)        # F9: Fox-Monette GVIF, the correct VIF for a Cox fit
})

halt <- function(section, ...)
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)

assert_n <- function(observed, expected, section, what) {
  if (!identical(as.integer(observed), as.integer(expected)))
    halt(section, what, ": expected ", expected, ", observed ", observed)
  message(sprintf("  ok  %-10s %-56s = %d", section, what, observed))
  invisible(TRUE)
}

COHORTS      <- c("TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-ESCA",
                  "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL")
META_ELIGIBLE <- c("TCGA-COAD" = TRUE, "TCGA-READ" = TRUE, "TCGA-STAD" = TRUE,
                   "TCGA-ESCA" = TRUE, "TCGA-PAAD" = TRUE, "TCGA-LIHC" = TRUE,
                   "TCGA-CHOL" = FALSE)            # Amendment 8, descriptive only
# FIXED by registration, not derived from the fit. analysis_plan.md:1186-1187:
# "Parameter counts are fixed in B.j (M1 = 1, M2 = 5, M3 = 6, M4 = 7, with
# stage_group contributing two parameters for its three levels)." The audit
# proposed deriving the count from the realised fit so that a cohort missing the
# `missing` stage level would count one fewer parameter; that would be a
# deviation from the registered denominator, and it would make the EPV band
# depend on realised level occupancy. Registered counts are used; `use_sex` is
# the one adjustment the plan itself specifies (B.j: sex dropped below 10 of
# either sex), and it is applied only where a sex term exists (see model_formula).
MODEL_PARAMS <- c(M1 = 1L, M2 = 5L, M3 = 6L, M4 = 7L)
B_RESAMPLES  <- 2000L
SEED_BASE    <- 20260731L
EPV_FIT      <- 10   # >= 10 fit and pool; 5-10 fit, pool, flag; < 5 do not fit
EPV_FLOOR    <- 5

MODEL_RHS <- list(
  M1 = "score",
  M2 = "score + age + sex + stage_group",
  M3 = "score + age + sex + stage_group + purity",
  M4 = "score + age + sex + stage_group + purity + stromal_score")

# ============================================================ PURE FUNCTIONS

#' Realised EPV band for one model (Amendment 8), applied to REALISED events.
epv_band <- function(events, params) {
  e <- events / params
  if (!is.finite(e)) halt("B.j", "EPV is not finite (events=", events,
                          ", params=", params, "); the band cannot be assigned")
  else if (FALSE) NA_character_
  else if (e >= EPV_FIT)   "fit_and_pool"
  else if (e >= EPV_FLOOR) "fit_pool_flag_LOO"
  else                     "do_not_fit"
}

#' Model formula, dropping `sex` where the cohort cannot support it (B.j: fewer
#' than 10 patients of either sex). Dropping a term changes the parameter count,
#' so the count is returned with the formula rather than assumed.
model_formula <- function(model, use_sex) {
  rhs <- MODEL_RHS[[model]]
  if (!use_sex) rhs <- sub("\\s*\\+\\s*sex", "", rhs)
  # F6: subtract only where a sex term was actually present. M1's RHS is `score`
  # alone, so an unconditional subtraction gave params = 0 -> EPV = Inf -> band NA
  # -> stop condition 5 silently skipped for M1 in any sex-dropped cohort.
  drop1 <- !use_sex && grepl("sex", MODEL_RHS[[model]], fixed = TRUE)
  list(formula = stats::as.formula(paste("Surv(time, event) ~", rhs)),
       params  = MODEL_PARAMS[[model]] - if (drop1) 1L else 0L)
}

#' Complete-case set across ALL model-4 covariates (B.j).
#'
#' Applied ONCE, to the M4 covariate set, so all four models are fitted on an
#' identical patient set and the attenuation is not confounded by changing n.
complete_case_set <- function(d, use_sex) {
  need <- c("time", "event", "score", "age", "stage_group", "purity", "stromal_score")
  if (use_sex) need <- c(need, "sex")
  miss <- setdiff(need, names(d))
  if (length(miss)) halt("B.j", "missing column(s): ", paste(miss, collapse = ", "))
  keep <- stats::complete.cases(d[, need, drop = FALSE])
  d[keep, , drop = FALSE]
}

#' Fit one model and return the SCORE coefficient only.
#'
#' Convergence and separation are checked rather than assumed: coxph can return a
#' finite object with a divergent coefficient, and stop condition 5 requires those
#' to surface.
fit_one <- function(d, model, use_sex, section = "B.j") {
  mf <- model_formula(model, use_sex)
  f  <- try(survival::coxph(mf$formula, data = d, x = TRUE, ties = "efron"),
            silent = TRUE)
  if (inherits(f, "try-error"))
    return(list(ok = FALSE, reason = paste("coxph error:", attr(f, "condition")$message)))
  b  <- stats::coef(f)["score"]
  se <- sqrt(diag(stats::vcov(f)))["score"]
  bad <- !is.finite(b) || !is.finite(se) || abs(b) > 20 || se > 20
  list(ok = !bad,
       reason = if (bad) "non-finite or divergent score coefficient (separation)" else NA_character_,
       beta = unname(b), se = unname(se), fit = f, params = mf$params,
       n = nrow(d), events = sum(d$event == 1L))
}

#' All four models on one cohort's complete-case set.
fit_cohort <- function(d, use_sex, section = "B.j") {
  lapply(setNames(names(MODEL_PARAMS), names(MODEL_PARAMS)),
         function(m) fit_one(d, m, use_sex, section))
}

#' The attenuation estimand and its components (B.k).
attenuation <- function(b2, b3, b4) {
  c(attenuation_total  = b2 - b4,
    attenuation_purity = b2 - b3,
    attenuation_stroma = b3 - b4,
    # SECONDARY: unstable when b2 is near zero, undefined on a sign flip.
    # B.k, analysis_plan.md:962-963 -- "reported only if beta2 != 0 AND
    # sign(beta2)==sign(beta4)". The sign half was missing: on the realised data
    # COAD (b2=+0.0046, b4=-0.0713) would print prop_attenuated = 16.52 and CHOL
    # 1.08, both from a sign flip the registration says makes it undefined.
    prop_attenuated    = if (is.finite(b2) && abs(b2) > 1e-8 &&
                             is.finite(b4) && sign(b2) == sign(b4))
                           (b2 - b4) / b2 else NA_real_)
}

#' Paired nonparametric bootstrap over patients within one cohort.
#'
#' All four models are refitted inside each resample, so b2 and b4 are paired and
#' their difference carries the correlation between them. Resamples in which any
#' model fails to fit are recorded and dropped, never silently treated as zero.
bootstrap_cohort <- function(d, use_sex, cohort_index, B = B_RESAMPLES,
                             section = "B.k") {
  n <- nrow(d)
  reps <- withr::with_seed(SEED_BASE + cohort_index, {
    lapply(seq_len(B), function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      db  <- d[idx, , drop = FALSE]
      # A resample can lose a factor level entirely; refitting with a one-level
      # factor errors, which is a dropped resample, not a result.
      fits <- try(lapply(names(MODEL_PARAMS), function(m) {
        mf <- model_formula(m, use_sex)
        survival::coxph(mf$formula, data = db, x = FALSE, ties = "efron")
      }), silent = TRUE)
      if (inherits(fits, "try-error")) return(NULL)
      bb <- vapply(fits, function(f) unname(stats::coef(f)["score"]), numeric(1))
      if (any(!is.finite(bb)) || any(abs(bb) > 20)) return(NULL)
      # Count resamples in which ANY coefficient (not just the score) diverges or
      # is aliased. In READ -- 165 patients, 36 events, and a `missing` stage cell
      # of 9 patients with 1 event -- coxph warns "coefficient may be infinite" on
      # nearly every resample, always for stage_groupmissing. Diagnosed rather
      # than suppressed: across 300 test resamples ZERO had any |coef| > 20 and
      # the score coefficient stayed within [-0.79, 0.90], so the estimand is
      # unaffected. The count is carried out so a future cohort where this touches
      # the score is visible instead of silent.
      allc <- unlist(lapply(fits, stats::coef))
      attr_flag <- any(is.na(allc)) || any(abs(allc[!is.na(allc)]) > 20)
      c(b1 = bb[1], b2 = bb[2], b3 = bb[3], b4 = bb[4],
        attenuation(bb[2], bb[3], bb[4]), nuisance_unstable = as.numeric(attr_flag))
    })
  })
  ok <- !vapply(reps, is.null, logical(1))
  R <- do.call(rbind, reps[ok])
  list(reps = R, n_ok = sum(ok), n_failed = sum(!ok),
       n_nuisance_unstable = if (is.null(R)) 0L else sum(R[, "nuisance_unstable"] > 0))
}

#' Percentile CI from bootstrap replicates.
boot_ci <- function(v, probs = c(0.025, 0.975)) {
  v <- v[is.finite(v)]
  if (length(v) < 100L) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, probs))
}

#' Meta-analysis for one model (B.l).
#'
#' HKSJ (`test = "knha"`) is the PRIMARY interval -- prespecified because with six
#' cohorts the normal approximation is anticonservative. The Wald interval is
#' reported alongside for comparability, and fixed-effect as a sensitivity.
meta_one <- function(yi, vi, labels) {
  ok <- is.finite(yi) & is.finite(vi) & vi > 0
  if (sum(ok) < 2L) return(NULL)
  yi <- yi[ok]; vi <- vi[ok]; labels <- labels[ok]
  hk <- metafor::rma(yi = yi, vi = vi, method = "REML", test = "knha")
  wd <- metafor::rma(yi = yi, vi = vi, method = "REML")
  fe <- metafor::rma(yi = yi, vi = vi, method = "FE")
  pr <- try(predict(hk), silent = TRUE)
  list(k = length(yi), cohorts = labels,
       est = as.numeric(hk$beta), se_hksj = hk$se,
       ci_lo_hksj = hk$ci.lb, ci_hi_hksj = hk$ci.ub, p_hksj = hk$pval,
       ci_lo_wald = wd$ci.lb, ci_hi_wald = wd$ci.ub, p_wald = wd$pval,
       fe_est = as.numeric(fe$beta), fe_ci_lo = fe$ci.lb, fe_ci_hi = fe$ci.ub,
       tau2 = hk$tau2, I2 = hk$I2, H2 = hk$H2,
       Q = hk$QE, Q_df = hk$k - 1L, Q_p = hk$QEp,
       pi_lo = if (!inherits(pr, "try-error")) pr$pi.lb else NA_real_,
       pi_hi = if (!inherits(pr, "try-error")) pr$pi.ub else NA_real_)
}

#' cox.zph for every model, and the registered score:log(time) sensitivity.
ph_check <- function(fits, d, use_sex) {
  do.call(rbind, lapply(names(fits), function(m) {
    f <- fits[[m]]
    if (!isTRUE(f$ok)) return(NULL)
    z <- try(survival::cox.zph(f$fit), silent = TRUE)
    if (inherits(z, "try-error")) return(NULL)
    tab <- z$table
    data.frame(model = m,
               p_score  = if ("score" %in% rownames(tab)) tab["score", "p"] else NA_real_,
               p_global = tab["GLOBAL", "p"],
               violated_score  = !is.na(tab["score", "p"]) && tab["score", "p"] < 0.05,
               violated_global = tab["GLOBAL", "p"] < 0.05,
               stringsAsFactors = FALSE)
  }))
}

#' The registered PH sensitivity: score x log(time) interaction via tt().
#'
#' log(t) GUARD: 47 patients in the analysis set have time = 0 (same-day event or
#' censoring), and log(0) = -Inf would propagate silently into the interaction
#' term -- coxph does not error on it, it returns a fit whose tt coefficient is
#' meaningless. cox.zph itself uses the ordered event times, not log(t), so this
#' affects only the sensitivity. Following the standard convention, time is
#' shifted by the smallest positive time unit before logging; the shift is applied
#' identically in every cohort and reported.
ph_sensitivity <- function(d, model, use_sex) {
  mf  <- model_formula(model, use_sex)
  rhs <- paste(deparse(mf$formula[[3]]), collapse = " ")
  f2  <- stats::as.formula(paste("Surv(time, event) ~", rhs, "+ tt(score)"))
  fit <- try(survival::coxph(f2, data = d, ties = "efron",
                             tt = function(x, t, ...) x * log(t + 1)), silent = TRUE)
  if (inherits(fit, "try-error"))
    return(data.frame(model = model, beta_score = NA_real_, se_score = NA_real_,
                      beta_tt = NA_real_, se_tt = NA_real_, p_tt = NA_real_,
                      failed = TRUE, reason = as.character(attr(fit, "condition")$message),
                      stringsAsFactors = FALSE))
  b  <- stats::coef(fit); v <- sqrt(diag(stats::vcov(fit)))
  i  <- grep("^tt\\(score\\)$", names(b))
  data.frame(model = model,
             beta_score = unname(b["score"]), se_score = unname(v["score"]),
             beta_tt = if (length(i)) unname(b[i]) else NA_real_,
             se_tt   = if (length(i)) unname(v[i]) else NA_real_,
             p_tt    = if (length(i)) 2 * stats::pnorm(-abs(b[i] / v[i])) else NA_real_,
             failed = FALSE, reason = NA_character_,
             stringsAsFactors = FALSE)
}

#' VIFs for M4, computed on the model matrix so factor terms are handled.
#' F9: a Cox coefficient's precision comes from the risk-set-weighted information
#' matrix, not from the unweighted covariance of X, so an OLS VIF on the design
#' matrix does not describe the collinearity inflating M4's standard errors. The
#' Fox-Monette GVIF is computed from the coefficient covariance, which is the
#' correct analogue and what car::vif implements for coxph.
vif_m4 <- function(fit) {
  # car returns a 3-column matrix (GVIF, Df, GVIF^(1/(2*Df))) when any term has
  # Df > 1, and a bare named vector otherwise. Columns are taken by POSITION: the
  # third column's name is the literal "GVIF^(1/(2*Df))", so matching on "GVIF"
  # picks the wrong column and matching exactly returns NULL.
  v <- suppressWarnings(try(car::vif(fit), silent = TRUE))
  if (inherits(v, "try-error")) return(NULL)
  if (is.matrix(v) && ncol(v) >= 3L)
    data.frame(term = rownames(v), GVIF = round(v[, 1], 4), Df = as.integer(v[, 2]),
               GVIF_adj = round(v[, 3], 4), stringsAsFactors = FALSE)
  else
    data.frame(term = names(v), GVIF = round(as.numeric(v), 4),
               Df = 1L, GVIF_adj = round(sqrt(as.numeric(v)), 4),
               stringsAsFactors = FALSE)
}

# ==================================================================== driver
if (sys.nframe() == 0L) {

  OUTDIR <- "output"; dir.create(OUTDIR, showWarnings = FALSE)
  message("\n== B.j/B.k/B.l  survival, attenuation, meta-analysis ==")

  clin  <- read.csv(file.path(OUTDIR, "clinical_analysis_set.csv"), stringsAsFactors = FALSE)
  pur   <- read.csv(file.path(OUTDIR, "purity_per_patient.csv"), stringsAsFactors = FALSE)
  scr   <- read.csv(file.path(OUTDIR, "scores_per_patient.csv"), stringsAsFactors = FALSE)
  for (nm in c("clin", "pur", "scr")) if (!nrow(get(nm))) halt("B.j", nm, " is empty")

  # Join on barcode -- the only key that is unique within cohort and carries the
  # aliquot, not just the patient.
  d <- merge(clin, pur[, c("barcode", "purity", "purity_estimate", "CPE",
                           "purity_source", "StromalScore")], by = "barcode")
  d <- merge(d, scr[, c("barcode", "score", "score_143", "stromal_score",
                        "stromal_score_subscore")], by = "barcode")
  assert_n(nrow(d), nrow(clin), "B.j", "joined rows == clinical rows")
  if (anyDuplicated(d$barcode)) halt("B.j", "duplicate barcode after join")
  d$stage_group <- factor(d$stage_group, levels = c("I/II", "III/IV", "missing"))
  if (anyNA(d$stage_group)) halt("B.j", "stage_group has NA outside its explicit level")
  d$sex <- factor(d$sex)

  # ---- per cohort ---------------------------------------------------------
  res <- lapply(seq_along(COHORTS), function(i) {
    cc <- COHORTS[i]
    dc <- d[d$cohort == cc, , drop = FALSE]
    use_sex <- isTRUE(unique(dc$use_sex))
    cc_set  <- complete_case_set(dc, use_sex)
    cc_set$stage_group <- droplevels(cc_set$stage_group)
    ev    <- sum(cc_set$event == 1L)
    # Bands BEFORE fitting: Amendment 8's EPV < 5 rule is "model NOT fitted for
    # that cohort and that model" (panel_definition.md:294-296). Previously the
    # band was computed, written to the CSV and never consulted, so a do_not_fit
    # model was still fitted and still pooled.
    bands <- vapply(names(MODEL_PARAMS), function(m)
      epv_band(ev, model_formula(m, use_sex)$params), character(1))
    fits <- fit_cohort(cc_set, use_sex)
    for (m in names(fits)) if (identical(bands[[m]], "do_not_fit"))
      fits[[m]] <- list(ok = FALSE, reason = "EPV < 5: not fitted per Amendment 8")
    list(cohort = cc, data = cc_set, use_sex = use_sex, fits = fits,
         n = nrow(cc_set), events = ev, bands = bands,
         n_dropped_incomplete = nrow(dc) - nrow(cc_set),
         ph = ph_check(fits, cc_set, use_sex))
  })
  names(res) <- COHORTS

  # ---- stop condition 5: a cohort the EPV rule says to fit must converge ----
  for (cc in COHORTS) for (m in names(MODEL_PARAMS)) {
    b <- res[[cc]]$bands[[m]]; f <- res[[cc]]$fits[[m]]
    if (b %in% c("fit_and_pool", "fit_pool_flag_LOO") && !isTRUE(f$ok))
      halt("B.j", "STOP CONDITION 5: ", cc, " ", m, " (band ", b,
           ") failed to fit: ", f$reason)
  }
  message("  ok  B.j        every model the EPV rule says to fit converged")

  per_cohort <- do.call(rbind, lapply(COHORTS, function(cc) {
    r <- res[[cc]]
    do.call(rbind, lapply(names(MODEL_PARAMS), function(m) {
      f <- r$fits[[m]]
      data.frame(cohort = cc, model = m, n = r$n, events = r$events,
                 params = model_formula(m, r$use_sex)$params,
                 epv = round(r$events / model_formula(m, r$use_sex)$params, 2),
                 band = r$bands[[m]], use_sex = r$use_sex,
                 fitted = isTRUE(f$ok),
                 beta = if (isTRUE(f$ok)) round(f$beta, 6) else NA_real_,
                 se   = if (isTRUE(f$ok)) round(f$se, 6) else NA_real_,
                 HR   = if (isTRUE(f$ok)) round(exp(f$beta), 4) else NA_real_,
                 HR_lo = if (isTRUE(f$ok)) round(exp(f$beta - 1.96 * f$se), 4) else NA_real_,
                 HR_hi = if (isTRUE(f$ok)) round(exp(f$beta + 1.96 * f$se), 4) else NA_real_,
                 p = if (isTRUE(f$ok)) signif(2 * stats::pnorm(-abs(f$beta / f$se)), 4) else NA_real_,
                 meta_eligible = unname(META_ELIGIBLE[[cc]]),
                 stringsAsFactors = FALSE)
    }))
  }))
  # ---------------------------------------------------------------- B.n
  # "The per-cohort estimates are secondary and descriptive. Where they are
  # tested, p-values are adjusted across the six meta-analysed cohorts by
  # Benjamini-Hochberg FDR at q = 0.05, and both raw and adjusted values are
  # reported in the same table."
  #
  # FAMILY: the six META-ANALYSED cohorts, WITHIN model. CHOL is excluded from the
  # family by B.n's own text ("CHOL is reported descriptively and is not counted in
  # the multiplicity family") and receives NA, not an adjusted value computed over
  # a family it does not belong to.
  #
  # NO CORRECTION ACROSS M1-M4. B.n: they "are a prespecified nested sequence
  # addressing a single question, not four independent hypotheses; no correction is
  # applied across them, and this is stated." Adjusting within model rather than
  # over the whole 4x6 grid is that statement made operational.
  #
  # The primary inference is unaffected: B.n applies no correction to the pooled
  # estimate, which is "one estimand, one test".
  per_cohort$p_adj_BH <- NA_real_
  per_cohort$multiplicity_family <- NA_character_
  for (m in names(MODEL_PARAMS)) {
    ix <- which(per_cohort$model == m & per_cohort$meta_eligible & !is.na(per_cohort$p))
    if (!length(ix)) next
    per_cohort$p_adj_BH[ix] <- signif(stats::p.adjust(per_cohort$p[ix], method = "BH"), 4)
    per_cohort$multiplicity_family[ix] <-
      sprintf("BH within %s, k=%d meta-analysed cohorts", m, length(ix))
  }
  per_cohort$multiplicity_family[!per_cohort$meta_eligible] <-
    "not in family (descriptive cohort, B.n)"

  # BH is monotone and never shrinks a p-value; an adjusted value below its raw
  # value would mean the family was mis-assembled.
  bad_bh <- which(!is.na(per_cohort$p_adj_BH) & per_cohort$p_adj_BH < per_cohort$p - 1e-9)
  if (length(bad_bh))
    halt("B.n", "BH-adjusted p is below the raw p for: ",
         paste(per_cohort$cohort[bad_bh], per_cohort$model[bad_bh], collapse = ", "))
  if (any(!is.na(per_cohort$p_adj_BH[!per_cohort$meta_eligible])))
    halt("B.n", "a descriptive cohort received an adjusted p-value; ",
         "CHOL is not in the multiplicity family")
  n_fam <- vapply(names(MODEL_PARAMS), function(m)
    sum(per_cohort$model == m & !is.na(per_cohort$p_adj_BH)), integer(1))
  message("  ok  B.n        BH applied within model, family sizes: ",
          paste(sprintf("%s=%d", names(n_fam), n_fam), collapse = " "),
          "  (CHOL excluded; no correction across M1-M4)")
  n_sig_raw <- sum(per_cohort$p < 0.05 & per_cohort$meta_eligible, na.rm = TRUE)
  n_sig_adj <- sum(per_cohort$p_adj_BH < 0.05, na.rm = TRUE)
  message(sprintf("  ..  B.n        per-cohort tests at 0.05: %d raw, %d after BH",
                  n_sig_raw, n_sig_adj))

  write.csv(per_cohort, file.path(OUTDIR, "survival_per_cohort.csv"), row.names = FALSE)

  # ---- attenuation + paired bootstrap (B.k) --------------------------------
  message("  ..  B.k        paired bootstrap, B = ", B_RESAMPLES, ", seed ", SEED_BASE)
  att <- do.call(rbind, lapply(seq_along(COHORTS), function(i) {
    cc <- COHORTS[i]; r <- res[[cc]]
    b <- vapply(names(MODEL_PARAMS), function(m)
      if (isTRUE(r$fits[[m]]$ok)) r$fits[[m]]$beta else NA_real_, numeric(1))
    if (any(!is.finite(b[c("M2", "M3", "M4")]))) return(NULL)
    a  <- attenuation(b[["M2"]], b[["M3"]], b[["M4"]])
    bs <- bootstrap_cohort(r$data, r$use_sex, i)
    ci <- lapply(c("attenuation_total", "attenuation_purity", "attenuation_stroma",
                   "prop_attenuated"),
                 function(k) boot_ci(bs$reps[, k]))
    names(ci) <- c("total", "purity", "stroma", "prop")
    data.frame(cohort = cc, n = r$n, events = r$events,
               beta_M1 = round(b[["M1"]], 6), beta_M2 = round(b[["M2"]], 6),
               beta_M3 = round(b[["M3"]], 6), beta_M4 = round(b[["M4"]], 6),
               attenuation_total  = round(a[["attenuation_total"]], 6),
               att_total_lo = round(ci$total[1], 6), att_total_hi = round(ci$total[2], 6),
               attenuation_purity = round(a[["attenuation_purity"]], 6),
               att_purity_lo = round(ci$purity[1], 6), att_purity_hi = round(ci$purity[2], 6),
               attenuation_stroma = round(a[["attenuation_stroma"]], 6),
               att_stroma_lo = round(ci$stroma[1], 6), att_stroma_hi = round(ci$stroma[2], 6),
               # F10: the per-cohort CI is a PERCENTILE interval while the
               # meta-analysis treats each cohort as Normal(est, SD^2). Both are
               # reported so the divergence is visible rather than implicit.
               att_total_norm_lo = round(a[["attenuation_total"]] -
                 1.96 * stats::sd(bs$reps[, "attenuation_total"], na.rm = TRUE), 6),
               att_total_norm_hi = round(a[["attenuation_total"]] +
                 1.96 * stats::sd(bs$reps[, "attenuation_total"], na.rm = TRUE), 6),
               prop_attenuated = round(a[["prop_attenuated"]], 4),
               prop_lo = round(ci$prop[1], 4), prop_hi = round(ci$prop[2], 4),
               boot_ok = bs$n_ok, boot_failed = bs$n_failed,
               # resamples where a NUISANCE coefficient diverged/aliased; the
               # score coefficient is screened separately and independently
               boot_nuisance_unstable = bs$n_nuisance_unstable,
               # Bootstrap SE of EACH paired difference, for the meta-analysis.
               # attenuation_purity (b2-b3) and attenuation_stroma (b3-b4) are
               # different quantities from attenuation_total (b2-b4) and have
               # different variances; reusing the total's SE as the weight for all
               # three would give wrong tau2, CI, Q and I2 for the components.
               att_total_se  = round(stats::sd(bs$reps[, "attenuation_total"],  na.rm = TRUE), 6),
               att_purity_se = round(stats::sd(bs$reps[, "attenuation_purity"], na.rm = TRUE), 6),
               att_stroma_se = round(stats::sd(bs$reps[, "attenuation_stroma"], na.rm = TRUE), 6),
               meta_eligible = unname(META_ELIGIBLE[[cc]]),
               stringsAsFactors = FALSE)
  }))
  write.csv(att, file.path(OUTDIR, "attenuation_per_cohort.csv"), row.names = FALSE)

  # ---- meta-analysis (B.l) -------------------------------------------------
  meta_rows <- function(tag, sub, value_col, se_col) {
    m <- meta_one(sub[[value_col]], sub[[se_col]]^2, sub$cohort)
    if (is.null(m)) return(NULL)
    data.frame(analysis = tag, k = m$k,
               cohorts = paste(sub("^TCGA-", "", m$cohorts), collapse = "+"),
               n_total = sum(sub$n[match(m$cohorts, sub$cohort)]),
               est = round(m$est, 6), se_hksj = round(m$se_hksj, 6),
               ci_lo = round(m$ci_lo_hksj, 6), ci_hi = round(m$ci_hi_hksj, 6),
               p = signif(m$p_hksj, 4),
               HR = round(exp(m$est), 4),
               HR_lo = round(exp(m$ci_lo_hksj), 4), HR_hi = round(exp(m$ci_hi_hksj), 4),
               ci_lo_wald = round(m$ci_lo_wald, 6), ci_hi_wald = round(m$ci_hi_wald, 6),
               fe_est = round(m$fe_est, 6),
               fe_ci_lo = round(m$fe_ci_lo, 6), fe_ci_hi = round(m$fe_ci_hi, 6),
               tau2 = round(m$tau2, 6), I2 = round(m$I2, 2),
               Q = round(m$Q, 4), Q_df = m$Q_df, Q_p = signif(m$Q_p, 4),
               pi_lo = round(m$pi_lo, 6), pi_hi = round(m$pi_hi, 6),
               stringsAsFactors = FALSE)
  }

  pc_meta <- per_cohort[per_cohort$meta_eligible & per_cohort$fitted &
                        per_cohort$band != "do_not_fit", ]
  meta_models <- do.call(rbind, lapply(names(MODEL_PARAMS), function(m)
    meta_rows(m, pc_meta[pc_meta$model == m, ], "beta", "se")))

  att_meta_set <- att[att$meta_eligible & is.finite(att$att_total_se) & att$att_total_se > 0, ]
  meta_att <- meta_rows("attenuation_total", att_meta_set, "attenuation_total", "att_total_se")
  att_p_set  <- att[att$meta_eligible & is.finite(att$att_purity_se) & att$att_purity_se > 0, ]
  att_s_set  <- att[att$meta_eligible & is.finite(att$att_stroma_se) & att$att_stroma_se > 0, ]
  meta_att_p <- meta_rows("attenuation_purity", att_p_set, "attenuation_purity", "att_purity_se")
  meta_att_s <- meta_rows("attenuation_stroma", att_s_set, "attenuation_stroma", "att_stroma_se")

  # F8: a pool that fails to form must not vanish silently -- rbind drops NULL
  # without comment, so the PRIMARY RESULT row would simply be absent.
  if (is.null(meta_att))
    halt("B.l", "the primary attenuation pool did not form: fewer than 2 cohorts ",
         "with a finite positive bootstrap SE")
  message("  ok  B.l        primary attenuation pool formed over k = ", meta_att$k,
          " cohorts (", meta_att$cohorts, ")")

  # F7: the registered leave-one-out sensitivities apply to the PRIMARY ESTIMAND
  # too, not only to the four model pools. Amendment 8 requires the drop-all-low-
  # EPV pool; B.l requires leave-one-cohort-out as a table.
  low_epv_att <- unique(pc_meta$cohort[pc_meta$band == "fit_pool_flag_LOO"])
  att_sens <- list()
  if (length(low_epv_att)) {
    sset <- att_meta_set[!att_meta_set$cohort %in% low_epv_att, ]
    att_sens$drop_low_epv <- meta_rows("attenuation_total_drop_lowEPV", sset,
                                       "attenuation_total", "att_total_se")
  }
  att_sens$loco <- do.call(rbind, lapply(att_meta_set$cohort, function(cc)
    meta_rows(paste0("attenuation_total_LOCO_minus_", sub("^TCGA-", "", cc)),
              att_meta_set[att_meta_set$cohort != cc, ],
              "attenuation_total", "att_total_se")))

  # ---- contributing sets: does M2 differ from M4? (Amendment 8 items 3-4) ---
  set_M2 <- sort(pc_meta$cohort[pc_meta$model == "M2"])
  set_M4 <- sort(pc_meta$cohort[pc_meta$model == "M4"])
  sets_differ <- !identical(set_M2, set_M4)
  matched <- NULL
  if (sets_differ) {
    inter <- intersect(set_M2, set_M4)
    message("  !!  B.l        contributing sets DIFFER between M2 and M4; ",
            "matched re-pool over the intersection is REQUIRED (Amendment 8 item 4)")
    matched <- do.call(rbind, lapply(c("M2", "M4"), function(m) {
      s <- pc_meta[pc_meta$model == m & pc_meta$cohort %in% inter, ]
      meta_rows(paste0("matched_", m), s, "beta", "se")
    }))
  } else {
    message("  ok  B.l        contributing cohort set is identical for M2 and M4 (",
            length(set_M2), " cohorts); no matched re-pool needed")
  }

  # ---- registered sensitivities -------------------------------------------
  low_epv <- unique(pc_meta$cohort[pc_meta$band == "fit_pool_flag_LOO"])
  sens <- list()
  # (a) leave-one-out omitting ALL low-EPV cohorts simultaneously
  if (length(low_epv)) {
    s <- pc_meta[!pc_meta$cohort %in% low_epv, ]
    sens$drop_low_epv <- do.call(rbind, lapply(names(MODEL_PARAMS), function(m)
      meta_rows(paste0("drop_lowEPV_", m), s[s$model == m, ], "beta", "se")))
  }
  # (b) leave-one-cohort-out, per model (influence)
  sens$loco <- do.call(rbind, lapply(names(MODEL_PARAMS), function(m) {
    s <- pc_meta[pc_meta$model == m, ]
    do.call(rbind, lapply(s$cohort, function(cc)
      meta_rows(paste0("LOCO_", m, "_minus_", sub("^TCGA-", "", cc)),
                s[s$cohort != cc, ], "beta", "se")))
  }))

  write.csv(do.call(rbind, list(meta_models, meta_att, meta_att_p, meta_att_s, matched,
                                sens$drop_low_epv, att_sens$drop_low_epv)),
            file.path(OUTDIR, "meta_analysis.csv"), row.names = FALSE)
  write.csv(rbind(sens$loco, att_sens$loco),
            file.path(OUTDIR, "meta_leave_one_cohort_out.csv"), row.names = FALSE)

  # ---- PH, VIF, CHOL descriptive ------------------------------------------
  ph_all <- do.call(rbind, lapply(COHORTS, function(cc) {
    p <- res[[cc]]$ph; if (is.null(p)) return(NULL)
    cbind(cohort = cc, p)
  }))
  write.csv(ph_all, file.path(OUTDIR, "ph_tests.csv"), row.names = FALSE)

  viol <- unique(ph_all$cohort[ph_all$violated_score | ph_all$violated_global])
  ph_sens <- if (length(viol)) do.call(rbind, lapply(viol, function(cc) {
    r <- res[[cc]]
    ms <- ph_all$model[ph_all$cohort == cc &
                       (ph_all$violated_score | ph_all$violated_global)]
    do.call(rbind, lapply(ms, function(m) {
      x <- ph_sensitivity(r$data, m, r$use_sex)
      if (is.null(x)) NULL else cbind(cohort = cc, x)
    }))
  })) else NULL
  if (!is.null(ph_sens))
    write.csv(ph_sens, file.path(OUTDIR, "ph_sensitivity_score_logtime.csv"), row.names = FALSE)

  vifs <- do.call(rbind, lapply(COHORTS, function(cc) {
    f <- res[[cc]]$fits$M4
    if (!isTRUE(f$ok)) return(NULL)
    v <- vif_m4(f$fit); if (is.null(v)) return(NULL)
    cbind(cohort = cc, v)
  }))
  write.csv(vifs, file.path(OUTDIR, "vif_m4.csv"), row.names = FALSE)

  # correlation of the main score with BOTH stromal variants (the fallback
  # replaced one with the other, so both are reported)
  strom_r <- do.call(rbind, lapply(COHORTS, function(cc) {
    x <- d[d$cohort == cc, ]
    data.frame(cohort = cc,
               r_estimate_stromal = round(cor(x$score, x$stromal_score), 4),
               r_panel_subscore   = round(cor(x$score, x$stromal_score_subscore), 4),
               r_score_purity     = round(cor(x$score, x$purity), 4),
               stringsAsFactors = FALSE)
  }))
  write.csv(strom_r, file.path(OUTDIR, "stromal_correlations.csv"), row.names = FALSE)

  # ---- registered sensitivities requiring a full refit ---------------------
  # Each swaps ONE input and re-runs the identical path: complete-case set, four
  # models, EPV bands, meta-analysis. No bootstrap (the primary interval is the
  # primary analysis's); point estimates and model-based SEs.
  refit_variant <- function(tag, transform) {
    out <- do.call(rbind, lapply(seq_along(COHORTS), function(i) {
      cc <- COHORTS[i]
      dc <- transform(d[d$cohort == cc, , drop = FALSE])
      if (is.null(dc) || !nrow(dc)) return(NULL)
      use_sex <- isTRUE(unique(dc$use_sex))
      s <- complete_case_set(dc, use_sex)
      s$stage_group <- droplevels(s$stage_group)
      if (nrow(s) < 20L) return(NULL)
      ev <- sum(s$event == 1L)
      fits <- fit_cohort(s, use_sex)
      for (m in names(fits))
        if (identical(epv_band(ev, model_formula(m, use_sex)$params), "do_not_fit"))
          fits[[m]] <- list(ok = FALSE, reason = "EPV < 5: not fitted per Amendment 8")
      do.call(rbind, lapply(names(MODEL_PARAMS), function(m) {
        f <- fits[[m]]; pr <- model_formula(m, use_sex)$params
        data.frame(variant = tag, cohort = cc, model = m, n = nrow(s), events = ev,
                   epv = round(ev / pr, 2), band = epv_band(ev, pr),
                   fitted = isTRUE(f$ok),
                   beta = if (isTRUE(f$ok)) round(f$beta, 6) else NA_real_,
                   se   = if (isTRUE(f$ok)) round(f$se, 6) else NA_real_,
                   HR   = if (isTRUE(f$ok)) round(exp(f$beta), 4) else NA_real_,
                   meta_eligible = unname(META_ELIGIBLE[[cc]]),
                   stringsAsFactors = FALSE)
      }))
    }))
    if (is.null(out)) return(NULL)
    # F13: record how many rows each variant actually dropped, so a near-identity
    # refit (no_redacted drops ONE LIHC patient) is legible as such.
    out$rows_dropped_vs_primary <- NA_integer_
    for (cc in unique(out$cohort))
      out$rows_dropped_vs_primary[out$cohort == cc] <-
        sum(d$cohort == cc) - out$n[out$cohort == cc][1]
    att_v <- do.call(rbind, lapply(unique(out$cohort), function(cc) {
      b <- setNames(out$beta[out$cohort == cc], out$model[out$cohort == cc])
      if (any(!is.finite(b[c("M2", "M4")]))) return(NULL)
      data.frame(variant = tag, cohort = cc,
                 attenuation_total = round(b[["M2"]] - b[["M4"]], 6),
                 stringsAsFactors = FALSE)
    }))
    list(per_cohort = out, attenuation = att_v,
         meta = do.call(rbind, lapply(names(MODEL_PARAMS), function(m) {
           s <- out[out$meta_eligible & out$fitted & out$band != "do_not_fit" & out$model == m, ]
           meta_rows(paste0(tag, "_", m), s, "beta", "se")
         })))
  }

  variants <- list(
    # alternative endpoint: the non-primary one, with its registered censoring
    # 05 names these sens_time/sens_event (verified against the committed
    # analysis set, not assumed).
    alt_endpoint = function(x) {
      x$time  <- x$sens_time; x$event <- x$sens_event; x },
    # ESTIMATE purity everywhere, in place of CPE where CPE was primary
    estimate_purity = function(x) { x$purity <- x$purity_estimate; x },
    # the 143-gene sensitivity score
    score_143 = function(x) { x$score <- x$score_143; x },
    # redaction-excluded refit (registered decision B)
    # 05 carries the raw CDR `redaction` column; a redacted patient has a
    # non-empty value there. Registered decision B: primary RETAINS them, this
    # sensitivity excludes them.
    no_redacted = function(x) x[is.na(x$redaction) | !nzchar(trimws(x$redaction)), , drop = FALSE]
  )
  have_sens_cols <- all(c("sens_time", "sens_event") %in% names(d))
  if (!have_sens_cols)
    halt("B.l", "the alternative-endpoint sensitivity is REGISTERED and its inputs ",
         "(sens_time/sens_event) are absent from the analysis set; it must not be ",
         "silently skipped")
  if (!"redaction" %in% names(d))
    halt("B.l", "the redaction-excluded sensitivity is registered but the ",
         "`redaction` column is absent from the analysis set")
  sens_out <- lapply(names(variants), function(v) {
    if (v == "alt_endpoint" && !have_sens_cols) return(NULL)
    message("  ..  B.l        sensitivity refit: ", v)
    refit_variant(v, variants[[v]])
  })
  names(sens_out) <- names(variants)
  # Amendment 8 item 4 applies WITHIN each sensitivity too: if a variant's M2 and
  # M4 pools rest on different cohorts, an apparent attenuation could be an
  # artefact of the differing sets. Realised: alt_endpoint drops READ from M4
  # (PFI->OS gives 25 events, EPV 3.57 < 5), so its matched re-pool is required.
  sens_matched <- do.call(rbind, lapply(names(sens_out), function(v) {
    z <- sens_out[[v]]; if (is.null(z)) return(NULL)
    o <- z$per_cohort
    poolable <- o$meta_eligible & o$fitted & o$band != "do_not_fit"
    s2 <- sort(o$cohort[poolable & o$model == "M2"])
    s4 <- sort(o$cohort[poolable & o$model == "M4"])
    if (identical(s2, s4)) return(NULL)
    inter <- intersect(s2, s4)
    message("  !!  B.l        ", v, ": M2 and M4 pools differ (M2 k=", length(s2),
            ", M4 k=", length(s4), "); matched re-pool over k=", length(inter))
    do.call(rbind, lapply(c("M2", "M4"), function(m)
      meta_rows(paste0(v, "_matched_", m),
                o[poolable & o$model == m & o$cohort %in% inter, ], "beta", "se")))
  }))
  if (!is.null(sens_matched))
    write.csv(sens_matched, file.path(OUTDIR, "sensitivity_matched_repool.csv"),
              row.names = FALSE)

  sens_pc   <- do.call(rbind, lapply(sens_out, function(z) if (is.null(z)) NULL else z$per_cohort))
  sens_meta <- do.call(rbind, lapply(sens_out, function(z) if (is.null(z)) NULL else z$meta))
  sens_att  <- do.call(rbind, lapply(sens_out, function(z) if (is.null(z)) NULL else z$attenuation))
  if (!is.null(sens_pc))
    write.csv(sens_pc, file.path(OUTDIR, "sensitivity_per_cohort.csv"), row.names = FALSE)
  if (!is.null(sens_meta))
    write.csv(sens_meta, file.path(OUTDIR, "sensitivity_meta.csv"), row.names = FALSE)
  if (!is.null(sens_att))
    write.csv(sens_att, file.path(OUTDIR, "sensitivity_attenuation.csv"), row.names = FALSE)

  # ---- CHOL, descriptively (Amendment 8: reported, never pooled) -----------
  chol <- per_cohort[per_cohort$cohort == "TCGA-CHOL", ]
  write.csv(chol, file.path(OUTDIR, "chol_descriptive.csv"), row.names = FALSE)

  message("\n-- per cohort, per model --");        print(per_cohort, row.names = FALSE)
  message("\n-- attenuation --");                  print(att, row.names = FALSE)
  message("\n-- meta-analysis --");                print(meta_models, row.names = FALSE)
  message("\n-- attenuation, meta-analysed (PRIMARY) --"); print(meta_att, row.names = FALSE)
  message("\n-- PH tests --");                     print(ph_all, row.names = FALSE)
  message("\n-- VIF (M4) --");                     print(vifs, row.names = FALSE)
  message("\n-- score/stromal/purity correlations --"); print(strom_r, row.names = FALSE)
  message("\nB.j/B.k/B.l complete.")
}
