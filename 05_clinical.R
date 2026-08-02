#!/usr/bin/env Rscript
# 05_clinical.R -- Part B.i: clinical assembly. NO score, NO purity, NO model.
#
# Builds the per-cohort analysis set: TCGA-CDR endpoints joined to the expression
# cohorts' sample metadata, filtered to one primary-tumour aliquot per patient,
# with the B.j covariates and Amendment 7's endpoint designation.
#
# SCOPE BOUNDARY: this script must not compute purity, construct a score, or fit
# a model. Those are 06 and 07. The driver asserts it below rather than merely
# announcing it -- an earlier header claimed enforcement that did not exist.
#
# STRUCTURE. Everything below the divider is a PURE FUNCTION of its arguments --
# no globals, no file reads inside the transforms, nothing that depends on which
# cohort ran before. The driver at the bottom is the only part that touches disk,
# and it is guarded by `if (sys.nframe() == 0L)` so that `source()`ing this file
# from 06/07 (or from a validation-cohort driver) yields the functions WITHOUT
# running anything. Discovery and validation cohorts therefore run the identical
# code path rather than a copy-pasted variant.
#
# CONVENTION: every function that filters returns BOTH the filtered data and the
# count dropped, so §3.3's per-step reporting is a property of the pipeline rather
# than something reconstructed afterwards.

suppressPackageStartupMessages({
  library(readxl)
  library(SummarizedExperiment)
})

for (p in c("readxl", "SummarizedExperiment", "S4Vectors", "tools")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("HALT [B.i]: required package '%s' is not installed", p), call. = FALSE)
}

halt <- function(section, ...)
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)

assert_n <- function(observed, expected, section, what) {
  if (!identical(as.integer(observed), as.integer(expected)))
    halt(section, what, ": expected ", expected, ", observed ", observed)
  message(sprintf("  ok  %-14s %-52s = %d", section, what, observed))
  invisible(TRUE)
}

# ---------------------------------------------------------------- registered
# All seven cohorts. Registered in analysis_plan.md B.i; panel_definition.md names
# only six of them (PAAD is absent there), so the plan is the authority.
COHORTS <- c("TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-ESCA",
             "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL")

# AMENDMENT 7. Primary endpoint = OS wherever CDR Table 3 marks OS usable without
# caution; where OS is cautioned, the CDR-preferred alternative. READ is the only
# one of the seven whose OS mark carries an asterisk, so READ is the only PFI
# cohort. Event counts played NO part in this designation and must not be used to
# revisit it. The non-primary endpoint is the prespecified sensitivity everywhere.
ENDPOINT_PRIMARY <- c("TCGA-COAD" = "OS",  "TCGA-READ" = "PFI", "TCGA-STAD" = "OS",
                      "TCGA-ESCA" = "OS",  "TCGA-PAAD" = "OS",  "TCGA-LIHC" = "OS",
                      "TCGA-CHOL" = "OS")

# AMENDMENT 8. A cohort whose Table 3 explanation declares the sample size too
# small for ALL endpoints is descriptive-only regardless of EPV. CHOL is currently
# the only cohort meeting it. Implemented as a FLAG, never by dropping the cohort:
# CHOL is still assembled, still reported, and merely carries meta_eligible=FALSE.
META_ELIGIBLE <- setNames(COHORTS != "TCGA-CHOL", COHORTS)

# Administrative censoring, applied identically across cohorts (B.i).
CENSOR_DAYS <- c(OS = 3650L, PFI = 1825L)   # 10 years OS, 5 years PFI

# Registered CDR 2018 snapshot (analysis_plan.md B.i, "descriptive context only").
# These describe the WORKBOOK, not the merged analysis set -- the plan states
# explicitly that "realised event counts in the merged analysis set will differ".
# They are asserted against the CDR source to prove it is the registered release,
# and NEVER against the realised set. See `check_cdr_snapshot()`.
CDR_SNAPSHOT <- data.frame(
  cohort     = COHORTS,
  n          = c(459L, 170L, 443L, 185L, 185L, 377L,  45L),
  os_events  = c(102L,  26L, 172L,  77L, 100L, 132L,  22L),
  pfi_events = c(123L,  39L, 143L,  87L, 110L, 185L,  23L),
  stringsAsFactors = FALSE)

SEX_MIN_PER_LEVEL <- 10L   # B.j: sex dropped in any cohort with <10 of either sex

# Estimated parameters per model (B.j). stage_group's explicit `missing` level
# makes it TWO parameters, not one; this is Amendment 8's EPV denominator.
MODEL_PARAMS <- c(M1 = 1L, M2 = 5L, M3 = 6L, M4 = 7L)

# AMENDMENT 8 EPV bands, applied PER MODEL: >=10 fitted and pooled; 5-10 fitted,
# pooled, flagged with a leave-one-out sensitivity; <5 not fitted for that cohort
# and that model.
epv_band <- function(epv) {
  if (is.na(epv)) NA_character_
  else if (epv >= 10) "fit_and_pool"
  else if (epv >= 5)  "fit_pool_flag_LOO"
  else                "do_not_fit"
}

# Provisional disposition recorded in the plan, for CONTRAST with the realised
# one. Divergence is reported, never silently accepted -- and never used to
# revise the rule.
PROVISIONAL_DISPOSITION <- c("TCGA-COAD" = "fitted", "TCGA-STAD" = "fitted",
                             "TCGA-ESCA" = "fitted", "TCGA-PAAD" = "fitted",
                             "TCGA-LIHC" = "fitted", "TCGA-READ" = "low_EPV",
                             "TCGA-CHOL" = "cohort_excluded")

# REDACTION (registered decision B, 2026-08-02): redacted samples are RETAINED in
# the primary analysis and their count reported per cohort; a prespecified
# sensitivity refits every model excluding them. TCGA redaction flags identify
# real problems (sample swaps, contamination, pathology failures), so exclusion is
# defensible -- but redaction is not registered in the plan, and a
# post-registration exclusion affecting the analysis set must be visible rather
# than folded into the primary. Both directions are therefore reported.

# CDR missing-value sentinels. TCGA encodes absent values as bracketed strings
# rather than NA, so a naive read leaves "[Not Available]" as a stage LEVEL.
# A9: "[Completed]" is NOT a missing-value sentinel in TCGA -- it is a real
# status value (treatment/follow-up completion) and was wrongly listed here.
CDR_NA_STRINGS <- c("[Not Available]", "[Not Applicable]", "[Unknown]",
                    "[Not Evaluated]", "[Discrepancy]", "")

# ============================================================ PURE FUNCTIONS
# Each takes its inputs explicitly and returns a value. No globals are read
# except the registered constants above, which are literals, not state.

#' Normalise TCGA's bracketed missing-value sentinels to NA.
na_sentinel <- function(x, sentinels = CDR_NA_STRINGS) {
  x <- trimws(as.character(x))
  x[x %in% sentinels] <- NA_character_
  x
}

#' Collapse an AJCC stage string to the registered three-level factor.
#'
#' B.j: "AJCC stage collapsed to I/II vs III/IV ... patients with missing stage
#' are retained with an explicit `missing` level rather than dropped." Three
#' levels, two estimated parameters -- the EPV denominator in Amendment 8 depends
#' on this, so the level count is asserted by the caller.
#'
#' Roman numerals are parsed from the stem, so "Stage IIIA" -> III/IV and
#' "Stage IIB" -> I/II. "Stage X" (extent cannot be assessed) is MISSING, not a
#' stage: it is an explicit statement that stage is unknown.
collapse_stage <- function(stage_raw) {
  s <- toupper(na_sentinel(stage_raw))
  s <- gsub("^STAGE\\s+", "", trimws(s))
  s <- gsub("[^IVX]", "", s)            # drop A/B/C substage and any stray text
  out <- rep(NA_character_, length(s))
  out[s %in% c("I", "II")]   <- "I/II"
  out[s %in% c("III", "IV")] <- "III/IV"
  # "X" and anything unparsed stay NA and become the explicit `missing` level.
  factor(ifelse(is.na(out), "missing", out), levels = c("I/II", "III/IV", "missing"))
}

#' Patient barcode (first 12 characters) from a full TCGA aliquot barcode.
#'
#' TCGA-W5-AA39-01A-11R-A41I-07 -> TCGA-W5-AA39. Asserted rather than assumed:
#' a truncated or reformatted barcode would join silently to nothing.
patient_of <- function(barcode, section = "B.i") {
  b <- as.character(barcode)
  if (any(is.na(b))) halt(section, "NA barcode in sample metadata")
  bad <- !grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-", b)
  if (any(bad))
    halt(section, "barcode does not match the TCGA aliquot format: ",
         paste(utils::head(b[bad], 3), collapse = ", "))
  substr(b, 1, 12)
}

#' Apply administrative censoring to one endpoint.
#'
#' Beyond the cutoff the patient is known to be alive/progression-free AT the
#' cutoff, so the event is set to 0 and the time truncated. Returns the count
#' affected so the driver can report it -- a silent truncation would change event
#' counts with no trace.
censor_admin <- function(time, event, cutoff_days) {
  keep <- !is.na(time)
  beyond <- keep & time > cutoff_days
  n_events_lost <- sum(beyond & !is.na(event) & event == 1L)
  # A2: preserve NA. Conditioning `beyond` on time alone would coerce a MISSING
  # event indicator to 0 -- silently reclassifying an unknown outcome as
  # known-event-free at the cutoff instead of letting finalise_cohort drop it.
  event[beyond & !is.na(event)] <- 0L
  time[beyond]  <- cutoff_days
  list(time = time, event = event,
       n_truncated = sum(beyond), n_events_censored = n_events_lost)
}

#' Read the CDR workbook. Pure w.r.t. its path argument; the only I/O in the
#' function layer, isolated here so every transform below is testable offline.
read_cdr <- function(path, section = "B.i") {
  if (!file.exists(path)) halt(section, "CDR workbook not found: ", path)
  # Column types are pinned rather than guessed. readxl infers from the first
  # 1000 rows, which types `residual_tumor` (col 25) as logical and then emits
  # ~1,277 coercion warnings when it meets "R0"/"R1"/"R2" further down. That
  # column is unused here, but a wall of type warnings on a clinical workbook is
  # exactly the noise a real problem would hide in, so every column is read as
  # text and the numeric fields are converted explicitly below.
  d <- suppressMessages(readxl::read_excel(path, sheet = "TCGA-CDR", col_types = "text"))
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  num_cols <- c("age_at_initial_pathologic_diagnosis", "OS", "OS.time",
                "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")
  for (n in intersect(num_cols, names(d))) {
    v <- suppressWarnings(as.numeric(na_sentinel(d[[n]])))
    if (all(is.na(v)) && any(!is.na(d[[n]])))
      halt(section, "column '", n, "' did not parse as numeric; the workbook layout has changed")
    d[[n]] <- v
  }
  need <- c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis",
            "gender", "ajcc_pathologic_tumor_stage", "OS", "OS.time",
            "PFI", "PFI.time", "Redaction")
  miss <- setdiff(need, names(d))
  if (length(miss))
    halt(section, "CDR workbook is missing expected column(s): ",
         paste(miss, collapse = ", "), ". This is not the registered release.")
  d$md5 <- NULL
  d
}

#' Assert the CDR workbook IS the registered 2018 release.
#'
#' The registered per-cohort N and event counts (B.i) describe THIS WORKBOOK. The
#' plan states realised counts in the merged set will differ, so the snapshot is
#' checked here, at the source, and never against the merged analysis set. A
#' mismatch means the workbook is a different CDR release and every downstream
#' number would be incomparable to the registration.
check_cdr_snapshot <- function(cdr, snapshot = CDR_SNAPSHOT, section = "B.i") {
  # A3: without this, an NA in `type` makes every sum() NA, the comparison NA,
  # which() empty, and the check passes vacuously.
  if (anyNA(cdr$type)) halt(section, "CDR 'type' column contains NA")
  code <- sub("^TCGA-", "", snapshot$cohort)
  obs <- data.frame(
    cohort     = snapshot$cohort,
    n          = vapply(code, function(k) sum(cdr$type == k), integer(1)),
    os_events  = vapply(code, function(k) sum(cdr$type == k & cdr$OS  == 1L, na.rm = TRUE), integer(1)),
    pfi_events = vapply(code, function(k) sum(cdr$type == k & cdr$PFI == 1L, na.rm = TRUE), integer(1)),
    stringsAsFactors = FALSE)
  bad <- which(obs$n != snapshot$n | obs$os_events != snapshot$os_events |
               obs$pfi_events != snapshot$pfi_events)
  if (length(bad))
    halt(section, "the CDR workbook does not reproduce the registered 2018 snapshot ",
         "(analysis_plan.md B.i). ",
         paste(sprintf("%s: observed n=%d/OS=%d/PFI=%d, registered n=%d/OS=%d/PFI=%d",
                       obs$cohort[bad], obs$n[bad], obs$os_events[bad], obs$pfi_events[bad],
                       snapshot$n[bad], snapshot$os_events[bad], snapshot$pfi_events[bad]),
               collapse = "; "),
         ". This is a DIFFERENT CDR release; do not proceed -- every registered ",
         "count downstream would be incomparable.")
  message("  ok  B.i            CDR workbook reproduces the registered 2018 snapshot (7 cohorts)")
  invisible(obs)
}

#' Sample-level metadata for one cohort, as a plain data.frame.
#'
#' Reads colData ONLY -- never an assay. 05 assembles clinical data; expression is
#' 06's business. Kept separate so the filter chain below is a pure function of a
#' plain frame and can be unit-tested without a SummarizedExperiment.
sample_meta <- function(se, section = "B.i") {
  cd <- SummarizedExperiment::colData(se)
  need <- c("barcode", "sample_type")
  miss <- setdiff(need, colnames(cd))
  if (length(miss)) halt(section, "colData lacks: ", paste(miss, collapse = ", "))
  data.frame(barcode     = as.character(cd$barcode),
             sample_type = as.character(cd$sample_type),
             stringsAsFactors = FALSE)
}

#' Decompose a TCGA aliquot barcode into its ordered fields.
#'
#' TCGA-W5-AA39-01A-11R-A41I-07
#'  1-12 patient | 14-15 sample | 16 vial | 18-19 portion | 20 analyte
#'  22-25 plate  | 27-28 centre
#' Asserted rather than assumed: a reformatted barcode would silently yield NA
#' sort keys and make the tie-break arbitrary.
barcode_fields <- function(barcode, section = "B.i") {
  b <- as.character(barcode)
  if (any(nchar(b) < 25L))
    halt(section, "aliquot barcode too short to carry vial/portion/plate: ",
         paste(utils::head(b[nchar(b) < 25L], 3), collapse = ", "))
  f <- data.frame(
    patient = substr(b,  1, 12),
    sample  = substr(b, 14, 15),
    vial    = substr(b, 16, 16),
    portion = suppressWarnings(as.integer(substr(b, 18, 19))),
    plate   = substr(b, 22, 25),
    stringsAsFactors = FALSE)
  if (any(!grepl("^[A-Z]$", f$vial)))
    halt(section, "vial field is not a single letter: ",
         paste(utils::head(b[!grepl("^[A-Z]$", f$vial)], 3), collapse = ", "))
  if (any(is.na(f$portion)))
    halt(section, "portion field is not numeric: ",
         paste(utils::head(b[is.na(f$portion)], 3), collapse = ", "))
  f
}

#' Filter to one primary-tumour aliquot per patient, recording every drop.
#'
#' Steps, in order, each counted:
#'   1. primary tumour only  (sample_type == "Primary Tumor")
#'   2. one aliquot per patient
#'
#' TIE-BREAK (registered decision A, 2026-08-02): earliest VIAL letter, then
#' lowest PORTION, then lowest PLATE. Vial letter orders the successive vials cut
#' from a specimen and portion orders the pieces taken from a vial, so this
#' prefers the earliest-derived material. It deliberately does NOT sort the full
#' barcode: that string leads with the plate and centre codes, which order on
#' processing batch and sequencing site -- arbitrary with respect to sample
#' quality, and liable to prefer a late re-extraction over the original.
#' Applies only to patients with several sequenced aliquots; `n_multi_aliquot_
#' patients` and `aliquot_count_dist` report how many and how deep.
filter_samples <- function(meta, section = "B.i") {
  n0 <- nrow(meta)
  if (n0 == 0L) halt(section, "no samples in cohort metadata")

  is_primary <- meta$sample_type == "Primary Tumor"
  m1 <- meta[which(is_primary), , drop = FALSE]
  n_not_primary <- n0 - nrow(m1)
  if (nrow(m1) == 0L) halt(section, "no Primary Tumor samples in this cohort")

  f <- barcode_fields(m1$barcode, section)
  m1$patient <- f$patient
  tab <- table(m1$patient)
  n_multi <- sum(tab > 1L)
  dist <- table(as.integer(tab[tab > 1L]))   # how many patients had 2, 3, ... aliquots

  # A5: method="radix" forces C-locale collation. order()'s default is
  # locale-dependent for character vectors, so the vial and plate keys -- and
  # therefore which aliquot is retained -- could differ between machines.
  ord <- order(f$patient, f$vial, f$portion, f$plate, method = "radix")
  m1  <- m1[ord, , drop = FALSE]
  m2  <- m1[!duplicated(m1$patient), , drop = FALSE]
  n_dup_aliquot <- nrow(m1) - nrow(m2)

  if (anyDuplicated(m2$patient)) halt(section, "deduplication left duplicate patients")
  list(meta = m2,
       aliquot_count_dist = dist,
       counts = c(n_samples_in            = n0,
                  n_dropped_not_primary   = n_not_primary,
                  n_dropped_dup_aliquot   = n_dup_aliquot,
                  n_multi_aliquot_patients= n_multi,
                  n_after_sample_filter   = nrow(m2)))
}

#' Join CDR endpoints to the filtered samples, recording unmatched patients.
join_cdr <- function(meta, cdr, cohort, section = "B.i") {
  code <- sub("^TCGA-", "", cohort)
  c1 <- cdr[which(cdr$type == code), , drop = FALSE]
  if (nrow(c1) == 0L) halt(section, cohort, ": no CDR rows for type '", code, "'")
  if (anyDuplicated(c1$bcr_patient_barcode))
    halt(section, cohort, ": CDR has duplicate patient barcodes")

  i <- match(meta$patient, c1$bcr_patient_barcode)
  n_unmatched <- sum(is.na(i))
  keep <- which(!is.na(i))
  out <- data.frame(
    cohort      = cohort,
    patient     = meta$patient[keep],
    barcode     = meta$barcode[keep],
    age         = as.numeric(c1$age_at_initial_pathologic_diagnosis[i[keep]]),
    sex         = na_sentinel(c1$gender[i[keep]]),
    stage_group = collapse_stage(c1$ajcc_pathologic_tumor_stage[i[keep]]),
    OS          = as.integer(c1$OS[i[keep]]),
    OS.time     = as.numeric(c1$OS.time[i[keep]]),
    PFI         = as.integer(c1$PFI[i[keep]]),
    PFI.time    = as.numeric(c1$PFI.time[i[keep]]),
    redaction   = na_sentinel(c1$Redaction[i[keep]]),
    stringsAsFactors = FALSE)
  out$sex <- factor(out$sex)
  list(data = out, counts = c(n_dropped_no_cdr_match = n_unmatched,
                              n_after_cdr_join       = nrow(out)))
}

#' Attach the endpoint designation and apply administrative censoring.
#'
#' `time`/`event` carry the PRIMARY endpoint; `sens_time`/`sens_event` the
#' prespecified sensitivity (the other endpoint), so downstream code never has to
#' re-derive which is which -- the single place that decision is made is
#' ENDPOINT_PRIMARY.
apply_endpoints <- function(df, cohort, endpoint_map = ENDPOINT_PRIMARY,
                            censor = CENSOR_DAYS, section = "B.i") {
  primary <- endpoint_map[[cohort]]
  if (is.null(primary) || !primary %in% c("OS", "PFI"))
    halt(section, cohort, ": no registered primary endpoint (Amendment 7)")
  sens <- setdiff(c("OS", "PFI"), primary)

  pick <- function(ep) list(time = df[[paste0(ep, ".time")]], event = df[[ep]])
  p <- pick(primary); s <- pick(sens)
  pc <- censor_admin(p$time, p$event, censor[[primary]])
  sc <- censor_admin(s$time, s$event, censor[[sens]])

  df$endpoint_primary <- primary
  df$endpoint_sens    <- sens
  df$time             <- pc$time
  df$event            <- pc$event
  df$sens_time        <- sc$time
  df$sens_event       <- sc$event

  n_no_time <- sum(is.na(df$time) | is.na(df$event))
  list(data = df,
       counts = c(n_primary_truncated  = pc$n_truncated,
                  n_primary_events_censored = pc$n_events_censored,
                  n_sens_truncated     = sc$n_truncated,
                  n_missing_endpoint   = n_no_time))
}

#' Drop patients with no usable primary endpoint, and set the per-cohort sex rule.
#'
#' B.j: sex is "dropped in any cohort with < 10 patients of either sex". That is a
#' COHORT-level decision, so it is recorded as a flag on every row rather than by
#' deleting the column -- 07 reads the flag when it assembles the model formula.
finalise_cohort <- function(df, cohort, meta_eligible = META_ELIGIBLE,
                            sex_min = SEX_MIN_PER_LEVEL, section = "B.i") {
  n_in <- nrow(df)
  ok <- !is.na(df$time) & !is.na(df$event) & df$time >= 0
  out <- df[which(ok), , drop = FALSE]
  n_dropped <- n_in - nrow(out)
  if (nrow(out) == 0L) halt(section, cohort, ": no patients with a usable primary endpoint")

  sex_tab <- table(droplevels(out$sex))
  use_sex <- length(sex_tab) >= 2L && all(sex_tab >= sex_min)
  out$use_sex <- use_sex
  out$meta_eligible <- unname(meta_eligible[[cohort]])

  list(data = out,
       counts = c(n_dropped_no_endpoint = n_dropped, n_final = nrow(out)),
       sex_counts = sex_tab, use_sex = use_sex)
}

#' Assemble one cohort end to end. Pure: all inputs explicit.
build_cohort_clinical <- function(se, cdr, cohort) {
  sec <- paste0("B.i/", sub("^TCGA-", "", cohort))
  f <- filter_samples(sample_meta(se, sec), sec)
  j <- join_cdr(f$meta, cdr, cohort, sec)
  e <- apply_endpoints(j$data, cohort, section = sec)
  z <- finalise_cohort(e$data, cohort, section = sec)

  # --- invariants, each halting -------------------------------------------
  if (anyDuplicated(z$data$patient)) halt(sec, "duplicate patients in the final set")
  if (!all(levels(z$data$stage_group) == c("I/II", "III/IV", "missing")))
    halt(sec, "stage_group levels are not the registered three (I/II, III/IV, missing); ",
         "Amendment 8's EPV denominator assumes exactly two parameters for this term")
  if (nlevels(z$data$stage_group) != 3L)
    halt(sec, "stage_group must retain all three levels even when one is unobserved")
  cap <- CENSOR_DAYS[[z$data$endpoint_primary[1]]]
  if (any(z$data$time > cap)) halt(sec, "administrative censoring did not apply")
  if (!all(z$data$event %in% c(0L, 1L))) halt(sec, "primary event is not 0/1")
  if (!identical(unname(z$data$endpoint_primary[1]), unname(ENDPOINT_PRIMARY[[cohort]])))
    halt(sec, "endpoint designation departs from Amendment 7")

  # Realised EPV per model from the realised event count (Amendment 8), computed
  # here rather than in 08 so the disposition is known before any model is fitted.
  ev <- sum(z$data$event == 1L)
  epv <- ev / MODEL_PARAMS
  epv_tab <- data.frame(cohort = cohort, model = names(MODEL_PARAMS),
                        params = as.integer(MODEL_PARAMS),
                        events = ev, epv = round(unname(epv), 2),
                        band = vapply(unname(epv), epv_band, character(1)),
                        stringsAsFactors = FALSE)

  counts <- c(f$counts, j$counts, e$counts, z$counts)
  list(data = z$data, counts = counts, use_sex = z$use_sex,
       sex_counts = z$sex_counts, epv = epv_tab,
       aliquot_count_dist = f$aliquot_count_dist,
       summary = data.frame(
         cohort            = cohort,
         endpoint_primary  = z$data$endpoint_primary[1],
         endpoint_sens     = z$data$endpoint_sens[1],
         n_final           = nrow(z$data),
         events_primary    = sum(z$data$event == 1L),
         events_sens       = sum(z$data$sens_event == 1L, na.rm = TRUE),
         use_sex           = z$use_sex,
         meta_eligible     = unname(META_ELIGIBLE[[cohort]]),
         n_stage_missing   = sum(z$data$stage_group == "missing"),
         n_redacted        = sum(!is.na(z$data$redaction)),
         stringsAsFactors  = FALSE))
}

# ==================================================================== driver
# Guarded: `source("05_clinical.R")` from 06/07 or a validation driver defines the
# functions above and runs NOTHING.
if (sys.nframe() == 0L) {

  OUTDIR <- "output"; dir.create(OUTDIR, showWarnings = FALSE)
  message("\n== B.i  clinical assembly (no score, no purity, no model) ==")

  cdr <- read_cdr("data/manual/TCGA-CDR.xlsx")
  check_cdr_snapshot(cdr)
  assert_n(length(COHORTS), 7L, "B.i", "registered cohorts")
  if (!setequal(names(ENDPOINT_PRIMARY), COHORTS))
    halt("B.i", "the endpoint map does not cover exactly the registered cohorts")
  if (!identical(sort(names(ENDPOINT_PRIMARY)[ENDPOINT_PRIMARY == "PFI"]), "TCGA-READ"))
    halt("B.i", "Amendment 7 designates READ as the only PFI cohort; the map says otherwise")
  if (!identical(names(META_ELIGIBLE)[!META_ELIGIBLE], "TCGA-CHOL"))
    halt("B.i", "Amendment 8 makes CHOL the only descriptive-only cohort")

  res <- lapply(COHORTS, function(cc) {
    se <- readRDS(sprintf("data/tcga/%s_se.rds", cc))
    on.exit(rm(se))
    build_cohort_clinical(se, cdr, cc)
  })
  names(res) <- COHORTS

  clinical <- do.call(rbind, lapply(res, `[[`, "data"))
  summary_tab <- do.call(rbind, lapply(res, `[[`, "summary"))
  flow <- do.call(rbind, lapply(COHORTS, function(cc)
    data.frame(cohort = cc, step = names(res[[cc]]$counts),
               n = as.integer(res[[cc]]$counts), stringsAsFactors = FALSE)))

  if (anyDuplicated(clinical$patient))
    halt("B.i", "a patient appears in more than one cohort")
  assert_n(sum(summary_tab$meta_eligible), 6L, "B.i", "cohorts eligible for meta-analysis")
  assert_n(nrow(summary_tab), 7L, "B.i", "cohorts assembled (CHOL flagged, not dropped)")

  epv_tab <- do.call(rbind, lapply(res, `[[`, "epv"))

  # -- STOP CONDITION 2: realised n across the six meta-eligible cohorts --------
  n_meta <- sum(summary_tab$n_final[summary_tab$meta_eligible])
  message(sprintf("\n  ..  n across the six meta-eligible cohorts = %d (plan expects ~1,480-1,500 pre-purity)", n_meta))
  if (n_meta < 1300L)
    halt("B.i", "STOP CONDITION 2: realised n across the six meta-eligible cohorts is ",
         n_meta, ", below the 1,300 floor. The plan expects ~1,480-1,500 before ",
         "purity attrition. Report before proceeding; do not relax the filters.")

  # -- STOP CONDITION 3: redacted samples as a share of any cohort --------------
  red <- within(summary_tab, pct <- 100 * n_redacted / n_final)
  message("  ..  redacted: ",
          paste(sprintf("%s %d (%.1f%%)", sub("^TCGA-", "", red$cohort), red$n_redacted, red$pct),
                collapse = " | "))
  if (any(red$pct > 2))
    halt("B.i", "STOP CONDITION 3: redacted samples exceed 2% in ",
         paste(red$cohort[red$pct > 2], collapse = ", "),
         ". Report before proceeding.")

  # -- realised vs provisional disposition (reported, never used to revise) -----
  realised <- vapply(COHORTS, function(cc) {
    if (!META_ELIGIBLE[[cc]]) return("cohort_excluded")
    b <- epv_tab$band[epv_tab$cohort == cc]
    if (any(b == "do_not_fit")) "some_models_not_fitted"
    else if (any(b == "fit_pool_flag_LOO")) "low_EPV"
    else "fitted"
  }, character(1))
  disp <- data.frame(cohort = COHORTS,
                     provisional = unname(PROVISIONAL_DISPOSITION[COHORTS]),
                     realised = unname(realised),
                     differs = unname(realised) != unname(PROVISIONAL_DISPOSITION[COHORTS]),
                     stringsAsFactors = FALSE)

  write.csv(clinical,    file.path(OUTDIR, "clinical_analysis_set.csv"), row.names = FALSE)
  write.csv(summary_tab, file.path(OUTDIR, "clinical_cohort_summary.csv"), row.names = FALSE)
  write.csv(flow,        file.path(OUTDIR, "clinical_filter_flow.csv"), row.names = FALSE)
  write.csv(epv_tab,     file.path(OUTDIR, "clinical_epv.csv"), row.names = FALSE)
  write.csv(disp,        file.path(OUTDIR, "clinical_disposition.csv"), row.names = FALSE)

  message("\n-- per cohort --")
  print(summary_tab, row.names = FALSE)
  message("\n-- realised EPV (Amendment 8 bands) --")
  print(epv_tab, row.names = FALSE)
  message("\n-- disposition vs the provisional table --")
  print(disp, row.names = FALSE)
  if (any(disp$differs))
    message("  !!  disposition DIFFERS from provisional for: ",
            paste(disp$cohort[disp$differs], collapse = ", "),
            " -- reported per instruction; the Amendment 8 rule is applied as registered.")
  message("\n-- aliquot multiplicity among multi-aliquot patients --")
  for (cc in COHORTS) {
    d <- res[[cc]]$aliquot_count_dist
    message(sprintf("  %-10s %s", sub("^TCGA-", "", cc),
      if (!length(d)) "none" else paste(sprintf("%s aliquots: %d patients", names(d), as.integer(d)), collapse=" | ")))
  }
  # A7: assert the scope boundary rather than announcing it. If a model had been
  # fitted, survival would be loaded; if a score or purity had been built, the
  # analysis set would carry the column.
  if ("survival" %in% loadedNamespaces())
    halt("B.i", "the survival package is loaded; 05 must not fit a model")
  if (any(c("score", "purity", "stromal_score") %in% names(clinical)))
    halt("B.i", "the analysis set carries a score or purity column; 05 must not build one")
  message("\nB.i complete. HARD STOP asserted: survival not loaded, no score or ",
          "purity column present.")
}
