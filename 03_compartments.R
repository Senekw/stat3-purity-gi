# 03_compartments.R ============================================================
# Part A (A.a-A.g) of analysis_plan.md v1.5, under panel_definition.md
# Amendments 1-10.  Written from the plan; no code reused from the unattributed
# file removed in 12e058e (see NOTES_FOR_REVIEW.md section 10).
#
# Fail-closed by construction.  Every count, label and invariant the plan names
# is asserted BEFORE the quantity that depends on it is computed.  A halt is a
# result, not a failure: it means the data are not what the registration says
# they are, and the correct response is an amendment, never a loosened check.
#
# Scope: A.a through A.g ONLY.  No score construction, no survival model, no
# outcome data.  Part B is not called from here.
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})
options(stringsAsFactors = FALSE, warn = 1)

for (p in c("Matrix", "rhdf5", "withr")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("HALT [deps]: package '", p, "' is required and not installed.", call. = FALSE)
}

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}

# assert_n: the plan's named assertion helper.  Reports BOTH numbers, always --
# the removed file's halt message quoted expected values it had not computed.
assert_n <- function(observed, expected, section, what) {
  if (!isTRUE(observed == expected))
    halt(section, what, ": expected ", expected, ", observed ", observed,
         ". Do not adjust this assertion; amend the plan or fix the input.")
  message(sprintf("  ok  %-28s %-46s = %s", section, what, observed))
  invisible(TRUE)
}

ROOT   <- normalizePath(".", mustWork = TRUE)
RAW    <- file.path(ROOT, "data", "raw")
GEO    <- file.path(ROOT, "data", "geo")
OUTDIR <- file.path(ROOT, "output")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

need_file <- function(path, section) {
  if (!file.exists(path)) halt(section, "required file not found: ", path)
  path
}

# ==============================================================================
# A.a  ATLAS SET  (Amendments 9 and 10)
# ==============================================================================
# THREE atlases.  GSE155698 (Amdt 9) and GSE183904 (Amdt 10) were removed after
# direct inspection showed neither deposits any cell-type annotation.
#
# `tissue` is NOT taken on trust.  The removed file hardcoded a tissue vector
# that was transposed -- GSE125449 tagged gastric, GSE183904 tagged
# liver_biliary -- and nothing in it could detect the error.  Here each atlas
# declares a `tissue_witness`: a marker gene set and a diagnostic label pattern
# that must both be present in that atlas's OWN content.  The tissue claim is
# verified against the data in verify_tissue() below, and a mismatch halts.

atlases <- data.frame(
  atlas  = c("GSE125449",     "GSE178341",  "Peng"),
  tissue = c("liver_biliary", "colorectal", "pancreatic"),
  source = c("Ma et al. 2019 (GEO)", "Pelka et al. 2021 (GEO)",
             "Peng et al. 2019, Besca reprocessing (Zenodo 10.5281/zenodo.3969339)"),
  stringsAsFactors = FALSE
)

assert_n(nrow(atlases), 3L, "A.a", "atlases in the registered set")
if ("GSE155698" %in% atlases$atlas) halt("A.a", "GSE155698 must not be in the set (Amendment 9)")
if ("GSE183904" %in% atlases$atlas) halt("A.a", "GSE183904 must not be in the set (Amendment 10)")
if (anyDuplicated(atlases$atlas))   halt("A.a", "duplicate atlas identifiers")
if (anyDuplicated(atlases$tissue))
  halt("A.a", "duplicate tissue labels: Amendment 10 states one atlas per tissue")

# Tissue witnesses.  Marker genes must be present in the atlas's gene universe
# AND detected; the label pattern must match at least one of its own cell-type
# labels.  Deliberately organ-specific: a transposition of two rows in the frame
# above cannot satisfy the other atlas's witness.
TISSUE_WITNESS <- list(
  liver_biliary = list(
    genes  = c("ALB", "APOA1", "KRT19", "SERPINA1"),   # hepatocyte + cholangiocyte
    labels = "HPC-like|Malignant cell",
    note   = "hepatic/biliary parenchymal programme"
  ),
  colorectal = list(
    genes  = c("CDX2", "CDX1", "LGALS4", "TFF3"),      # intestinal identity
    labels = "^Epi$",
    note   = "intestinal CDX2-positive epithelium"
  ),
  pancreatic = list(
    genes  = c("PRSS1", "CPA1", "CTRB1", "KRT19"),     # acinar + ductal
    labels = "pancreatic|acinar|ductal|stellate",
    note   = "exocrine pancreatic programme"
  )
)
if (!setequal(names(TISSUE_WITNESS), atlases$tissue))
  halt("A.a", "tissue witnesses do not cover exactly the declared tissues")

# verify_tissue: the claim in `atlases$tissue` is checked against the atlas's own
# gene universe and label vocabulary.  Called once per atlas on open.
verify_tissue <- function(atlas_id, claimed_tissue, gene_universe, labels) {
  w <- TISSUE_WITNESS[[claimed_tissue]]
  if (is.null(w)) halt("A.a", atlas_id, ": no tissue witness for '", claimed_tissue, "'")

  present <- w$genes[w$genes %in% gene_universe]
  if (length(present) < 2L)
    halt("A.a", atlas_id, " claims tissue '", claimed_tissue, "' but only ",
         length(present), " of ", length(w$genes), " witness genes (",
         paste(w$genes, collapse = ", "), ") are in its gene universe. ",
         "Check the tissue assignment in `atlases` against the atlas source.")

  if (!any(grepl(w$labels, labels, ignore.case = TRUE)))
    halt("A.a", atlas_id, " claims tissue '", claimed_tissue,
         "' but no cell-type label matches the expected ", w$note,
         " (pattern: ", w$labels, "). Observed labels: ",
         paste(sort(unique(labels)), collapse = ", "))

  # Cross-check: the claimed tissue must fit BETTER than any other declared
  # tissue.  Catches a transposition that happens to share a marker.
  score <- vapply(names(TISSUE_WITNESS), function(tt) {
    ww <- TISSUE_WITNESS[[tt]]
    sum(ww$genes %in% gene_universe) + 10L * any(grepl(ww$labels, labels, ignore.case = TRUE))
  }, integer(1))
  if (score[[claimed_tissue]] < max(score))
    halt("A.a", atlas_id, " claims tissue '", claimed_tissue,
         "' but its content fits '", names(which.max(score)),
         "' better (scores: ", paste(names(score), score, sep = "=", collapse = ", "),
         "). Suspect a transposed tissue label.")

  message(sprintf("  ok  A.a  %-12s tissue '%s' verified against atlas content",
                  atlas_id, claimed_tissue))
  invisible(TRUE)
}

# ==============================================================================
# LOCKED PANEL  (panel_definition.md, locked 2026-07-31 at ac9c5e0)
# ==============================================================================
panel_file <- need_file(file.path(ROOT, "data", "panel", "panel_locked.csv"), "panel")
panel      <- read.csv(panel_file, check.names = FALSE)
if (!"gene" %in% names(panel)) halt("panel", "panel_locked.csv has no `gene` column")
assert_n(nrow(panel), 152L, "panel", "locked panel genes")

PANEL_GENES <- trimws(as.character(panel$gene))
if (any(!nzchar(PANEL_GENES))) halt("panel", "blank gene symbol in the locked panel")
if (anyDuplicated(PANEL_GENES)) halt("panel", "duplicate gene symbols in the locked panel")
assert_n(length(PANEL_GENES), 152L, "panel", "unique locked gene symbols")

# The origin six, reported as a labelled subset per section 4 of the
# prespecification.  BCL2, MMP9 and HGF are non-qualifying under Amendment 2
# (human-only ChIP-seq) and are reported separately, never silently dropped.
ORIGIN_SIX      <- c("SOCS3", "BCL2", "MYC", "MMP9", "HGF", "IL6")
ORIGIN_IN_PANEL <- intersect(ORIGIN_SIX, PANEL_GENES)          # expect SOCS3, MYC, IL6
ORIGIN_NONQUAL  <- setdiff(ORIGIN_SIX, PANEL_GENES)            # expect BCL2, MMP9, HGF
if (!setequal(ORIGIN_IN_PANEL, c("SOCS3", "MYC", "IL6")))
  halt("panel", "origin-six membership changed: in-panel = ",
       paste(ORIGIN_IN_PANEL, collapse = ", "), " (expected SOCS3, MYC, IL6)")
if (!setequal(ORIGIN_NONQUAL, c("BCL2", "MMP9", "HGF")))
  halt("panel", "origin-six non-qualifying set changed: ",
       paste(ORIGIN_NONQUAL, collapse = ", "), " (expected BCL2, MMP9, HGF)")

# ==============================================================================
# A.g INVARIANT, VERIFIED BEFORE ANY OBSERVED DATA IS USED
# ==============================================================================
# Three atlases => 3^3 = 27 evidence patterns over {dominant, not, NA}.
# Asserts k_all3 <= k_evalall <= k, and that the three do not collapse into
# each other (a sensitivity that can never differ is not a sensitivity).

verify_k_ordering <- function(n_atlas = 3L) {
  grid <- expand.grid(rep(list(c(TRUE, FALSE, NA)), n_atlas), KEEP.OUT.ATTRS = FALSE)
  assert_n(nrow(grid), 3L^n_atlas, "A.g", "evidence patterns enumerated")
  seen <- character(0)
  for (i in seq_len(nrow(grid))) {
    d       <- as.logical(unlist(grid[i, ]))
    n_dom   <- sum(d, na.rm = TRUE)
    n_eval  <- sum(!is.na(d))
    in_k    <- n_dom >= 2L
    in_ev   <- n_dom >= 2L && n_eval == n_atlas
    in_all  <- n_dom == n_atlas && n_eval == n_atlas
    if (!(in_all <= in_ev && in_ev <= in_k))
      halt("A.g", "k ordering violated at pattern ",
           paste(ifelse(is.na(d), "NA", d), collapse = ","),
           ": k_all3=", in_all, " k_evalall=", in_ev, " k=", in_k)
    seen <- c(seen, paste(in_k, in_ev, in_all))
  }
  n_distinct <- length(unique(seen))
  if (n_distinct < 3L)
    halt("A.g", "k, k_evalall and k_all3 collapse: only ", n_distinct,
         " distinct outcome triples over ", nrow(grid), " patterns. ",
         "A sensitivity that cannot differ from the primary is not informative.")
  message(sprintf("  ok  A.g  ordering k_all3 <= k_evalall <= k verified over %d patterns (%d distinct triples)",
                  nrow(grid), n_distinct))
  invisible(TRUE)
}
verify_k_ordering(3L)

# ==============================================================================
# A.c  COMPARTMENT HARMONISATION MAP  (Amendment 4)
# ==============================================================================
# EXACT MATCH ON ENUMERATED LABELS -- deliberately not regex.
#
# The removed file mapped by unanchored case-insensitive substring in six
# sequential passes, later passes overwriting earlier ones.  Because "T cell"
# matches inside "MalignanT CELL" and the lymphoid pass ran after the epithelial
# pass, GSE125449's `Malignant cell` mapped to LYMPHOID: every malignant cell in
# the biliary atlas counted as a lymphocyte, inverting the study's own estimand,
# and no assertion could catch it because `lymphoid` is a valid compartment.
#
# Here every source label is written out and looked up by exact string match.
# A label not in the table halts.  There is no fallback and no default.

COMPARTMENTS <- c("epithelial", "fibroblast_stromal", "myeloid",
                  "lymphoid", "endothelial", "other")

# --- GSE125449 (Ma): complete and verified against the committed per-cell table
MAP_GSE125449 <- c(
  "Malignant cell" = "epithelial",          # judgement: malignant hepatocyte/cholangiocyte
  "HPC-like"       = "epithelial",          # hepatic progenitor-like
  "CAF"            = "fibroblast_stromal",
  "TEC"            = "endothelial",         # NB: fibroblast_stromal in the pilot's 5-compartment scheme
  "TAM"            = "myeloid",
  "T cell"         = "lymphoid",
  "B cell"         = "lymphoid",
  "unclassified"   = "other"
)

# --- GSE178341 (Pelka): clTopLevel, read from the cluster file
MAP_GSE178341 <- c(
  "Epi"    = "epithelial",
  "Strom"  = "fibroblast_stromal",          # includes fibroblast/pericyte/endothelial subsets at cl295 level
  "Myeloid"= "myeloid",
  "Mast"   = "myeloid",                     # judgement: granulocyte lineage, per A.c rule
  "TNKILC" = "lymphoid",
  "Plasma" = "lymphoid",
  "B"      = "lymphoid"
)

# --- Peng (Besca celltype0): lineage level, read from obs/celltype0
MAP_Peng <- c(
  "epithelial cell"    = "epithelial",      # ductal + acinar + endocrine, per A.c judgement call
  "fibroblast"         = "fibroblast_stromal",
  "endothelial cell"   = "endothelial",
  "hematopoietic cell" = NA_character_,     # RESOLVED BELOW -- see note
  "neural cell"        = "other"            # neural crest lineage, not stromal (A.c)
)

COMPARTMENT_MAPS <- list(GSE125449 = MAP_GSE125449,
                         GSE178341 = MAP_GSE178341,
                         Peng      = MAP_Peng)

# `hematopoietic cell` at celltype0 collapses myeloid AND lymphoid into one
# level, so it cannot be mapped at this granularity without discarding the
# myeloid/lymphoid distinction that Amendment 4 requires as separate
# compartments.  Peng therefore uses celltype1 for the haematopoietic branch
# only; celltype0 supplies the rest.  This is a documented refinement of A.c's
# "level-1 label" wording, not a change of estimand: the epithelial /
# non-epithelial split -- the only split A.e's f(pi) depends on -- is identical
# either way.
MAP_Peng_ct1 <- c(
  "myeloid leukocyte"       = "myeloid",
  "T cell"                  = "lymphoid",
  "lymphocyte of B lineage" = "lymphoid"
)

# map_labels: exact lookup, halt on anything unmapped.  No partial matching, no
# regex, no silent NA.
map_labels <- function(labels, map, atlas_id) {
  lab <- as.character(labels)
  unmapped <- setdiff(unique(lab), names(map))
  if (length(unmapped))
    halt("A.c", atlas_id, " unmapped level-1 labels: ",
         paste(sort(unmapped), collapse = ", "),
         ". Extend the map in analysis_plan.md as a dated amendment, then re-run. ",
         "Do not add a catch-all.")
  out <- unname(map[lab])
  if (anyNA(out))
    halt("A.c", atlas_id, " produced NA compartments after mapping; ",
         "every enumerated label must resolve to one of: ",
         paste(COMPARTMENTS, collapse = ", "))
  bad <- setdiff(unique(out), COMPARTMENTS)
  if (length(bad))
    halt("A.c", atlas_id, " mapped to unknown compartment(s): ", paste(bad, collapse = ", "))
  out
}

# Regression guard for the exact defect that removed the previous file.
# Runs before any atlas is opened; failure here means the map has drifted back
# toward substring semantics.
local({
  chk <- c("Malignant cell" = "epithelial", "HPC-like" = "epithelial",
           "T cell" = "lymphoid", "TEC" = "endothelial", "CAF" = "fibroblast_stromal")
  got <- MAP_GSE125449[names(chk)]
  if (!identical(unname(got), unname(chk)))
    halt("A.c", "GSE125449 map regression: 'Malignant cell' must map to epithelial, ",
         "got ", paste(names(chk), unname(got), sep = "->", collapse = ", "))
  message("  ok  A.c  map regression guard passed ('Malignant cell' -> epithelial)")
})

# ==============================================================================
# A.a / A.b  ATLAS LOADERS
# ==============================================================================
# Each returns list(X = genes x cells raw UMI counts, cells = data.frame(
#   cell_id, patient_id, compartment), atlas, tissue).
# RAW counts only: no normalisation anywhere in this file (A.d).

# ------------------------------------------------------------------ GSE125449
load_GSE125449 <- function() {
  sec <- "A.a/GSE125449"
  dir <- file.path(GEO, "GSE125449")
  sets <- c("Set1", "Set2")

  meta <- do.call(rbind, lapply(sets, function(s) {
    f <- need_file(file.path(dir, sprintf("GSE125449_%s_samples.txt.gz", s)), sec)
    d <- read.delim(gzfile(f), check.names = FALSE)
    if (!all(c("Sample", "Cell Barcode", "Type") %in% names(d)))
      halt(sec, s, " schema changed; expected Sample / Cell Barcode / Type")
    d$set <- s
    d
  }))

  # All 19 samples are tumours (10 iCCA, 9 HCC): no tumour/normal filter needed.
  # Recorded as an explicit no-op with the premise asserted, per A.a.
  assert_n(length(unique(meta$Sample)), 19L, sec, "samples (all tumour, no filter)")

  mats <- lapply(sets, function(s) {
    m <- Matrix::readMM(gzfile(need_file(file.path(dir, sprintf("GSE125449_%s_matrix.mtx.gz", s)), sec)))
    g <- read.delim(gzfile(need_file(file.path(dir, sprintf("GSE125449_%s_genes.tsv.gz", s)), sec)),
                    header = FALSE)
    b <- read.delim(gzfile(need_file(file.path(dir, sprintf("GSE125449_%s_barcodes.tsv.gz", s)), sec)),
                    header = FALSE)
    rownames(m) <- g[[ncol(g)]]      # symbol column
    colnames(m) <- b[[1]]
    m
  })
  names(mats) <- sets

  # Genes must agree across sets before cbind, else the row space is incoherent.
  if (!identical(rownames(mats$Set1), rownames(mats$Set2)))
    halt(sec, "Set1 and Set2 gene rows differ; cannot combine without a join")

  X <- cbind(mats$Set1, mats$Set2)
  meta$cell_id <- meta$`Cell Barcode`
  keep <- meta$cell_id %in% colnames(X)
  if (!any(keep)) halt(sec, "no metadata barcode matches a matrix column")
  meta <- meta[keep, , drop = FALSE]
  X    <- X[, meta$cell_id, drop = FALSE]

  verify_tissue("GSE125449", "liver_biliary", rownames(X), meta$Type)

  cells <- data.frame(
    cell_id     = meta$cell_id,
    patient_id  = meta$Sample,                       # sample == patient here
    compartment = map_labels(meta$Type, MAP_GSE125449, "GSE125449"),
    stringsAsFactors = FALSE
  )
  assert_n(length(unique(cells$patient_id)), 19L, sec, "bootstrap units (patients)")
  list(X = X, cells = cells, atlas = "GSE125449", tissue = "liver_biliary")
}

# ------------------------------------------------------------------ GSE178341
# Three counts asserted (Amendment 10), because they count different things:
#   129 GEO GSM channels with specimen_type == "T"
#   128 channels carrying >= 1 cell in the deposited matrices
#    62 unique patients -- the bootstrap unit
# Plus the IDENTITY of the single absent channel, so that a *different* channel
# dropping in a future release halts instead of passing on an unchanged total.
#
# PID is the patient identifier.  PatientTypeID is NOT: it is patient x specimen
# and yields 64 over tumour cells, because C130 and C171 each contributed two
# spatially distinct tumour specimens (_TA/_TB).  Using it as the bootstrap unit
# would resample 64 pseudo-patients as independent and understate the interval.

GSE178341_ABSENT_CHANNEL <- "C144_T_1_1_12_c1_v2"   # GSM5388094, CD45pCD3nCD19nMACS

load_GSE178341 <- function() {
  sec <- "A.a/GSE178341"
  f_meta <- need_file(file.path(RAW, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"), sec)
  f_clus <- need_file(file.path(RAW, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"), sec)
  f_h5   <- need_file(file.path(RAW, "GSE178341_crc10x_full_c295v4_submit.h5"), sec)

  meta <- read.csv(gzfile(f_meta), check.names = FALSE)
  clus <- read.csv(gzfile(f_clus), check.names = FALSE)

  for (nm in c("cellID", "SPECIMEN_TYPE", "PID")) {
    if (!nm %in% names(meta))
      halt(sec, "per-cell metatable has no `", nm, "` column. ",
           if (nm == "PID")
             "PatientTypeID is NOT a substitute: it is patient x specimen and yields 64, not 62."
           else "Join via GSM before filtering; do not proceed.")
  }
  for (nm in c("sampleID", "batchID", "clTopLevel")) {
    if (!nm %in% names(clus)) halt(sec, "cluster file has no `", nm, "` column")
  }
  assert_n(nrow(meta), 370115L, sec, "cells in metatable")
  assert_n(nrow(clus), 370115L, sec, "cells in cluster file")
  if (anyDuplicated(meta$cellID)) halt(sec, "cellID is not unique in the metatable")
  if (!setequal(meta$cellID, clus$sampleID))
    halt(sec, "metatable cellID and cluster sampleID are not the same cell set")

  keep_cell <- trimws(meta$SPECIMEN_TYPE) == "T"
  if (!any(keep_cell)) halt(sec, "no tumour cells after SPECIMEN_TYPE filter")

  # (1) GEO-side channel count, parsed from the series metadata cached in A.b.
  f_geo <- file.path(RAW, "GSE178341_geo_tumour_channels.txt")
  if (!file.exists(f_geo))
    halt(sec, "missing ", basename(f_geo), " -- the GEO tumour channel list. ",
         "A.b requires the 129-channel GEO count to be asserted against a ",
         "committed artefact, not re-fetched at analysis time.")
  geo_channels <- unique(trimws(readLines(f_geo)))
  geo_channels <- geo_channels[nzchar(geo_channels)]
  assert_n(length(geo_channels), 129L, sec, "GEO GSM channels, specimen_type == T")

  # (2) channels actually carrying cells, via batchID -- NOT via row counts
  b_of_cell <- clus$batchID[match(meta$cellID[keep_cell], clus$sampleID)]
  in_data_channels <- sort(unique(b_of_cell[!is.na(b_of_cell)]))
  assert_n(length(in_data_channels), 128L, sec, "channels with >=1 cell, tumour")

  # (3) patients -- the bootstrap unit
  pid <- meta$PID[keep_cell]
  assert_n(length(unique(pid)), 62L, sec, "unique patients (PID), tumour")

  # PatientTypeID must NOT be used; assert the 64 explicitly so the distinction
  # is visible in the log rather than rediscovered as a bug.
  if ("PatientTypeID" %in% names(meta)) {
    n_ptid <- length(unique(meta$PatientTypeID[keep_cell]))
    if (n_ptid != 64L)
      halt(sec, "PatientTypeID yields ", n_ptid, " over tumour cells, expected 64 ",
           "(60 single-specimen patients + C130/C171 with _TA and _TB each). ",
           "The specimen structure has changed; re-check which patients are split.")
    message("  ok  A.a/GSE178341   PatientTypeID = 64 (patient x specimen) confirmed distinct from PID = 62")
  }

  # identity of the absent channel, not merely the count
  missing_ch <- setdiff(geo_channels, in_data_channels)
  if (!identical(sort(missing_ch), GSE178341_ABSENT_CHANNEL))
    halt(sec, "expected exactly ", GSE178341_ABSENT_CHANNEL,
         " to be absent from the deposited matrices; observed absent: ",
         if (length(missing_ch)) paste(sort(missing_ch), collapse = ", ") else "(none)")
  message("  ok  A.a/GSE178341   absent channel is ", GSE178341_ABSENT_CHANNEL, " as registered")

  lab <- clus$clTopLevel[match(meta$cellID[keep_cell], clus$sampleID)]
  X   <- read_h5_counts(f_h5, sec)
  X   <- X[, meta$cellID[keep_cell], drop = FALSE]

  verify_tissue("GSE178341", "colorectal", rownames(X), lab)

  cells <- data.frame(
    cell_id     = meta$cellID[keep_cell],
    patient_id  = pid,                                   # PID, never PatientTypeID
    compartment = map_labels(lab, MAP_GSE178341, "GSE178341"),
    stringsAsFactors = FALSE
  )
  list(X = X, cells = cells, atlas = "GSE178341", tissue = "colorectal")
}

# ------------------------------------------------------------------ helpers
# 10x-style HDF5 -> dgCMatrix, raw counts.
read_h5_counts <- function(path, sec) {
  ls_h5 <- rhdf5::h5ls(path, recursive = TRUE)
  grp   <- if ("matrix" %in% ls_h5$name) "matrix" else ls_h5$group[[2]]
  gx    <- function(n) rhdf5::h5read(path, paste0(grp, "/", n))
  need  <- c("data", "indices", "indptr", "shape", "barcodes")
  have  <- ls_h5$name[ls_h5$group == paste0("/", grp)]
  if (!all(need %in% have))
    halt(sec, "HDF5 group '", grp, "' lacks CSC fields; found: ", paste(have, collapse = ", "))
  shp <- as.integer(gx("shape"))
  fg  <- ls_h5[ls_h5$group == paste0("/", grp, "/features"), "name"]
  gsym <- if (length(fg)) as.character(gx("features/name")) else as.character(gx("genes"))
  M <- Matrix::sparseMatrix(i = as.integer(gx("indices")) + 1L,
                            p = as.integer(gx("indptr")),
                            x = as.numeric(gx("data")),
                            dims = shp, index1 = FALSE)
  rownames(M) <- gsym
  colnames(M) <- as.character(gx("barcodes"))
  M
}

# h5ad categorical -> character, for the Besca-style layout used by Peng.
h5ad_factor <- function(path, col) {
  codes <- rhdf5::h5read(path, paste0("obs/", col))
  cats  <- try(rhdf5::h5read(path, paste0("obs/__categories/", col)), silent = TRUE)
  if (inherits(cats, "try-error")) return(as.character(codes))
  as.character(cats)[as.integer(codes) + 1L]
}

# ----------------------------------------------------------------------- Peng
# Besca-reprocessed release, Zenodo 10.5281/zenodo.3969339, reprocessed from
# BIGD PRJCA001063.  Annotations are Besca-derived, NOT the authors' own (A.b).
# Raw GSA CRA001160 was rejected: it carries no cell-type annotation.
# The "CRC" in the filename is a workflow naming artifact -- the data are
# pancreatic, and the script halts if a colorectal vocabulary is seen instead.
PENG_MD5 <- "41fb7b9f27b7bb613ff979baaac5272f"

load_Peng <- function() {
  sec <- "A.a/Peng"
  fp  <- need_file(file.path(RAW, "StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"), sec)

  if (requireNamespace("tools", quietly = TRUE)) {
    md5 <- unname(tools::md5sum(fp))
    if (!identical(md5, PENG_MD5))
      halt(sec, "md5 mismatch: expected ", PENG_MD5, ", observed ", md5,
           ". This is not the registered Besca release.")
    message("  ok  A.a/Peng        md5 ", PENG_MD5, " verified")
  }

  obs <- rhdf5::h5ls(fp, recursive = TRUE)
  present <- paste0(sub("^/", "", obs$group), "/", obs$name)
  for (nm in c("obs/CONDITION", "obs/Patient", "obs/celltype0", "obs/celltype1")) {
    if (!nm %in% present) halt(sec, "required field missing: ", nm)
  }

  cond  <- h5ad_factor(fp, "CONDITION")
  pat   <- h5ad_factor(fp, "Patient")
  ct0   <- h5ad_factor(fp, "celltype0")
  ct1   <- h5ad_factor(fp, "celltype1")
  assert_n(length(ct0), 57423L, sec, "cells in the Besca release")
  assert_n(length(unique(pat)), 35L, sec, "patients (24 tumour + 11 normal)")

  # Explicit anti-colorectal guard, per the A.b verification requirement.
  if (any(grepl("colorect|colonocyte|enterocyte|goblet|CDX2", ct1, ignore.case = TRUE)))
    halt(sec, "cell-type vocabulary looks COLORECTAL, not pancreatic: ",
         paste(sort(unique(ct1)), collapse = ", "),
         ". The 'CRC' in the filename should be a naming artifact only. STOP.")
  if (!any(grepl("pancrea|acinar|stellate|ductal", ct1, ignore.case = TRUE)))
    halt(sec, "no pancreatic cell type found in celltype1: ",
         paste(sort(unique(ct1)), collapse = ", "))

  keep <- trimws(cond) == "T"
  assert_n(length(unique(pat[keep])), 24L, sec, "tumour patients")
  assert_n(length(unique(pat[!keep])), 11L, sec, "normal patients (excluded)")

  # celltype0 for the epithelial/stromal/endothelial branches; celltype1 splits
  # `hematopoietic cell` into myeloid vs lymphoid, which celltype0 collapses.
  lab <- ct0
  hem <- ct0 == "hematopoietic cell"
  comp <- character(length(lab))
  comp[!hem] <- map_labels(lab[!hem], MAP_Peng[names(MAP_Peng) != "hematopoietic cell"], "Peng/celltype0")
  comp[hem]  <- map_labels(ct1[hem],  MAP_Peng_ct1, "Peng/celltype1(haematopoietic)")

  X <- read_h5ad_counts(fp, sec)
  X <- X[, keep, drop = FALSE]

  verify_tissue("Peng", "pancreatic", rownames(X), ct1)

  cells <- data.frame(
    cell_id     = colnames(X),
    patient_id  = pat[keep],
    compartment = comp[keep],
    stringsAsFactors = FALSE
  )
  list(X = X, cells = cells, atlas = "Peng", tissue = "pancreatic")
}

# h5ad X -> genes x cells raw counts.  AnnData stores cells x genes CSR; this
# returns the transpose without normalising anything.
read_h5ad_counts <- function(path, sec) {
  src <- if (length(rhdf5::h5ls(path, recursive = TRUE)$name[
               rhdf5::h5ls(path, recursive = TRUE)$group == "/raw"]) > 0) "raw/X" else "X"
  ls_h5 <- rhdf5::h5ls(path, recursive = TRUE)
  fields <- ls_h5$name[ls_h5$group == paste0("/", src)]
  if (!all(c("data", "indices", "indptr") %in% fields))
    halt(sec, "h5ad '", src, "' is not sparse CSR/CSC; found: ", paste(fields, collapse = ", "))
  vg <- if (src == "raw/X") "raw/var" else "var"
  gsym <- as.character(rhdf5::h5read(path, paste0(vg, "/index")))
  bc   <- as.character(rhdf5::h5read(path, "obs/index"))
  M <- Matrix::sparseMatrix(j = as.integer(rhdf5::h5read(path, paste0(src, "/indices"))) + 1L,
                            p = as.integer(rhdf5::h5read(path, paste0(src, "/indptr"))),
                            x = as.numeric(rhdf5::h5read(path, paste0(src, "/data"))),
                            dims = c(length(bc), length(gsym)), index1 = FALSE)
  M <- Matrix::t(M)                 # -> genes x cells
  rownames(M) <- gsym
  colnames(M) <- bc
  M
}

# ==============================================================================
# A.d  PSEUDOBULK  -- RAW counts summed; normalisation strictly AFTER summing
# ==============================================================================
# FORBIDDEN: CP10K, log, scaling, per-compartment means, or any library-size
# correction BEFORE the rowSums.  The estimand is each compartment's
# contribution to a bulk library = cells x transcripts per cell.  Normalising
# first discards the transcripts-per-cell term and reweights every compartment
# to equal RNA content -- the error that made the HPA pilot measure the wrong
# quantity (feasibility_assessment.md addendum).

pseudobulk_raw <- function(X, cells, genes) {
  S <- matrix(0, nrow = length(genes), ncol = length(COMPARTMENTS),
              dimnames = list(genes, COMPARTMENTS))
  n_cells <- setNames(integer(length(COMPARTMENTS)), COMPARTMENTS)
  have <- intersect(genes, rownames(X))
  for (cc in COMPARTMENTS) {
    idx <- which(cells$compartment == cc)
    n_cells[cc] <- length(idx)
    if (!length(idx) || !length(have)) next
    S[have, cc] <- as.numeric(Matrix::rowSums(X[have, idx, drop = FALSE]))  # RAW
  }
  # Genes absent from this atlas are NA, never 0 (A.d).
  S[setdiff(genes, have), ] <- NA_real_
  list(counts = S, n_cells = n_cells)
}

EVIDENCE_MIN <- 20L   # prespecified; genes below this are insufficient-evidence

# ==============================================================================
# A.e  PURITY SWEEP
# ==============================================================================
PI_GRID <- seq(0.30, 0.70, by = 0.01)
if (length(PI_GRID) != 41L) halt("A.e", "purity grid must have 41 points")

# I[g,c] = S[g,c] / n_c  -- mean counts per cell, retaining between-compartment
# differences in total RNA content (intended: a neutrophil and a carcinoma cell
# do not contribute equal transcript mass).
#
# w_epi = pi ; w_c = (1-pi) * n_c / sum(n_c over non-epithelial)
# f(pi) = w_epi*I[,epi] / sum_c w_c*I[,c]
sweep_f <- function(S, n_cells, pi_grid = PI_GRID) {
  I <- sweep(S, 2, pmax(n_cells, 1L), "/")
  I[, n_cells == 0L] <- 0
  non_epi <- setdiff(COMPARTMENTS, "epithelial")
  n_non   <- sum(n_cells[non_epi])
  if (n_non == 0L) halt("A.e", "no non-epithelial cells; f(pi) is undefined")
  share   <- n_cells[non_epi] / n_non

  vapply(pi_grid, function(pi) {
    w <- c(setNames(pi, "epithelial"), setNames((1 - pi) * share, non_epi))[COMPARTMENTS]
    num <- I[, "epithelial"] * w[["epithelial"]]
    den <- as.numeric(I %*% w)
    ifelse(den > 0, num / den, NA_real_)
  }, numeric(nrow(S)))
}

# f(pi) is monotone increasing in pi.  Asserted numerically, never assumed.
assert_monotone <- function(F, atlas_id) {
  ok <- apply(F, 1, function(r) {
    if (all(is.na(r))) return(TRUE)
    all(diff(r[!is.na(r)]) >= -1e-9)
  })
  if (!all(ok))
    halt("A.e", atlas_id, ": f(pi) is not monotone increasing for ", sum(!ok),
         " gene(s), e.g. ", paste(head(rownames(F)[!ok], 5), collapse = ", "),
         ". Check the weights -- do not relax this check.")
  message("  ok  A.e  ", atlas_id, " f(pi) monotone over all ", ncol(F), " grid points")
  invisible(TRUE)
}

# Dominance: f(pi) > 0.50 at EVERY grid point.  Equivalent to f(0.30) > 0.50 by
# monotonicity, but the full grid is evaluated and reported.
dominance_from_F <- function(F, evidence_ok) {
  d <- apply(F, 1, function(r) if (all(is.na(r))) NA else all(r > 0.50))
  d[!evidence_ok] <- NA
  d
}

# ==============================================================================
# A.f  PATIENT-LEVEL BOOTSTRAP
# ==============================================================================
# Resampling unit: PATIENT.  Not cell, not GSM, not channel.  Cells are not
# independent within a patient and several GSE178341 GSMs share a specimen, so
# resampling cells would give intervals far too narrow.
# The whole A.d -> A.e chain is recomputed inside each resample, not just the
# final ratio.
B_RESAMPLES <- 2000L
SEED_BASE   <- 20260731L

one_atlas_stats <- function(X, cells, genes) {
  pb <- pseudobulk_raw(X, cells, genes)
  tot <- rowSums(pb$counts, na.rm = TRUE)
  evidence_ok <- !is.na(pb$counts[, 1]) & tot >= EVIDENCE_MIN
  F  <- sweep_f(pb$counts, pb$n_cells)
  rownames(F) <- genes
  list(pb = pb, F = F, evidence_ok = evidence_ok,
       dom = dominance_from_F(F, evidence_ok),
       dom50 = { d <- F[, which.min(abs(PI_GRID - 0.50))] > 0.50
                 d[!evidence_ok] <- NA; d })
}

bootstrap_atlas <- function(X, cells, genes, atlas_index, B = B_RESAMPLES) {
  pats <- unique(cells$patient_id)
  by_pat <- split(seq_len(nrow(cells)), cells$patient_id)
  withr::with_seed(SEED_BASE + atlas_index, {
    reps <- lapply(seq_len(B), function(b) {
      draw <- sample(pats, length(pats), replace = TRUE)
      idx  <- unlist(by_pat[draw], use.names = FALSE)
      st   <- one_atlas_stats(X[, idx, drop = FALSE], cells[idx, , drop = FALSE], genes)
      list(f30 = st$F[, 1], dom = st$dom, dom50 = st$dom50)
    })
  })
  list(
    f30_lo = apply(vapply(reps, `[[`, numeric(length(genes)), "f30"), 1, quantile, 0.025, na.rm = TRUE),
    f30_hi = apply(vapply(reps, `[[`, numeric(length(genes)), "f30"), 1, quantile, 0.975, na.rm = TRUE),
    dom_reps   = vapply(reps, `[[`, logical(length(genes)), "dom"),
    dom50_reps = vapply(reps, `[[`, logical(length(genes)), "dom50")
  )
}

# ==============================================================================
# A.g  k AND ITS VARIANTS  (Amendments 5, 9, 10)
# ==============================================================================
# k         : epithelial-dominant in >= 2 of the 3 atlases, any GI tissue
# k_all3    : dominant in ALL three (and evaluable in all three)
# k_evalall : k restricted to genes evaluable in all three
# k_50      : as k, but dominance evaluated at pi = 0.50 only (Amendment 6)
compute_k <- function(dominance) {
  if (ncol(dominance) != 3L) halt("A.g", "dominance matrix must have exactly 3 atlas columns")
  n_dom  <- rowSums(dominance, na.rm = TRUE)
  n_eval <- rowSums(!is.na(dominance))
  k         <- sum(n_dom >= 2L)
  k_all3    <- sum(n_dom == 3L & n_eval == 3L)
  k_evalall <- sum(n_dom >= 2L & n_eval == 3L)
  stopifnot(k_all3 <= k_evalall, k_evalall <= k)
  list(k = k, k_all3 = k_all3, k_evalall = k_evalall,
       n_dom = n_dom, n_eval = n_eval)
}

# Amendment 3 branch, on primary k against the 152-gene panel.
# Wording is Amendment 3's own, not a paraphrase: the middle band is EXPLORATORY
# ONLY, which is materially weaker than "reported with caveats".  Both bounds are
# proportions of the final panel; k >= 8 is the amendment's additional floor on
# the top band, and coincides with ceiling(5% of 152) at the bottom band.
branch_of_k <- function(k, n_panel = 152L) {
  hi <- ceiling(0.20 * n_panel)   # 31 at n_panel = 152
  lo <- ceiling(0.05 * n_panel)   # 8  at n_panel = 152
  if (k >= hi && k >= 8L)
    "epithelial subscore BUILT; subscore survival models run as SECONDARY"
  else if (k >= lo)
    "subscore reported as EXPLORATORY only"
  else
    "decomposition DESCRIPTIVE only"
}

# ==============================================================================
# DRIVER
# ==============================================================================
message("\n== A.a/A.b  loading the three registered atlases ==")
LOADERS <- list(GSE125449 = load_GSE125449,
                GSE178341 = load_GSE178341,
                Peng      = load_Peng)
if (!setequal(names(LOADERS), atlases$atlas))
  halt("A.a", "loader set does not match the registered atlas set")

dat <- lapply(atlases$atlas, function(a) LOADERS[[a]]())
names(dat) <- atlases$atlas

message("\n== A.d/A.e  pseudobulk and purity sweep ==")
stats <- lapply(names(dat), function(a) {
  st <- one_atlas_stats(dat[[a]]$X, dat[[a]]$cells, PANEL_GENES)
  assert_monotone(st$F, a)
  message(sprintf("  ..  %-12s evaluable genes = %d / 152", a, sum(st$evidence_ok)))
  st
})
names(stats) <- names(dat)

dominance    <- vapply(stats, `[[`, logical(152), "dom")
dominance_50 <- vapply(stats, `[[`, logical(152), "dom50")
rownames(dominance) <- rownames(dominance_50) <- PANEL_GENES
colnames(dominance) <- colnames(dominance_50) <- names(stats)

message("\n== A.g  k and variants ==")
K   <- compute_k(dominance)
K50 <- compute_k(dominance_50)
message(sprintf("  k = %d | k_all3 = %d | k_evalall = %d | k_50 = %d",
                K$k, K$k_all3, K$k_evalall, K50$k))
message("  Amendment 3 branch (on primary k): ", branch_of_k(K$k))

message("\n== A.f  patient-level bootstrap (B = ", B_RESAMPLES, ", unit = patient) ==")
boots <- lapply(seq_along(dat), function(i)
  bootstrap_atlas(dat[[i]]$X, dat[[i]]$cells, PANEL_GENES, atlas_index = i))
names(boots) <- names(dat)

k_reps <- vapply(seq_len(B_RESAMPLES), function(b) {
  dm <- vapply(boots, function(bt) bt$dom_reps[, b], logical(152))
  rownames(dm) <- PANEL_GENES
  unlist(compute_k(dm)[c("k", "k_all3", "k_evalall")])
}, numeric(3))
k_ci <- apply(k_reps, 1, quantile, c(0.025, 0.975), na.rm = TRUE)

# ==============================================================================
# OUTPUTS  (A.g required reporting)
# ==============================================================================
f30 <- vapply(stats, function(s) s$F[, 1], numeric(152))
dom_long <- do.call(rbind, lapply(names(stats), function(a) data.frame(
  gene        = PANEL_GENES,
  atlas       = a,
  tissue      = atlases$tissue[atlases$atlas == a],
  f_at_0.30   = stats[[a]]$F[, 1],
  f_at_0.50   = stats[[a]]$F[, which.min(abs(PI_GRID - 0.50))],
  f_at_0.70   = stats[[a]]$F[, length(PI_GRID)],
  ci_lo_0.30  = boots[[a]]$f30_lo,
  ci_hi_0.30  = boots[[a]]$f30_hi,
  evaluable   = stats[[a]]$evidence_ok,
  dominant    = stats[[a]]$dom,
  dominant_50 = stats[[a]]$dom50,
  origin_six  = PANEL_GENES %in% ORIGIN_SIX,
  stringsAsFactors = FALSE)))

write.csv(dom_long, file.path(OUTDIR, "compartment_dominance_matrix.csv"), row.names = FALSE)

evaluability <- table(factor(K$n_eval, levels = 3:0))
write.csv(data.frame(atlases_evaluable = names(evaluability),
                     n_genes = as.integer(evaluability)),
          file.path(OUTDIR, "evaluability_distribution.csv"), row.names = FALSE)

per_tissue <- vapply(names(stats), function(a) sum(stats[[a]]$dom, na.rm = TRUE), integer(1))
write.csv(data.frame(atlas = names(per_tissue),
                     tissue = atlases$tissue[match(names(per_tissue), atlases$atlas)],
                     n_dominant = as.integer(per_tissue)),
          file.path(OUTDIR, "per_tissue_dominance.csv"), row.names = FALSE)

k_tab <- data.frame(
  quantity = c("k", "k_all3", "k_evalall", "k_50"),
  estimate = c(K$k, K$k_all3, K$k_evalall, K50$k),
  ci_lo    = c(k_ci[1, "k"], k_ci[1, "k_all3"], k_ci[1, "k_evalall"], NA),
  ci_hi    = c(k_ci[2, "k"], k_ci[2, "k_all3"], k_ci[2, "k_evalall"], NA),
  stringsAsFactors = FALSE)
write.csv(k_tab, file.path(OUTDIR, "k_estimates.csv"), row.names = FALSE)

# Origin six as a labelled subset: SOCS3/MYC/IL6 in panel, BCL2/MMP9/HGF not.
origin_tab <- dom_long[dom_long$gene %in% ORIGIN_SIX, ]
origin_tab$qualifying <- origin_tab$gene %in% ORIGIN_IN_PANEL
write.csv(origin_tab, file.path(OUTDIR, "origin_six_compartment.csv"), row.names = FALSE)

writeLines(c(
  paste0("run_date: ", Sys.Date()),
  paste0("plan_version: 1.5"),
  paste0("atlases: ", paste(atlases$atlas, collapse = ", ")),
  paste0("panel_genes: ", length(PANEL_GENES)),
  paste0("bootstrap: B=", B_RESAMPLES, " unit=patient seed_base=", SEED_BASE),
  paste0("k=", K$k, " k_all3=", K$k_all3, " k_evalall=", K$k_evalall, " k_50=", K50$k),
  paste0("branch: ", branch_of_k(K$k))
), file.path(OUTDIR, "compartment_provenance.txt"))

message("\nPart A complete (A.a-A.g). HARD STOP: no gene list locked, no score ",
        "constructed, no survival model fitted, no outcome data opened.")
sessionInfo()
