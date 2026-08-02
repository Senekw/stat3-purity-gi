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

for (p in c("Matrix", "rhdf5", "withr", "tools")) {   # F12: tools is load-bearing (md5)
  if (!requireNamespace(p, quietly = TRUE))
    stop("HALT [deps]: package '", p, "' is required and not installed.", call. = FALSE)
}

halt <- function(section, ...) {
  stop(paste0("HALT [", section, "]: ", paste0(c(...), collapse = "")), call. = FALSE)
}

# F13: every loader returns a matrix and a cells frame that must stay aligned --
# columns of X are the cells, in the same order.  Each loader builds that
# correspondence differently (GSE125449 by composite name, GSE178341 by name from
# a deposit-wide vector, Peng positionally from a logical), so the invariant is
# asserted once, here, rather than trusted three times.
assert_cells_aligned <- function(X, cells, sec) {
  if (anyDuplicated(colnames(X)))
    halt(sec, "duplicate cell ids in X: ",
         paste(utils::head(unique(colnames(X)[duplicated(colnames(X))]), 5), collapse = ", "))
  assert_n(ncol(X), nrow(cells), sec, "matrix columns == cells rows")
  if (!identical(colnames(X), cells$cell_id))
    halt(sec, "X columns and cells$cell_id are not aligned. Every downstream ",
         "quantity indexes cells positionally against X; a misalignment would ",
         "attribute one cell's counts to another compartment or patient.")
  invisible(TRUE)
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

# Tissue witnesses.  Its PURPOSE is to catch a transposed tissue assignment, and
# it does that; it is not a general tissue classifier.
#
# What is actually checked (S2 -- the previous comment here overstated it by
# saying marker genes must be "detected"):
#   (i)  marker-gene PRESENCE in the atlas's gene universe -- membership only, no
#        expression is read.  For a whole-transcriptome human reference all three
#        witnesses' genes are present, so this term is effectively constant and
#        contributes no discrimination between tissues.  It is a precondition,
#        not evidence.
#   (ii) the LABEL pattern, which must match at least one of the atlas's own
#        cell-type labels.  This term carries the whole discrimination: the
#        patterns are mutually exclusive across the three declared tissues, so a
#        transposed row cannot satisfy another atlas's witness.
# Scoring weights the label term 10x the gene term for exactly that reason.
# Verified on real labels: GSE125449 scores liver_biliary = 14 against
# colorectal = 4 and pancreatic = 2, so claiming either other tissue halts.
# Note the tie-break is `<` against the max, so a tie passes -- adequate here
# because the label patterns cannot both match, and deliberately left alone.
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

  # Guard against being handed a ROW-SUBSET matrix's rownames instead of the
  # atlas's full gene universe.  No witness gene is a panel gene, so a subset
  # would fail the marker test for reasons that have nothing to do with tissue.
  # A real transcriptome has thousands of genes; the reporting set has 155.
  if (length(gene_universe) < 1000L)
    halt("A.a", atlas_id, ": verify_tissue received ", length(gene_universe),
         " genes, which looks like a row-subset rather than the atlas's full gene ",
         "universe. Witness genes are not panel genes, so this would halt spuriously. ",
         "Pass the complete feature list.")

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

# --- S1: the REPORTING set (155) vs the INFERENTIAL set (152) -----------------
# Prespecification section 4 and Amendment 2 both commit to reporting BCL2, MMP9
# and HGF as a LABELLED NON-QUALIFYING SUBSET -- never silently dropped.  A
# previous version filtered the output table from the 152-gene panel, so those
# three could not appear at all and the `qualifying` column was TRUE in every row.
#
# Two sets are therefore carried, and the distinction is load-bearing:
#   PANEL_GENES  (152) -- the locked panel. The ONLY set any k variant sees.
#   REPORT_GENES (155) -- panel + the three non-qualifying, for compartment
#                         fractions and the origin-six table ONLY.
# The three extra genes are reported with fractions but must never enter k, k_50,
# k_all3, k_evalall, the dominance matrix or the evaluability distribution.
# assert_inferential_set() below enforces that at every point where it matters.
# F18: no duplicate check here.  ORIGIN_NONQUAL is setdiff(ORIGIN_SIX,
# PANEL_GENES) and PANEL_GENES is asserted unique above, so REPORT_GENES cannot
# contain a duplicate by construction -- a guard on it would read as protection
# while being unreachable.  The length assertion below is the real check.
REPORT_GENES <- c(PANEL_GENES, ORIGIN_NONQUAL)
assert_n(length(REPORT_GENES), 155L, "panel", "reporting genes (152 panel + 3 non-qualifying)")

# The inferential firewall.  Called on every gene-indexed object that feeds a k
# variant.  Cheap, and it makes the 152/155 boundary impossible to cross silently.
assert_inferential_set <- function(x, what) {
  g <- if (is.null(dim(x))) names(x) else rownames(x)
  if (is.null(g)) halt("A.g", what, " has no gene names; cannot verify the inferential set")
  if (!identical(g, PANEL_GENES))
    halt("A.g", what, " is not exactly the 152-gene locked panel in order (n = ", length(g),
         if (length(setdiff(g, PANEL_GENES)))
           paste0("; extra: ", paste(utils::head(setdiff(g, PANEL_GENES), 5), collapse = ", "))
         else "",
         if (length(setdiff(PANEL_GENES, g)))
           paste0("; missing: ", paste(utils::head(setdiff(PANEL_GENES, g), 5), collapse = ", "))
         else "",
         "). The non-qualifying subset must never enter a k variant.")
  invisible(TRUE)
}

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

# Amendment 3's bands are verified here too -- before any atlas is opened -- but
# the definition lives with branch_of_k in A.g.  Deferred to a call site after
# that definition; see `verify_branch_bands(length(PANEL_GENES))` below.

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

# --- Peng (Besca celltype1): AMENDMENT 11, 2026-08-01.
# Compartment mapping uses celltype1, not celltype0.  celltype0 collapses myeloid
# and lymphoid into a single `hematopoietic cell` category and therefore cannot
# express Amendment 4's six-compartment scheme at all; celltype1 separates
# myeloid leukocyte, T lineage and B lineage and maps onto the six directly.
#
# celltype1 is a strict refinement of celltype0 (verified: every celltype1 label
# nests inside exactly one celltype0 label), so this subdivides the non-epithelial
# mass without moving any cell across the epithelial boundary.  Amendment 11's
# stated direction of bias -- none -- is verified numerically in A.c below rather
# than argued.
MAP_Peng <- c(
  "pancreatic ductal cell"           = "epithelial",          # ductal
  "pancreatic acinar cell"           = "epithelial",          # acinar, per A.c judgement call
  "enteroendocrine cell"             = "epithelial",          # endocrine/islet, per A.c judgement call
  "fibroblast"                       = "fibroblast_stromal",
  "pancreatic stellate cell"         = "fibroblast_stromal",  # stellate = pancreatic CAF lineage
  "blood vessel endothelial cell"    = "endothelial",
  "myeloid leukocyte"                = "myeloid",
  "T cell"                           = "lymphoid",
  "lymphocyte of B lineage"          = "lymphoid",
  "neural cell"                      = "other"                # neural crest lineage, not stromal (A.c)
)

# The celltype0 map is retained ONLY to verify Amendment 11's bias claim (A.c
# below).  It is not used to build compartments.  `hematopoietic cell` is
# deliberately unmappable at this level -- that is precisely why Amendment 11
# moves to celltype1.
MAP_Peng_ct0_epi_only <- c(
  "epithelial cell"    = "epithelial",
  "fibroblast"         = "fibroblast_stromal",
  "endothelial cell"   = "endothelial",
  "hematopoietic cell" = "hematopoietic_UNRESOLVED",
  "neural cell"        = "other"
)

COMPARTMENT_MAPS <- list(GSE125449 = MAP_GSE125449,
                         GSE178341 = MAP_GSE178341,
                         Peng      = MAP_Peng)

# verify_amendment11_bias: Amendment 11 states its direction of bias is "none,
# and this is provable rather than argued" -- f(pi) is numerically identical under
# celltype0 and celltype1 because subdividing the non-epithelial mass changes
# neither the numerator nor the denominator of the estimand.
#
# That claim rests on ONE structural fact: the epithelial cell set must be exactly
# the same under both annotation levels.  Asserted here on the real labels, not
# assumed.  If it fails, the amendment's reasoning is falsified and the run halts.
verify_amendment11_bias <- function(ct0, ct1) {
  bad0 <- setdiff(unique(ct0), names(MAP_Peng_ct0_epi_only))
  if (length(bad0))
    halt("A.c", "Peng celltype0 has unmapped label(s): ", paste(sort(bad0), collapse = ", "))
  # F14: check celltype1 coverage too.  Without this, an unmapped ct1 label makes
  # MAP_Peng[ct1] return NA, epi_ct1 carries NA, identical() is FALSE, and the run
  # halts with "AMENDMENT 11 BIAS CLAIM FALSIFIED" -- blaming the amendment for
  # what is really an incomplete map.
  bad1 <- setdiff(unique(ct1), names(MAP_Peng))
  if (length(bad1))
    halt("A.c", "Peng celltype1 has unmapped label(s): ", paste(sort(bad1), collapse = ", "),
         ". Extend MAP_Peng as a dated amendment; this is a map-completeness ",
         "problem, not a falsification of Amendment 11.")
  epi_ct0 <- unname(MAP_Peng_ct0_epi_only[ct0]) == "epithelial"
  epi_ct1 <- unname(MAP_Peng[ct1]) == "epithelial"
  if (!identical(epi_ct0, epi_ct1))
    halt("A.c", "AMENDMENT 11 BIAS CLAIM FALSIFIED: the epithelial cell set differs ",
         "between celltype0 and celltype1 (", sum(epi_ct0), " vs ", sum(epi_ct1),
         " cells, ", sum(xor(epi_ct0, epi_ct1)), " disagreeing). Amendment 11 states the ",
         "direction of bias is none BECAUSE subdividing non-epithelial mass cannot move ",
         "the epithelial/non-epithelial boundary. It has moved. STOP and report.")

  # celltype1 must also be a strict refinement: each celltype1 label nests inside
  # exactly one celltype0 label.  Otherwise "subdivision" is the wrong description.
  cross <- table(ct1, ct0)
  multi <- rownames(cross)[rowSums(cross > 0) > 1L]
  if (length(multi))
    halt("A.c", "celltype1 is not a strict refinement of celltype0: label(s) ",
         paste(multi, collapse = ", "), " span more than one celltype0 category.")

  message(sprintf("  ok  A.c  Amendment 11 bias claim verified: epithelial set identical under celltype0/celltype1 (%d cells); celltype1 is a strict refinement",
                  sum(epi_ct1)))
  invisible(TRUE)
}

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
    rownames(m) <- g[[ncol(g)]]            # symbol column
    # F01: Set1 and Set2 are independent 10x runs drawing from the same barcode
    # whitelist.  Bare barcodes as column names would let cbind produce duplicate
    # colnames, and name-based column selection returns the FIRST match for every
    # request -- each colliding Set2 cell would silently receive a Set1 cell's
    # expression profile while keeping its own Sample and Type.  Measured on this
    # deposit the collision count is 0 (the sets carry distinct -1/-2 GEM-well
    # suffixes), so this is a latent defect, not an active one; the composite key
    # makes it unreachable regardless of how a future deposit is barcoded.
    colnames(m) <- paste0(s, "|", b[[1]])
    m
  })
  names(mats) <- sets

  # AMENDMENT 12, 2026-08-01: the two sets have DIFFERENT gene universes
  # (Set1 20,124 rows; Set2 19,572), so they are combined on the INTERSECTION.
  # Union with 0-fill would record an unmeasured gene as measured-at-zero,
  # violating A.d's "absence is NA, never 0"; union with NA-fill would give the
  # two sets different denominators. Intersection keeps every retained value a
  # measured one.
  assert_n(nrow(mats$Set1), 20124L, sec, "Set1 gene rows as deposited")
  assert_n(nrow(mats$Set2), 19572L, sec, "Set2 gene rows as deposited")

  # F02: the duplicate-symbol check must run on the SOURCE row names, not on the
  # intersection.  base::intersect() calls unique() on both arguments, so
  # anyDuplicated(intersect(...)) is structurally always 0 and the guard that used
  # to sit here could never fire.  The hazard is real but lives one line earlier:
  # CellRanger symbol columns repeat symbols, and name-based row selection would
  # silently take only the first.  Checked per set, restricted to the reporting
  # genes -- a duplicate elsewhere in the transcriptome cannot affect this analysis.
  for (s in sets) {
    rn  <- rownames(mats[[s]])
    dup <- unique(rn[duplicated(rn) & rn %in% REPORT_GENES])
    if (length(dup))
      halt(sec, s, " has duplicate row(s) for reporting gene symbol(s): ",
           paste(utils::head(sort(dup), 10), collapse = ", "),
           ". Name-based row selection would take only the first and undercount the ",
           "rest. Decide the aggregation rule as a dated amendment; do not sum silently.")
  }
  common <- intersect(rownames(mats$Set1), rownames(mats$Set2))
  assert_n(length(common), 18367L, sec, "genes in the Set1/Set2 intersection (Amendment 12)")

  # The six panel genes Amendment 12 names must be exactly the ones lost, so a
  # future deposit revision cannot silently change coverage.
  AMDT12_EXCLUDED <- c("CCL7", "CRLF2", "CSF2", "IL9R", "ITGB3", "LEP")
  in_either <- union(rownames(mats$Set1), rownames(mats$Set2))
  lost <- sort(setdiff(intersect(PANEL_GENES, in_either), common))
  if (!identical(lost, sort(AMDT12_EXCLUDED)))
    halt(sec, "Amendment 12 names CCL7, CRLF2, CSF2, IL9R, ITGB3, LEP as the panel genes ",
         "excluded by the intersection, but the observed set is ",
         if (length(lost)) paste(lost, collapse = ", ") else "(none)",
         ". The deposit's gene universes have changed; re-check coverage before proceeding.")
  # F04: test the origin-six claim over the REPORTING set, not the panel.  `lost`
  # above is derived from intersect(PANEL_GENES, ...) and BCL2/MMP9/HGF are by
  # definition not in PANEL_GENES, so the previous form of this guard could only
  # ever have caught SOCS3, MYC or IL6 -- half the claim its own message makes.
  lost_report <- sort(setdiff(intersect(REPORT_GENES, in_either), common))
  if (length(intersect(lost_report, ORIGIN_SIX)))
    halt(sec, "an origin-six gene is excluded by the Amendment 12 intersection: ",
         paste(sort(intersect(lost_report, ORIGIN_SIX)), collapse = ", "),
         ". Amendment 12 states no origin-six gene is affected.")
  assert_n(sum(PANEL_GENES %in% common), 143L, sec, "panel genes retained in GSE125449 (Amendment 12)")
  assert_n(sum(REPORT_GENES %in% common), 146L, sec, "reporting genes retained in GSE125449")
  message("  ok  A.a/GSE125449   Amendment 12 intersection: excluded ",
          paste(lost, collapse = ", "), " (panel coverage 143/152)")

  X <- cbind(mats$Set1[common, , drop = FALSE], mats$Set2[common, , drop = FALSE])
  if (anyDuplicated(colnames(X)))
    halt(sec, "duplicate cell ids after cbind: ",
         paste(utils::head(unique(colnames(X)[duplicated(colnames(X))]), 5), collapse = ", "),
         ". Name-based column selection would return the first match for every ",
         "request, assigning one cell's counts to another.")

  # composite key, matching the set-prefixed colnames above (F01)
  meta$cell_id <- paste0(meta$set, "|", meta$`Cell Barcode`)
  if (anyDuplicated(meta$cell_id))
    halt(sec, "duplicate set|barcode in the per-cell metadata: ",
         paste(utils::head(unique(meta$cell_id[duplicated(meta$cell_id)]), 5), collapse = ", "))
  keep <- meta$cell_id %in% colnames(X)
  if (!any(keep)) halt(sec, "no metadata barcode matches a matrix column")
  # assert the JOIN is complete, not merely non-empty: a silent drop here would
  # shrink the cell set without any count changing visibly.
  assert_n(sum(keep), nrow(meta), sec, "metadata cells found in the matrix")
  meta <- meta[keep, , drop = FALSE]
  X    <- X[, meta$cell_id, drop = FALSE]

  verify_tissue("GSE125449", "liver_biliary", rownames(X), meta$Type)

  cells <- data.frame(
    cell_id     = meta$cell_id,
    patient_id  = meta$Sample,                       # sample == patient here
    compartment = map_labels(meta$Type, MAP_GSE125449, "GSE125449"),
    # retained for Amendment 4's malignant-restricted sensitivity (F07): this is
    # the only atlas whose labels separate CNV-called malignant cells from other
    # epithelium, so the raw label is needed downstream
    cell_type_raw = as.character(meta$Type),
    stringsAsFactors = FALSE
  )
  assert_n(length(unique(cells$patient_id)), 19L, sec, "bootstrap units (patients)")
  assert_cells_aligned(X, cells, sec)
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
  # The file carries a provenance header (source URL, download date, derivation).
  # Strip comments and blanks before counting -- the header is documentation, not
  # data. Channel ids never begin with '#'.
  geo_channels <- trimws(readLines(f_geo))
  geo_channels <- unique(geo_channels[nzchar(geo_channels) & !startsWith(geo_channels, "#")])
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
  # Row-subset reader: the full 43,113 x 370,115 matrix (764M nonzeros) does not
  # fit in 16 GB, and only the reporting genes are ever read.
  X   <- read_h5_counts_subset(f_h5, sec, REPORT_GENES)
  X   <- X[, meta$cellID[keep_cell], drop = FALSE]

  # verify_tissue needs the atlas's FULL gene universe, not the row-subset X
  # returned by read_h5_counts_subset -- none of the colorectal witness genes is a
  # panel gene, so passing rownames(X) would halt spuriously.  The feature list is
  # cheap metadata and is read directly.
  full_gene_universe <- as.character(rhdf5::h5read(f_h5, "matrix/features/name"))
  verify_tissue("GSE178341", "colorectal", full_gene_universe, lab)

  cells <- data.frame(
    cell_id     = meta$cellID[keep_cell],
    patient_id  = pid,                                   # PID, never PatientTypeID
    compartment = map_labels(lab, MAP_GSE178341, "GSE178341"),
    stringsAsFactors = FALSE
  )
  assert_cells_aligned(X, cells, sec)
  list(X = X, cells = cells, atlas = "GSE178341", tissue = "colorectal")
}

# ------------------------------------------------------------------ helpers
# assert_raw_counts: A.d requires RAW UMI counts.  Raw counts are integer-valued;
# any library-size correction, CP10K or log transform makes them fractional.
# This is the guard the estimand rests on: pseudobulking normalised values is the
# error that made the HPA pilot measure the wrong quantity, and nothing
# downstream -- monotonicity, dominance, the k ordering -- can detect it, because
# normalised values are positive and plausibly scaled.
assert_raw_counts <- function(M, sec, src) {
  v <- M@x
  if (!length(v)) halt(sec, "matrix '", src, "' has no non-zero entries")
  if (any(!is.finite(v)))
    halt(sec, "matrix '", src, "' contains non-finite values")
  if (any(v < 0))
    halt(sec, "matrix '", src, "' contains negative values; raw UMI counts cannot be negative")
  if (!all(abs(v - round(v)) < 1e-8))
    halt(sec, "matrix '", src, "' is not integer-valued and is therefore NOT raw counts ",
         "(", sum(abs(v - round(v)) >= 1e-8), " of ", length(v), " non-zero entries fractional; ",
         "e.g. ", paste(utils::head(v[abs(v - round(v)) >= 1e-8], 3), collapse = ", "), "). ",
         "A.d forbids CP10K, log or any library-size correction before the rowSums. ",
         "Point the loader at the raw layer; do not normalise here and do not relax this check.")
  message(sprintf("  ok  %-28s %-46s = %s", sec, paste0("'", src, "' integer-valued (raw counts)"),
                  format(length(v), big.mark = ",")))
  invisible(TRUE)
}

# 10x-style HDF5 -> dgCMatrix, raw counts.
# read_h5_counts_subset: chunked CSC reader that materialises ONLY the requested
# gene rows.  GSE178341's matrix is 43,113 x 370,115 with 764,460,511 nonzeros --
# the data vector alone is ~6 GB as doubles, and building the full dgCMatrix
# exhausts 16 GB before any analysis runs.  Only 155 of 43,113 rows are ever used,
# so the file is streamed in column blocks and non-panel rows are discarded as
# they are read.
#
# Verified bitwise identical to the direct full-matrix construction on a
# 3,000-cell slice (identical() TRUE, max abs diff 0, same nnz), comparing BY ROW
# POSITION rather than by name so duplicate symbols cannot mask a mismatch.
#
# Duplicate gene symbols: this deposit carries pseudoautosomal `_PAR_Y` feature
# rows that duplicate a symbol (CRLF2, CSF2RA, IL3RA, IL9R).  A full-file scan
# confirmed ALL FOUR carry zero nonzeros and zero total counts, so dropping the
# empty copy is numerically identical to summing it -- no quantitative decision is
# being made here and no amendment is required.  A duplicate symbol whose second
# row is NON-empty is a different matter and still halts (B3, pseudobulk_raw).
read_h5_counts_subset <- function(path, sec, genes, chunk_cells = 20000L) {
  grp   <- "matrix"
  gx    <- function(n, idx = NULL) rhdf5::h5read(path, paste0(grp, "/", n), index = idx)
  gsym  <- as.character(gx("features/name"))
  gid   <- as.character(gx("features/id"))
  bc    <- as.character(gx("barcodes"))
  shp   <- as.integer(gx("shape"))
  indptr <- as.numeric(gx("indptr"))
  if (length(gsym) != shp[1] || length(bc) != shp[2])
    halt(sec, "HDF5 feature/barcode lengths disagree with declared shape")

  keep <- which(gsym %in% genes)
  if (!length(keep)) halt(sec, "no requested gene is present in this atlas")

  # Drop duplicate-symbol rows ONLY when the extra copies are entirely empty.
  dup_syms <- unique(gsym[keep][duplicated(gsym[keep])])
  if (length(dup_syms)) {
    nnz_of <- vapply(keep, function(r) 0L, integer(1))   # filled during the scan
    names(nnz_of) <- as.character(keep)
  }
  row_map <- integer(shp[1]); row_map[keep] <- seq_along(keep)

  I <- integer(0); J <- integer(0); XV <- numeric(0)
  cols <- seq_len(shp[2])
  for (s in seq(1L, length(cols), by = chunk_cells)) {
    cc <- cols[s:min(s + chunk_cells - 1L, length(cols))]
    a  <- indptr[cc[1]] + 1; b <- indptr[cc[length(cc)] + 1]
    if (b < a) next
    ii  <- as.integer(gx("indices", list(a:b)))
    sel <- row_map[ii + 1L] > 0L
    if (any(sel)) {
      dd    <- as.numeric(gx("data", list(a:b)))[sel]
      colid <- rep(cc, times = diff(indptr[c(cc, cc[length(cc)] + 1)]))[sel]
      I <- c(I, row_map[ii[sel] + 1L]); J <- c(J, colid); XV <- c(XV, dd)
    }
    rm(ii); gc(FALSE)
  }
  M <- Matrix::sparseMatrix(i = I, j = J, x = XV,
                            dims = c(length(keep), shp[2]))
  rownames(M) <- gsym[keep]; colnames(M) <- bc

  if (length(dup_syms)) {
    per_row <- Matrix::rowSums(M != 0)
    drop <- integer(0)
    for (g in dup_syms) {
      w <- which(rownames(M) == g)
      empty <- w[per_row[w] == 0]
      if (length(w) - length(empty) != 1L)
        halt(sec, "gene symbol '", g, "' occupies ", length(w), " feature rows of which ",
             length(w) - length(empty), " carry counts (ids: ",
             paste(gid[keep][w], collapse = ", "), "). Only all-but-one-empty duplicates ",
             "can be resolved without a quantitative decision; this one cannot. ",
             "Decide the aggregation rule as a dated amendment.")
      drop <- c(drop, empty)
    }
    message("  ok  ", sec, "  dropped ", length(drop),
            " empty duplicate feature row(s) (", paste(sort(dup_syms), collapse = ", "),
            "; all zero-count PAR_Y copies)")
    if (length(drop)) M <- M[-drop, , drop = FALSE]
  }
  assert_raw_counts(M, sec, "matrix/data (row-subset)")
  M
}

# G08: read_h5_counts() -- the whole-matrix 10x HDF5 reader -- was REMOVED here.
# It had no call sites once GSE178341 moved to read_h5_counts_subset(), and it
# could not be used on this data anyway: materialising all 43,113 x 370,115
# columns exhausts 16 GB before any analysis starts, which is why the subset
# reader exists. A superseded reader kept "just in case" is a hazard, not a
# safety net -- someone reaches for the familiar name and reintroduces the
# memory failure. The 0-based-index note that lived in its body is preserved in
# read_h5_counts_subset(), which does the same sparseMatrix construction.

# h5ad categorical -> character, for the Besca-style layout used by Peng.
h5ad_factor <- function(path, col) {
  codes <- rhdf5::h5read(path, paste0("obs/", col))
  cats  <- try(rhdf5::h5read(path, paste0("obs/__categories/", col)), silent = TRUE)
  # F16: halt rather than silently returning integer codes as strings.  On the
  # AnnData >= 0.8 layout (obs/<col>/categories + obs/<col>/codes) the old
  # fallback would have decoded CONDITION to "0"/"1", making trimws(cond) == "T"
  # uniformly FALSE -- the tumour filter would select nothing and the failure
  # would surface far downstream as an unrelated count assertion.
  if (inherits(cats, "try-error"))
    halt("A.a/Peng", "obs/__categories/", col, " not found; this h5ad is not the ",
         "legacy Besca categorical layout the loader assumes. Handle the ",
         "AnnData >= 0.8 obs/<col>/categories layout explicitly before proceeding.")
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

  # F12: unconditional.  This was wrapped in requireNamespace("tools"), so on a
  # machine without `tools` the identity check on the one atlas whose counts are
  # RECONSTRUCTED rather than read would have been skipped silently.  `tools` is
  # now in the dependency loop at the top, so absence halts there.
  md5 <- unname(tools::md5sum(fp))
  if (!identical(md5, PENG_MD5))
    halt(sec, "md5 mismatch: expected ", PENG_MD5, ", observed ", md5,
         ". This is not the registered Besca release.")
  message("  ok  A.a/Peng        md5 ", PENG_MD5, " verified")

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
  # Amendment 11: compartments come from celltype1 throughout -- not celltype0
  # with a haematopoietic patch.  The bias claim is verified against the real
  # labels first; if the epithelial set is not identical across the two levels,
  # the amendment's reasoning is falsified and this halts.
  verify_amendment11_bias(ct0, ct1)
  comp <- map_labels(ct1, MAP_Peng, "Peng/celltype1")

  X <- read_h5ad_counts(fp, sec)
  X <- X[, keep, drop = FALSE]

  verify_tissue("Peng", "pancreatic", rownames(X), ct1)

  cells <- data.frame(
    cell_id     = colnames(X),
    patient_id  = pat[keep],
    compartment = comp[keep],
    stringsAsFactors = FALSE
  )
  assert_cells_aligned(X, cells, sec)
  list(X = X, cells = cells, atlas = "Peng", tissue = "pancreatic")
}

# h5ad X -> genes x cells raw counts.  AnnData stores cells x genes CSR; this
# returns the transpose without normalising anything.
read_h5ad_counts <- function(path, sec) {
  ls_h5 <- rhdf5::h5ls(path, recursive = TRUE)

  # `raw/X` is REQUIRED, not preferred.  In a Besca-processed h5ad, `X` is
  # log-normalised; the previous version fell back to it silently whenever /raw
  # was absent or empty, which would pseudobulk normalised values -- exactly what
  # A.d forbids.  A missing raw layer means this is not the registered release,
  # so it halts rather than substituting a different quantity.
  src <- "raw/X"
  if (!length(ls_h5$name[ls_h5$group == "/raw"]))
    halt(sec, "h5ad has no /raw group, so no raw count layer. A.d requires RAW UMI counts ",
         "and `X` in this release is log-normalised. Falling back to `X` would measure the ",
         "wrong quantity. Verify the file is the registered Besca release (md5 ", PENG_MD5, ").")
  fields <- ls_h5$name[ls_h5$group == paste0("/", src)]
  if (!all(c("data", "indices", "indptr") %in% fields))
    halt(sec, "h5ad '", src, "' is not sparse CSR/CSC; found: ", paste(fields, collapse = ", "))

  # AnnData stores obs x var CSR.  Read the declared encoding rather than assuming
  # it: a CSC file read as CSR transposes the meaning of every index.
  enc <- tryCatch(as.character(rhdf5::h5readAttributes(path, src)[["encoding-type"]]),
                  error = function(e) character(0))
  if (length(enc) && !identical(enc, "csr_matrix"))
    halt(sec, "h5ad '", src, "' declares encoding-type '", enc,
         "', not 'csr_matrix'. The reader below assumes CSR (obs x var); a CSC ",
         "matrix read as CSR would transpose genes and cells. Handle explicitly.")
  vg <- "raw/var"
  gsym <- as.character(rhdf5::h5read(path, paste0(vg, "/index")))
  bc   <- as.character(rhdf5::h5read(path, "obs/index"))
  # index1 = FALSE declares `j` 0-based; the CSR `indices` field IS 0-based and is
  # passed through unmodified.  See the note in read_h5_counts_subset: adding 1L here
  # corrects the same offset twice and silently shifts every gene by one column
  # pre-transpose.
  idx0 <- as.integer(rhdf5::h5read(path, paste0(src, "/indices")))
  if (length(idx0) && (min(idx0) < 0L || max(idx0) > length(gsym) - 1L))
    halt(sec, "CSR column indices are not 0-based within [0, n_var-1]: observed range [",
         min(idx0), ", ", max(idx0), "] against n_var = ", length(gsym),
         ". Do not add an offset here; establish the file's indexing convention first.")
  M <- Matrix::sparseMatrix(j = idx0,
                            p = as.integer(rhdf5::h5read(path, paste0(src, "/indptr"))),
                            x = as.numeric(rhdf5::h5read(path, paste0(src, "/data"))),
                            dims = c(length(bc), length(gsym)), index1 = FALSE)
  M <- Matrix::t(M)                 # -> genes x cells
  rownames(M) <- gsym
  colnames(M) <- bc

  # ---- AMENDMENT 13: recover raw counts by inverting the deposited normalisation
  # raw/X is log1p(CP10K), not counts -- the deposit has no count layer at all.
  # counts = round(expm1(raw/X) * n_counts / 1e4).  The three assertions below are
  # the amendment's own, computed over ALL non-zero entries, each halting on
  # failure. They establish EXACT recovery, not approximation.
  n_counts <- as.numeric(rhdf5::h5read(path, "obs/n_counts"))
  if (length(n_counts) != ncol(M))
    halt(sec, "obs/n_counts has ", length(n_counts), " entries for ", ncol(M), " cells")
  if (any(!is.finite(n_counts)) || any(n_counts <= 0))
    halt(sec, "obs/n_counts contains non-positive or non-finite values")

  CP <- M
  CP@x <- expm1(CP@x)                                   # log1p -> CP10K

  # (i) per-cell sum of expm1(raw/X) == 10000, confirming CP10K
  cp_sums <- Matrix::colSums(CP)
  rel_i <- max(abs(cp_sums - 1e4) / 1e4)
  if (!(rel_i < 1e-6))
    halt(sec, "Amendment 13 assertion (i) FAILED: per-cell sum of expm1(raw/X) is not ",
         "10000 within 1e-6 relative tolerance (max relative deviation ",
         format(rel_i, digits = 4), ", worst cell ", colnames(M)[which.max(abs(cp_sums - 1e4))],
         "). raw/X is not CP10K-normalised as Amendment 13 states; the recovery is invalid.")

  R <- CP
  R@x <- CP@x * rep(n_counts, diff(CP@p)) / 1e4          # -> counts, pre-rounding

  # (iii) every recovered value within 0.1 of an integer BEFORE rounding
  dev_iii <- max(abs(R@x - round(R@x)))
  if (!(dev_iii < 0.1))
    halt(sec, "Amendment 13 assertion (iii) FAILED: a recovered value deviates ",
         format(dev_iii, digits = 4), " from the nearest integer, exceeding the 0.1 ",
         "tolerance. Rounding could return a wrong integer, so recovery is not exact ",
         "and the matrix must not be used.")

  # F03: assert_raw_counts() below runs AFTER the rounding, so on this matrix its
  # integrality test is a tautology and cannot fail.  The load-bearing check is
  # assertion (iii) above, on the PRE-rounded values -- that is what establishes
  # round() returns the original integer rather than manufacturing one.
  #
  # NO SECOND, TIGHTER BOUND IS IMPOSED HERE.  An earlier revision added one at
  # 1e-2, citing "observed max 2.9e-3" -- but that figure came from a 400-cell
  # sample taken before the recovery was implemented.  The full-corpus maximum
  # over all 139,415,620 entries is 1.17e-2, which that bound would have FAILED,
  # halting a pipeline that had already passed Amendment 13's own criterion.
  # Amendment 13 sets the tolerance at 0.1 and is the registered authority; the
  # measured 1.17e-2 sits an order of magnitude inside it and far below the 0.5
  # at which rounding could recover a wrong integer.  Tightening a registered
  # tolerance after seeing the data is exactly what preregistration forbids.
  R@x <- round(R@x)

  # (ii) per-cell sum of recovered counts == obs/n_counts
  rec_sums <- Matrix::colSums(R)
  rel_ii <- max(abs(rec_sums - n_counts) / n_counts)
  if (!(rel_ii < 1e-6))
    halt(sec, "Amendment 13 assertion (ii) FAILED: per-cell sum of recovered counts does ",
         "not equal obs/n_counts within 1e-6 relative tolerance (max relative deviation ",
         format(rel_ii, digits = 4), ", worst cell ",
         colnames(M)[which.max(abs(rec_sums - n_counts) / n_counts)],
         "). n_counts is a library size over a different gene set than raw/X; the ",
         "recovery is invalid.")

  message(sprintf("  ok  %-28s Amendment 13 recovery over %s non-zero entries:", sec,
                  format(length(R@x), big.mark = ",")))
  message(sprintf("      (i)   max rel. deviation of per-cell sum(expm1) from 10000 = %.3e  (tol 1e-6)", rel_i))
  message(sprintf("      (ii)  max rel. deviation of recovered per-cell sum from n_counts = %.3e  (tol 1e-6)", rel_ii))
  message(sprintf("      (iii) max deviation from integer before rounding = %.3e  (tol 0.1)", dev_iii))

  # A.d's integrality check is NOT relaxed: it runs unchanged on the recovered
  # matrix and must pass on genuine integers.
  assert_raw_counts(R, sec, paste0(src, " (Amendment 13 recovered counts)"))
  R
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

pseudobulk_raw <- function(X, cells, genes, detection = FALSE) {
  S <- matrix(0, nrow = length(genes), ncol = length(COMPARTMENTS),
              dimnames = list(genes, COMPARTMENTS))
  # F17: fraction of cells in each compartment with >=1 count for each gene.
  # Prespecification exclusion rule 2 is defined on this quantity, and B.h assigns
  # its computation to Part A.  Computed only when asked: inside the 2000-resample
  # bootstrap it is never needed, and the `> 0` pass is not free.
  D <- if (detection) matrix(NA_real_, nrow = length(genes), ncol = length(COMPARTMENTS),
                             dimnames = list(genes, COMPARTMENTS)) else NULL
  n_cells <- setNames(integer(length(COMPARTMENTS)), COMPARTMENTS)
  have <- intersect(genes, rownames(X))

  # Row selection by NAME returns only the FIRST match.  CellRanger feature lists
  # repeat gene symbols (distinct Ensembl IDs collapsing to one symbol), so a
  # duplicated panel symbol would be silently undercounted by whatever its second
  # row carries -- a per-gene, per-atlas bias in f(pi) that no downstream check
  # can see.  Halt rather than pick a summing rule here: aggregating duplicate
  # rows is a quantitative decision and belongs in an amendment, not in a helper.
  dup_panel <- unique(rownames(X)[duplicated(rownames(X)) & rownames(X) %in% genes])
  if (length(dup_panel))
    halt("A.d", "panel gene symbol(s) appear on more than one matrix row: ",
         paste(utils::head(sort(dup_panel), 10), collapse = ", "),
         if (length(dup_panel) > 10) paste0(" (+", length(dup_panel) - 10, " more)") else "",
         ". Name-based row selection would take only the first and undercount the rest. ",
         "Decide the aggregation rule as a dated amendment; do not silently sum here.")

  for (cc in COMPARTMENTS) {
    idx <- which(cells$compartment == cc)
    n_cells[cc] <- length(idx)
    if (!length(idx) || !length(have)) next
    S[have, cc] <- as.numeric(Matrix::rowSums(X[have, idx, drop = FALSE]))  # RAW
    if (detection)
      D[have, cc] <- as.numeric(Matrix::rowSums(X[have, idx, drop = FALSE] > 0)) / length(idx)
  }
  # Genes absent from this atlas are NA, never 0 (A.d).
  S[setdiff(genes, have), ] <- NA_real_
  if (detection) D[setdiff(genes, have), ] <- NA_real_
  list(counts = S, n_cells = n_cells, detection = D)
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

one_atlas_stats <- function(X, cells, genes, detection = FALSE) {
  pb <- pseudobulk_raw(X, cells, genes, detection = detection)
  tot <- rowSums(pb$counts, na.rm = TRUE)
  # by name, not position: the previous `pb$counts[, 1]` silently depended on
  # COMPARTMENTS[1] being "epithelial"
  evidence_ok <- !is.na(pb$counts[, "epithelial"]) & tot >= EVIDENCE_MIN
  F  <- sweep_f(pb$counts, pb$n_cells)
  rownames(F) <- genes
  # F17: the RAW compartment share A.d defines, S / rowSums(S).  Distinct from
  # f(pi), which reweights compartments by an assumed purity; this is the
  # unweighted share of each gene's pseudobulk signal contributed by each
  # compartment, and it is what prespecification exclusion rule 2 refers to.
  frac <- pb$counts / rowSums(pb$counts, na.rm = FALSE)
  list(pb = pb, F = F, evidence_ok = evidence_ok, frac = frac,
       dom = dominance_from_F(F, evidence_ok),
       dom50 = { d <- F[, which.min(abs(PI_GRID - 0.50))] > 0.50
                 d[!evidence_ok] <- NA; d })
}

bootstrap_atlas <- function(X, cells, genes, atlas_index, B = B_RESAMPLES) {
  pats <- unique(cells$patient_id)
  by_pat <- split(seq_len(nrow(cells)), cells$patient_id)

  # R10 performance change.  pseudobulk_raw only ever reads rows in `genes`, so
  # restricting X to those rows once here is numerically identical to carrying the
  # full matrix through every resample.  It exists because the full matrix would
  # otherwise be duplicated on each of 2000 iterations.
  #
  # Verified on real GSE125449 data (2026-08-01): 20 paired resamples under the
  # same seed, identical() TRUE across the entire nested output -- f30, dom,
  # dom50, F, evidence_ok, pseudobulk counts and n_cells -- at 47.5x the speed.
  # The benchmark ran as a standalone comparison script, not as part of the Part A
  # run, so output/partA_run.log does NOT contain it; the result is recorded in
  # NOTES_FOR_REVIEW.md section 14 and in the commit that introduced this line.
  X <- X[intersect(rownames(X), genes), , drop = FALSE]
  withr::with_seed(SEED_BASE + atlas_index, {
    reps <- lapply(seq_len(B), function(b) {
      draw <- sample(pats, length(pats), replace = TRUE)
      idx  <- unlist(by_pat[draw], use.names = FALSE)
      st   <- one_atlas_stats(X[, idx, drop = FALSE], cells[idx, , drop = FALSE], genes)
      list(f30 = st$F[, 1], dom = st$dom, dom50 = st$dom50)
    })
  })
  f30_mat    <- vapply(reps, `[[`, numeric(length(genes)), "f30")
  dom_reps   <- vapply(reps, `[[`, logical(length(genes)), "dom")
  dom50_reps <- vapply(reps, `[[`, logical(length(genes)), "dom50")
  rownames(f30_mat) <- rownames(dom_reps) <- rownames(dom50_reps) <- genes

  # F06: per-gene per-atlas dominance RATE -- the proportion of resamples in which
  # the gene was epithelial-dominant, with its percentile interval.  The
  # replicates were already being computed and stored; previously they were
  # consumed only for the aggregate k intervals, so no per-gene stability measure
  # reached the output at all.
  dom_rate   <- rowMeans(dom_reps,   na.rm = TRUE)
  dom50_rate <- rowMeans(dom50_reps, na.rm = TRUE)
  ci <- function(M) t(apply(M, 1, function(r)
    if (all(is.na(r))) c(NA_real_, NA_real_) else quantile(r, c(0.025, 0.975), na.rm = TRUE)))
  dr_ci <- ci(dom_reps)

  list(
    f30_lo = apply(f30_mat, 1, quantile, 0.025, na.rm = TRUE),
    f30_hi = apply(f30_mat, 1, quantile, 0.975, na.rm = TRUE),
    dom_rate      = dom_rate,
    dom50_rate    = dom50_rate,
    dom_rate_lo   = dr_ci[, 1],
    dom_rate_hi   = dr_ci[, 2],
    # G04: the DENOMINATOR behind each rate.  evidence_ok is recomputed inside
    # every resample, so a gene near the 20-count floor is evaluable in some draws
    # and NA in others; rowMeans(na.rm = TRUE) then divides by a per-gene count
    # that varies.  Reporting the rate without its denominator would let a gene
    # dominant in 50 of 50 evaluable draws read identically to one dominant in
    # 2000 of 2000.
    dom_n_eval   = rowSums(!is.na(dom_reps)),
    dom50_n_eval = rowSums(!is.na(dom50_reps)),
    # G02: retain the per-resample f(0.30) matrix.  Without it the saved object
    # cannot support a later interval recomputation, which is the stated purpose
    # of persisting the replicates at all.
    f30_reps   = f30_mat,
    dom_reps   = dom_reps,
    dom50_reps = dom50_reps
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
  # S1: every k variant is computed over exactly the 152 locked panel genes.
  assert_inferential_set(dominance, "dominance matrix passed to compute_k")
  n_dom  <- rowSums(dominance, na.rm = TRUE)
  n_eval <- rowSums(!is.na(dominance))
  k         <- sum(n_dom >= 2L)
  k_all3    <- sum(n_dom == 3L & n_eval == 3L)
  k_evalall <- sum(n_dom >= 2L & n_eval == 3L)
  stopifnot(k_all3 <= k_evalall, k_evalall <= k)
  list(k = k, k_all3 = k_all3, k_evalall = k_evalall,
       n_dom = n_dom, n_eval = n_eval)
}

# Amendment 3 branch.  Wording is Amendment 3's own, not a paraphrase: the middle
# band is EXPLORATORY ONLY, which is materially weaker than "reported with
# caveats".
#
# Both bounds are PROPORTIONS of the final panel, so they are computed from
# n_panel rather than written as literals:
#   BUILT       iff  k/n >= 20%  AND  k >= 8   (the amendment's absolute floor)
#   EXPLORATORY iff  5% <= k/n < 20%
#   DESCRIPTIVE iff  k/n < 5%
# For integer k, `k/n >= p` is exactly `k >= ceiling(p*n)` -- including when p*n
# is itself an integer -- so ceiling() implements the proportion faithfully at
# any n, not merely at n = 152.  verify_branch_bands() below asserts this rather
# than leaving it as a claim in a comment.
#
# NB the two lower bounds are DISTINCT quantities that happen to coincide at
# n = 152: ceiling(0.05 * 152) = 8 and the amendment's absolute floor is also 8.
# They diverge elsewhere (n = 140 -> 7, n = 161 -> 9).  The assertion pins the
# coincidence so a future panel-size change cannot silently conflate them.
BRANCH_ABS_FLOOR <- 8L   # Amendment 3's absolute floor on the top band

branch_of_k <- function(k, n_panel = length(PANEL_GENES)) {
  hi <- ceiling(0.20 * n_panel)
  lo <- ceiling(0.05 * n_panel)
  if (k >= hi && k >= BRANCH_ABS_FLOOR)
    "epithelial subscore BUILT; subscore survival models run as SECONDARY"
  else if (k >= lo)
    "subscore reported as EXPLORATORY only"
  else
    "decomposition DESCRIPTIVE only"
}

# verify_branch_bands: the bands are checked against Amendment 3 stated as
# PROPORTIONS, over every k for a range of panel sizes.  This is what makes the
# construction correct by definition rather than by arithmetic luck at n = 152.
verify_branch_bands <- function(n_panel) {
  band <- function(s) if (grepl("BUILT", s)) "BUILT"
                      else if (grepl("EXPLORATORY", s)) "EXPLORATORY" else "DESCRIPTIVE"
  for (n in unique(c(n_panel, 20L, 140L, 152L, 160L, 161L, 200L, 400L))) {
    for (k in 0:n) {
      want <- if (k >= 0.20 * n && k >= BRANCH_ABS_FLOOR) "BUILT"
              else if (k >= 0.05 * n) "EXPLORATORY" else "DESCRIPTIVE"
      got  <- band(branch_of_k(k, n))
      if (!identical(want, got))
        halt("A.g", "Amendment 3 band mismatch at n_panel = ", n, ", k = ", k,
             ": proportion rule gives ", want, ", branch_of_k gives ", got, ".")
    }
  }
  # The two lower bounds are separate constructions; record their relationship at
  # the live panel size instead of assuming they always agree.
  lo <- ceiling(0.05 * n_panel)
  message(sprintf("  ok  A.g  Amendment 3 bands verified as proportions (n=%d: DESCRIPTIVE k<%d, EXPLORATORY %d-%d, BUILT k>=%d%s)",
                  n_panel, lo, lo, ceiling(0.20 * n_panel) - 1L, ceiling(0.20 * n_panel),
                  if (lo == BRANCH_ABS_FLOOR)
                    sprintf("; 5%% bound and absolute floor coincide at %d", lo)
                  else
                    sprintf("; 5%% bound %d differs from absolute floor %d", lo, BRANCH_ABS_FLOOR)))
  invisible(TRUE)
}

# ==============================================================================
# DRIVER
# ==============================================================================
# Pre-data invariant, same class as verify_k_ordering(): runs before any atlas
# file is opened, so a band-definition error is caught without touching data.
verify_branch_bands(length(PANEL_GENES))

message("\n== A.a/A.b  loading the three registered atlases ==")
LOADERS <- list(GSE125449 = load_GSE125449,
                GSE178341 = load_GSE178341,
                Peng      = load_Peng)
if (!setequal(names(LOADERS), atlases$atlas))
  halt("A.a", "loader set does not match the registered atlas set")

dat <- lapply(atlases$atlas, function(a) LOADERS[[a]]())
names(dat) <- atlases$atlas

message("\n== A.d/A.e  pseudobulk and purity sweep ==")
# Statistics are computed over the 155-gene REPORTING set so the three
# non-qualifying origin genes get compartment fractions (S1).  Everything
# inferential is restricted back to the 152-gene panel immediately below.
stats <- lapply(names(dat), function(a) {
  st <- one_atlas_stats(dat[[a]]$X, dat[[a]]$cells, REPORT_GENES, detection = TRUE)
  assert_monotone(st$F, a)
  message(sprintf("  ..  %-12s evaluable genes = %d / %d panel (+%d / %d non-qualifying)",
                  a, sum(st$evidence_ok[PANEL_GENES]), length(PANEL_GENES),
                  sum(st$evidence_ok[ORIGIN_NONQUAL]), length(ORIGIN_NONQUAL)))
  st
})
names(stats) <- names(dat)

# THE INFERENTIAL FIREWALL (S1).  Dominance is subset to the locked 152 BEFORE
# any k variant is computed; the three non-qualifying genes stop here.
dominance    <- vapply(stats, function(s) s$dom[PANEL_GENES],   logical(length(PANEL_GENES)))
dominance_50 <- vapply(stats, function(s) s$dom50[PANEL_GENES], logical(length(PANEL_GENES)))
rownames(dominance) <- rownames(dominance_50) <- PANEL_GENES
colnames(dominance) <- colnames(dominance_50) <- names(stats)
assert_inferential_set(dominance,    "dominance matrix")
assert_inferential_set(dominance_50, "dominance matrix at pi = 0.50")
if (any(ORIGIN_NONQUAL %in% rownames(dominance)))
  halt("A.g", "non-qualifying gene(s) reached the dominance matrix: ",
       paste(intersect(ORIGIN_NONQUAL, rownames(dominance)), collapse = ", "))

message("\n== A.g  k and variants ==")
# Reproducibility gate.  The guard repairs applied after the first complete run
# (F01-F04, F18, plus F12-F16) were all intended to change NO number: they close
# latent hazards, replace checks that could never fire, and make skippable
# verifications unconditional.  This asserts that intent against the realised
# values rather than assuming it.  A mismatch means a repair altered a result and
# must be reported before anything downstream is read.
K_EXPECTED <- c(k = 46L, k_all3 = 24L, k_evalall = 43L, k_50 = 96L)
K   <- compute_k(dominance)
K50 <- compute_k(dominance_50)
message(sprintf("  k = %d | k_all3 = %d | k_evalall = %d | k_50 = %d",
                K$k, K$k_all3, K$k_evalall, K50$k))
local({
  obs <- c(k = K$k, k_all3 = K$k_all3, k_evalall = K$k_evalall, k_50 = K50$k)
  if (!identical(as.integer(obs), as.integer(K_EXPECTED)))
    halt("A.g", "POINT ESTIMATES CHANGED after the post-audit guard repairs. ",
         "observed ", paste(names(obs), obs, sep = "=", collapse = " "),
         " vs expected ", paste(names(K_EXPECTED), K_EXPECTED, sep = "=", collapse = " "),
         ". Every repair in that set was intended to be number-neutral, so a ",
         "difference is a finding: report it before reading anything downstream. ",
         "Do not update K_EXPECTED to match.")

  # G06: the four counts above are aggregates -- two genes could swap dominance
  # and leave every one of them unchanged.  The claim being asserted is that no
  # NUMBER changed, so compare the underlying matrices, not their sums.  Recorded
  # as a per-gene-per-atlas signature over both dominance matrices.
  sig <- paste(sum(dominance, na.rm = TRUE), sum(!is.na(dominance)),
               paste(which(!is.na(dominance) & dominance), collapse = ","),
               sum(dominance_50, na.rm = TRUE), sum(!is.na(dominance_50)),
               paste(which(!is.na(dominance_50) & dominance_50), collapse = ","),
               sep = "|")
  sig_file <- file.path(OUTDIR, "dominance_signature.txt")
  if (file.exists(sig_file)) {
    prev <- readLines(sig_file, warn = FALSE)[1]
    if (!identical(prev, sig))
      halt("A.g", "The DOMINANCE MATRIX changed even though k, k_all3, k_evalall ",
           "and k_50 are unchanged -- individual gene calls have moved while the ",
           "aggregates coincide. This is exactly what the aggregate gate cannot ",
           "see. Compare output/dominance_signature.txt against the previous run ",
           "and report before proceeding. Do not delete the signature file.")
    message("  ok  A.g  dominance matrix identical to the previous run (per-gene signature match)")
  } else {
    message("  ..  A.g  no previous dominance signature on disk; recording one for future runs")
  }
  writeLines(sig, sig_file)
  message("  ok  A.g  point estimates reproduce the pre-repair run exactly (k=46, k_all3=24, k_evalall=43, k_50=96)")
})
message("  Amendment 3 branch (on primary k): ", branch_of_k(K$k))

message("\n== A.f  patient-level bootstrap (B = ", B_RESAMPLES, ", unit = patient) ==")
boots <- lapply(seq_along(dat), function(i)
  bootstrap_atlas(dat[[i]]$X, dat[[i]]$cells, REPORT_GENES, atlas_index = i))
names(boots) <- names(dat)

# Persist the resample objects.  The bootstrap is the expensive step (~2 h at
# B=2000 over three atlases); saving the replicate matrices means a later
# reporting addition can be derived without re-running it.  Seed and B are stored
# alongside so a reload can be checked against the run that produced it.
# G02: the object must carry everything a later reporting addition needs, or the
# 2 h it exists to avoid gets spent anyway.  boots[[a]] now includes f30_reps
# (per-resample f(0.30)), dom_reps and dom50_reps, so any interval or rate can be
# recomputed; the scalars below let a reload be checked against the run that
# produced it rather than assumed compatible.
saveRDS(list(boots        = boots,
             seed_base    = SEED_BASE,
             B            = B_RESAMPLES,
             genes        = REPORT_GENES,
             panel        = PANEL_GENES,
             atlases      = names(dat),
             pi_grid      = PI_GRID,
             evidence_min = EVIDENCE_MIN,
             n_patients   = vapply(dat, function(d) length(unique(d$cells$patient_id)), integer(1)),
             script_md5   = unname(tools::md5sum("03_compartments.R")),
             run_date     = Sys.Date(),
             R_version    = R.version.string),
        file.path(OUTDIR, "bootstrap_replicates.rds"))
message("  ..  saved output/bootstrap_replicates.rds (B = ", B_RESAMPLES,
        ", seed_base = ", SEED_BASE, ")")

# Amendment 6 requires k_50 to be reported WITH its bootstrap interval, on the
# same footing as k.  bootstrap_atlas already returns dom50_reps -- the pi = 0.50
# dominance calls for every resample -- so k_50's interval is derived from the
# SAME 2000 patient-level resamples as the others.  A previous version consumed
# only dom_reps here and wrote NA into k_50's interval columns, which left an
# Amendment 6 requirement unmet even though the replicates existed.
k_reps <- vapply(seq_len(B_RESAMPLES), function(b) {
  # rows are the 155 reporting genes; subset to the locked 152 before compute_k
  dm   <- vapply(boots, function(bt) bt$dom_reps[PANEL_GENES,   b], logical(length(PANEL_GENES)))
  dm50 <- vapply(boots, function(bt) bt$dom50_reps[PANEL_GENES, b], logical(length(PANEL_GENES)))
  rownames(dm) <- rownames(dm50) <- PANEL_GENES
  c(unlist(compute_k(dm)[c("k", "k_all3", "k_evalall")]),
    k_50 = compute_k(dm50)$k)
}, numeric(4))
if (!identical(rownames(k_reps), c("k", "k_all3", "k_evalall", "k_50")))
  halt("A.f", "k_reps rows are ", paste(rownames(k_reps), collapse = ", "),
       "; expected k, k_all3, k_evalall, k_50. The CI columns below index by name.")
k_ci <- apply(k_reps, 1, quantile, c(0.025, 0.975), na.rm = TRUE)

# ==============================================================================
# OUTPUTS  (A.g required reporting)
# ==============================================================================
# Reported over all 155 genes (S1): the three non-qualifying origin genes carry
# compartment fractions and are labelled, never omitted.  `in_panel` marks which
# rows were eligible for k.
dom_long <- do.call(rbind, lapply(names(stats), function(a) data.frame(
  gene        = REPORT_GENES,
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
  # F06: per-gene per-atlas bootstrap stability of the dominance call
  dom_rate    = boots[[a]]$dom_rate,
  dom_rate_lo = boots[[a]]$dom_rate_lo,
  dom_rate_hi = boots[[a]]$dom_rate_hi,
  dom_n_eval  = boots[[a]]$dom_n_eval,      # G04: denominator behind dom_rate
  dom50_rate  = boots[[a]]$dom50_rate,
  origin_six  = REPORT_GENES %in% ORIGIN_SIX,
  in_panel    = REPORT_GENES %in% PANEL_GENES,
  # F05: values for the non-qualifying three are RETAINED, not blanked.
  # Amendment 2's disclosure requires readers to see whether excluding MMP9 and
  # HGF was consequential; suppressing their fractions would destroy exactly the
  # evidence the amendment exists to expose.  `qualifying` marks which rows were
  # eligible for a k variant.
  qualifying  = REPORT_GENES %in% PANEL_GENES,
  stringsAsFactors = FALSE)))
assert_n(sum(!dom_long$in_panel), 3L * length(stats), "S1",
         "non-qualifying rows present in the reported table")

write.csv(dom_long, file.path(OUTDIR, "compartment_dominance_matrix.csv"), row.names = FALSE)

# Header note travels with the file rather than living only in a commit message.
writeLines(c(
  "# compartment_dominance_matrix.csv -- readme",
  "#",
  "# 155 reporting genes x 3 atlases = 465 rows.",
  "#",
  "# qualifying / in_panel = FALSE marks BCL2, MMP9 and HGF: the labelled",
  "# non-qualifying subset (prespecification section 4, Amendment 2). Their",
  "# compartment fractions and dominance calls ARE reported, deliberately -- they",
  "# let a reader judge whether excluding MMP9 and HGF was consequential. They",
  "# enter NO k variant: k, k_50, k_all3 and k_evalall are computed over exactly",
  "# the 152 locked panel genes, enforced by assert_inferential_set().",
  "#",
  "# dominant     : f(pi) > 0.50 at EVERY point of the 41-point grid pi in [0.30, 0.70]",
  "# dominant_50  : f(0.50) > 0.50 only (Amendment 6)",
  "# evaluable    : >= 20 summed counts in that atlas (below this, no fraction is given)",
  "#",
  "# dom_rate     : proportion of the 2000 patient-level resamples in which the gene",
  "#                was epithelial-dominant. A stability measure, not a probability.",
  "# dom_n_eval   : the DENOMINATOR behind dom_rate -- how many of the 2000 resamples",
  "#                the gene was evaluable in. evidence_ok is recomputed inside every",
  "#                resample, so a gene near the 20-count floor is evaluable in only",
  "#                some draws. Read dom_rate WITH this column: 1.00 from 50 draws is",
  "#                not the same evidence as 1.00 from 2000.",
  "# dom_rate_lo/hi: A.f's percentile bootstrap CI (2.5/97.5) on the dominance",
  "#                INDICATOR. This is the registered statistic and is reported as",
  "#                specified, but note it is degenerate by construction: the",
  "#                indicator is 0/1, so for any gene dominant in more than 97.5% or",
  "#                fewer than 2.5% of resamples the interval collapses to [1,1] or",
  "#                [0,0]. dom_rate and dom_n_eval carry the usable information.",
  "# dom50_rate   : as dom_rate, at pi = 0.50 only (Amendment 6)",
  "#",
  "# ENDOTHELIAL COMPARTMENT, GSE178341: that atlas's clTopLevel vocabulary is",
  "# B | Epi | Mast | Myeloid | Plasma | Strom | TNKILC -- it has NO endothelial",
  "# label, and its endothelial cells sit inside Strom. Amendment 4 maps each",
  "# atlas's OWN level-1 labels, so the code matches the registration and no",
  "# amendment is required. The estimand is unaffected: subdividing non-epithelial",
  "# mass changes neither the numerator nor the denominator of f(pi). It is a",
  "# REPORTING limitation only -- per-compartment breakdowns are not comparable",
  "# across atlases for the endothelial and fibroblast_stromal rows."
), file.path(OUTDIR, "compartment_dominance_matrix.readme.txt"))

# --- F17: compartment fraction S/rowSums(S), and per-compartment detection rate
# Both over the 155 reporting genes. Prespecification exclusion rule 2 is defined
# on the detection rate, and B.h assigns its computation to Part A.
frac_long <- do.call(rbind, lapply(names(stats), function(a) {
  fr <- stats[[a]]$frac
  data.frame(gene = rep(rownames(fr), ncol(fr)),
             atlas = a,
             tissue = atlases$tissue[atlases$atlas == a],
             compartment = rep(colnames(fr), each = nrow(fr)),
             fraction = as.numeric(fr),
             n_cells = rep(stats[[a]]$pb$n_cells[colnames(fr)], each = nrow(fr)),
             qualifying = rep(rownames(fr) %in% PANEL_GENES, ncol(fr)),
             stringsAsFactors = FALSE)
}))
write.csv(frac_long, file.path(OUTDIR, "compartment_fractions.csv"), row.names = FALSE)

det_long <- do.call(rbind, lapply(names(stats), function(a) {
  D <- stats[[a]]$pb$detection
  if (is.null(D)) halt("A.d", a, ": detection rate was not computed")
  data.frame(gene = rep(rownames(D), ncol(D)),
             atlas = a,
             tissue = atlases$tissue[atlases$atlas == a],
             compartment = rep(colnames(D), each = nrow(D)),
             pct_cells_detected = as.numeric(D),
             n_cells = rep(stats[[a]]$pb$n_cells[colnames(D)], each = nrow(D)),
             qualifying = rep(rownames(D) %in% PANEL_GENES, ncol(D)),
             stringsAsFactors = FALSE)
}))
write.csv(det_long, file.path(OUTDIR, "detection_rate_by_compartment.csv"), row.names = FALSE)
message(sprintf("  ..  wrote compartment_fractions.csv (%d rows) and detection_rate_by_compartment.csv (%d rows)",
                nrow(frac_long), nrow(det_long)))

# G07: the endothelial caveat governs THESE two files most of all -- they are the
# per-compartment breakdown -- so it travels with them rather than only with the
# dominance matrix.
writeLines(c(
  "# compartment_fractions.csv and detection_rate_by_compartment.csv -- readme",
  "#",
  "# Both cover the 155 reporting genes x 3 atlases x 6 compartments (2,790 rows).",
  "# qualifying = FALSE marks BCL2, MMP9, HGF: reported, but in no k variant.",
  "#",
  "# fraction           : S[g,c] / sum_c S[g,c] -- the RAW share of a gene's",
  "#                      pseudobulk signal contributed by compartment c. This is",
  "#                      NOT f(pi): no purity reweighting is applied. A.d defines",
  "#                      this quantity and prespecification exclusion rule 2 is",
  "#                      stated on it.",
  "# pct_cells_detected : fraction of cells in compartment c with >= 1 count for",
  "#                      gene g. Independent of expression magnitude.",
  "# n_cells            : cells in that compartment in that atlas.",
  "#",
  "# COMPARTMENTS ABSENT FROM AN ATLAS. A fraction of 0 with n_cells = 0 means the",
  "# compartment DOES NOT EXIST in that atlas's annotation -- not that it",
  "# contributes nothing. GSE178341 is the case that matters: its clTopLevel",
  "# vocabulary (B | Epi | Mast | Myeloid | Plasma | Strom | TNKILC) has no",
  "# endothelial label, so its endothelial cells sit inside Strom and are counted",
  "# as fibroblast_stromal. Amendment 4 maps each atlas's OWN level-1 labels, so",
  "# this matches the registration and no amendment is required; the estimand is",
  "# unaffected because subdividing non-epithelial mass changes neither the",
  "# numerator nor the denominator of f(pi).",
  "#",
  "# CONSEQUENCE: the endothelial and fibroblast_stromal rows are NOT comparable",
  "# across atlases. 'Endothelial contribution is lower in colorectal' would be an",
  "# artefact of annotation granularity. Epithelial, myeloid and lymphoid rows are",
  "# comparable; 'other' is atlas-specific by construction."
), file.path(OUTDIR, "compartment_fractions.readme.txt"))

evaluability <- table(factor(K$n_eval, levels = 3:0))
write.csv(data.frame(atlases_evaluable = names(evaluability),
                     n_genes = as.integer(evaluability)),
          file.path(OUTDIR, "evaluability_distribution.csv"), row.names = FALSE)

# panel genes only -- per-tissue dominance is an inferential quantity (S1)
per_tissue <- vapply(names(stats),
                     function(a) sum(stats[[a]]$dom[PANEL_GENES], na.rm = TRUE), integer(1))
write.csv(data.frame(atlas = names(per_tissue),
                     tissue = atlases$tissue[match(names(per_tissue), atlases$atlas)],
                     n_dominant = as.integer(per_tissue)),
          file.path(OUTDIR, "per_tissue_dominance.csv"), row.names = FALSE)

k_tab <- data.frame(
  quantity = c("k", "k_all3", "k_evalall", "k_50"),
  estimate = c(K$k, K$k_all3, K$k_evalall, K50$k),
  ci_lo    = c(k_ci[1, "k"], k_ci[1, "k_all3"], k_ci[1, "k_evalall"], k_ci[1, "k_50"]),
  ci_hi    = c(k_ci[2, "k"], k_ci[2, "k_all3"], k_ci[2, "k_evalall"], k_ci[2, "k_50"]),
  stringsAsFactors = FALSE)
if (anyNA(k_tab$ci_lo) || anyNA(k_tab$ci_hi))
  halt("A.g", "k_estimates has a missing interval bound for: ",
       paste(k_tab$quantity[is.na(k_tab$ci_lo) | is.na(k_tab$ci_hi)], collapse = ", "),
       ". Amendment 6 requires k_50 to be reported with its bootstrap interval.")

# --- F09: the two prespecified disclosures, COMPUTED rather than left to a human
# reading a CSV.  Both are decided by the numbers, so deferring them to after the
# numbers are seen is exactly the discretion preregistration removes.
branch_lo <- branch_of_k(floor(k_ci[1, "k"]))
branch_hi <- branch_of_k(ceiling(k_ci[2, "k"]))
k_ci_spans_branch <- !identical(branch_lo, branch_hi)
k50_branch_differs <- !identical(branch_of_k(K$k), branch_of_k(K50$k))
k_tab$branch <- c(branch_of_k(K$k), NA, NA, branch_of_k(K50$k))
write.csv(k_tab, file.path(OUTDIR, "k_estimates.csv"), row.names = FALSE)

writeLines(c(
  "# Amendment 6 disclosures, computed from the run (F09)",
  "",
  paste0("k                    : ", K$k, "  (95% CI ", k_ci[1, "k"], " - ", k_ci[2, "k"], ")"),
  paste0("branch at k          : ", branch_of_k(K$k)),
  paste0("branch at CI lower   : ", branch_lo),
  paste0("branch at CI upper   : ", branch_hi),
  paste0("k_ci_spans_branch    : ", k_ci_spans_branch),
  if (k_ci_spans_branch)
    "  -> DISCLOSURE REQUIRED: k's interval crosses an Amendment 3 band boundary; the branch is not robust to sampling variation and must be reported as such."
  else
    "  -> no disclosure required on this ground: the whole interval selects one branch.",
  "",
  paste0("k_50                 : ", K50$k, "  (95% CI ", k_ci[1, "k_50"], " - ", k_ci[2, "k_50"], ")"),
  paste0("branch at k_50       : ", branch_of_k(K50$k)),
  paste0("k50_branch_differs   : ", k50_branch_differs),
  if (k50_branch_differs)
    "  -> DISCLOSURE REQUIRED: the full-band rule and the pi=0.50 rule select DIFFERENT Amendment 3 branches. Amendment 6 records that the full-band rule is the stricter one and leans toward this study's thesis; the divergence must be reported."
  else
    "  -> no disclosure required on this ground: both rules select the same branch."
), file.path(OUTDIR, "amendment6_disclosures.txt"))
message("  ..  Amendment 6: k_ci_spans_branch = ", k_ci_spans_branch,
        " | k50_branch_differs = ", k50_branch_differs)

# Origin six as a labelled subset (prespecification section 4): SOCS3/MYC/IL6
# qualify, BCL2/MMP9/HGF do not.  All six appear, with `qualifying` and a reason;
# a previous version drew from a 152-gene table in which the non-qualifying three
# could not appear, so `qualifying` was TRUE in every row.
origin_tab <- dom_long[dom_long$gene %in% ORIGIN_SIX, ]
origin_tab$non_qualifying_reason <- ifelse(
  origin_tab$gene %in% PANEL_GENES, "",
  ifelse(origin_tab$gene == "BCL2", "no ChIP-seq evidence in any group (fails Criterion B)",
                                    "mouse ChIP-seq only (excluded by Amendment 2, human-only)"))
assert_n(length(unique(origin_tab$gene)), 6L, "S1", "origin-six genes in the reported table")
if (!setequal(unique(origin_tab$gene[!origin_tab$qualifying]), ORIGIN_NONQUAL))
  halt("S1", "origin-six table's non-qualifying rows are ",
       paste(sort(unique(origin_tab$gene[!origin_tab$qualifying])), collapse = ", "),
       "; expected BCL2, MMP9, HGF.")
write.csv(origin_tab, file.path(OUTDIR, "origin_six_compartment.csv"), row.names = FALSE)

# --- F07: Amendment 4's malignant-restricted sensitivity, GSE125449 only
# Amendment 4 changed the estimand from malignant-epithelial to lineage-level
# epithelial, and names a sensitivity analysis restricting epithelial to
# CNV-called malignant cells.  GSE125449 is the ONLY atlas that can support it:
# Ma's `Type` column carries inferCNV-derived calls, so `Malignant cell` and
# `HPC-like` are separable.  Pelka's clTopLevel `Epi` and Peng's celltype1 carry
# no malignant/normal split, so this is a single-atlas sensitivity and cannot be
# replicated -- stated in the output, not left for a reader to infer.
# POINT ESTIMATES ONLY: no bootstrap, so no interval is reported.
local({
  a  <- "GSE125449"
  d  <- dat[[a]]
  hp <- d$cells$compartment == "epithelial" & d$cells$cell_type_raw == "HPC-like"
  if (!any(hp)) halt("A.g/sensitivity", a, ": no HPC-like cells found; the raw label ",
                     "column is missing or has changed, so the sensitivity cannot be computed.")
  cells_m <- d$cells
  cells_m$compartment[hp] <- "other"     # malignant-restricted epithelial
  st_m <- one_atlas_stats(d$X, cells_m, REPORT_GENES)
  st_0 <- stats[[a]]

  out <- data.frame(
    gene              = REPORT_GENES,
    qualifying        = REPORT_GENES %in% PANEL_GENES,
    f30_lineage       = st_0$F[, 1],
    f30_malignant     = st_m$F[, 1],
    delta_f30         = st_m$F[, 1] - st_0$F[, 1],
    dominant_lineage  = st_0$dom,
    dominant_malignant= st_m$dom,
    evaluable_lineage = st_0$evidence_ok,
    evaluable_malig   = st_m$evidence_ok,
    stringsAsFactors  = FALSE)
  out$call_changed <- !is.na(out$dominant_lineage) & !is.na(out$dominant_malignant) &
                       out$dominant_lineage != out$dominant_malignant
  write.csv(out, file.path(OUTDIR, "gse125449_malignant_sensitivity.csv"), row.names = FALSE)

  p <- out[out$qualifying, ]
  writeLines(c(
    "# Amendment 4 sensitivity: malignant-restricted epithelial, GSE125449 only",
    "#",
    "# Primary estimand (Amendment 4) is LINEAGE-level epithelial: Malignant cell +",
    "# HPC-like. This sensitivity restricts epithelial to CNV-called Malignant cell",
    "# and reassigns HPC-like to `other`.",
    "#",
    "# SINGLE-ATLAS LIMITATION: only GSE125449 deposits malignancy calls (Ma,",
    "# inferCNV: score > 80th pct AND corr > 0.4). GSE178341's clTopLevel `Epi` and",
    "# Peng's celltype1 carry no malignant/normal split, so this cannot be",
    "# replicated in the other two atlases and no k variant is recomputed from it.",
    "# It bears on ONE of the three atlases feeding the two-of-three rule.",
    "#",
    "# POINT ESTIMATES ONLY -- no bootstrap, no intervals.",
    "#",
    "# DISPOSITION OF NON-MALIGNANT EPITHELIUM (G05). Amendment 4 fixes the",
    "# numerator -- epithelial becomes CNV-called Malignant cell only -- but does",
    "# not say what becomes of the non-malignant epithelium. Two implementations",
    "# satisfy its wording and give different denominators: reassign HPC-like to a",
    "# non-epithelial compartment (done here: they move to `other`, so they stay in",
    "# the denominator), or drop those cells entirely (which would shrink it).",
    "# Retention is the conservative choice for this study's thesis: it keeps",
    "# non-malignant transcripts in the denominator and so cannot inflate the",
    "# epithelial fraction. Recorded because the amendment is silent, not because",
    "# the choice is forced.",
    "#",
    paste0("HPC-like cells reassigned : ", sum(hp)),
    paste0("epithelial cells, lineage : ", sum(d$cells$compartment == "epithelial")),
    paste0("epithelial cells, malignant-restricted: ", sum(cells_m$compartment == "epithelial")),
    paste0("panel genes evaluable in both: ", sum(p$evaluable_lineage & p$evaluable_malig)),
    paste0("panel genes dominant, lineage: ", sum(p$dominant_lineage, na.rm = TRUE)),
    paste0("panel genes dominant, malignant-restricted: ", sum(p$dominant_malignant, na.rm = TRUE)),
    paste0("panel genes whose dominance call CHANGED: ", sum(p$call_changed, na.rm = TRUE)),
    paste0("median delta f(0.30), panel: ",
           format(median(p$delta_f30, na.rm = TRUE), digits = 4))
  ), file.path(OUTDIR, "gse125449_malignant_sensitivity.readme.txt"))
  message(sprintf("  ..  A.4 sensitivity (GSE125449 only): %d HPC-like reassigned, %d/%d panel dominance calls changed",
                  sum(hp), sum(p$call_changed, na.rm = TRUE), sum(p$evaluable_lineage & p$evaluable_malig)))
})

# --- F19: GSE125449 dual compartment grouping
# The pilot used a five-compartment scheme with TEC in fibroblast_stromal;
# Amendment 4's six-compartment scheme puts TEC in endothelial.  f(pi) is
# unchanged by the regrouping -- it moves cells between two NON-epithelial
# compartments, so neither the numerator nor the denominator moves, the same
# algebra as Amendment 11.  Reported so pilot and Part A numbers are comparable.
local({
  a <- "GSE125449"
  cl <- dat[[a]]$cells
  six  <- table(factor(cl$compartment, levels = COMPARTMENTS))
  five <- six
  five["fibroblast_stromal"] <- five["fibroblast_stromal"] + five["endothelial"]
  five["endothelial"] <- 0L
  write.csv(data.frame(
    compartment          = COMPARTMENTS,
    n_cells_six_amdt4    = as.integer(six),
    n_cells_five_pilot   = as.integer(five),
    note                 = ifelse(COMPARTMENTS == "endothelial",
                            "TEC: endothelial under Amendment 4; folded into fibroblast_stromal in the pilot's 5-compartment scheme",
                            "")),
    file.path(OUTDIR, "gse125449_compartment_counts_dual.csv"), row.names = FALSE)

  # G03: cell counts alone do not let a reader reconcile a PILOT STROMAL FRACTION
  # with a Part A one -- the incomparability is in the fraction, because TEC's
  # transcripts moved out of the stromal denominator.  The folded share is a pure
  # column sum of the six-compartment fractions and is REPORTING ONLY: f(pi) is
  # untouched, since regrouping two non-epithelial compartments changes neither
  # its numerator nor its denominator.
  fr <- stats[[a]]$frac
  write.csv(data.frame(
    gene                       = rownames(fr),
    qualifying                 = rownames(fr) %in% PANEL_GENES,
    frac_epithelial            = fr[, "epithelial"],
    frac_stromal_six_amdt4     = fr[, "fibroblast_stromal"],
    frac_endothelial_six_amdt4 = fr[, "endothelial"],
    frac_stromal_five_pilot    = fr[, "fibroblast_stromal"] + fr[, "endothelial"],
    stringsAsFactors = FALSE),
    file.path(OUTDIR, "gse125449_stromal_fraction_dual.csv"), row.names = FALSE)
  message("  ..  wrote gse125449_stromal_fraction_dual.csv (pilot-comparable stromal share)")
})

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
