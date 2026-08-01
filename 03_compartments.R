# 03_compartments.R -----------------------------------------------------------
# Amendment 9 / analysis-plan v1.4.  Local files only.  This script is
# intentionally fail-closed: an absent annotation, unexpected count, or
# unmapped label stops the run before any compartment fraction or k is written.
# Part B is not called here.
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)
if (!requireNamespace("Matrix", quietly = TRUE)) stop("HALT: Matrix is required")
if (!requireNamespace("rhdf5", quietly = TRUE)) stop("HALT: rhdf5 is required")

halt <- function(...) stop(paste0("HALT: ", paste0(..., collapse = "")), call. = FALSE)
root <- normalizePath(".", mustWork = TRUE)
raw <- file.path(root, "data", "raw")
geo <- file.path(root, "data", "geo")

# A.a: the registered four-atlas set.  GSE155698 is deliberately not opened.
atlases <- data.frame(
  atlas = c("GSE125449", "GSE178341", "GSE183904", "Peng"),
  tissue = c("gastric", "colorectal", "liver_biliary", "pancreatic"),
  stringsAsFactors = FALSE
)
stopifnot(nrow(atlases) == 4, !"GSE155698" %in% atlases$atlas)

panel_file <- file.path(root, "data", "panel", "panel_locked.csv")
if (!file.exists(panel_file)) halt("missing locked panel: ", panel_file)
panel <- read.csv(panel_file, check.names = FALSE)
if (nrow(panel) != 152L || !"gene" %in% names(panel))
  halt("locked panel is not the expected 152-gene table")
genes <- unique(trimws(panel$gene))
if (length(genes) != 152L || any(!nzchar(genes))) halt("invalid locked gene list")

# A.g invariant: exhaustively check all 3^4 combinations of TRUE/FALSE/NA
# evidence before any observed dominance matrix is used.
evidence <- expand.grid(rep(list(c(FALSE, TRUE, NA)), 4), KEEP.OUT.ATTRS = FALSE)
for (i in seq_len(nrow(evidence))) {
  d <- as.logical(evidence[i, ])
  n_dom <- sum(d, na.rm = TRUE)
  n_eval <- sum(!is.na(d))
  vals <- c(n_dom == 4L && n_eval == 4L,
            n_dom >= 2L && n_eval == 4L,
            n_dom >= 2L && n_eval >= 3L,
            n_dom >= 2L)
  if (!all(vals[1] <= vals[2], vals[2] <= vals[3], vals[3] <= vals[4]))
    halt("k ordering failed in exhaustive 3^4 verification")
}
stopifnot(nrow(evidence) == 81L)

# A.c: fixed vocabulary mapping, with explicit judgement calls recorded.
map_compartment <- function(x) {
  y <- trimws(as.character(x))
  z <- rep(NA_character_, length(y))
  z[grepl("epithelial|malignant|tumou?r cell|ductal|enterocyte|goblet|paneth|gastric|pit|chief|parietal|acinar|endocrine|islet|hepatocyte|cholangiocyte|HPC-like", y, ignore.case = TRUE)] <- "epithelial"
  z[grepl("fibroblast|CAF|myofibroblast|stellate|pericyte|smooth muscle|mesothelial", y, ignore.case = TRUE)] <- "fibroblast_stromal"
  z[grepl("macrophage|TAM|monocyte|neutrophil|granulocyte|mast|dendritic|myeloid", y, ignore.case = TRUE)] <- "myeloid"
  z[grepl("T cell|CD4|CD8|Treg|NK|B cell|plasma|lymphocyte", y, ignore.case = TRUE)] <- "lymphoid"
  z[grepl("endothelial|TEC|lymphatic|vascular", y, ignore.case = TRUE)] <- "endothelial"
  z[grepl("unclassified|doublet|ambiguous|heterogeneous|enteric glial|Schwann", y, ignore.case = TRUE)] <- "other"
  z
}

# GSE125449: committed labels and expected sample structure.
sfiles <- list.files(file.path(geo, "GSE125449"), pattern = "samples\\.txt\\.gz$", full.names = TRUE)
if (length(sfiles) != 2L) halt("GSE125449 sample tables missing")
g125_meta <- do.call(rbind, lapply(sfiles, function(f) read.delim(gzfile(f), check.names = FALSE)))
if (!all(c("Sample", "Cell Barcode", "Type") %in% names(g125_meta))) halt("GSE125449 schema changed")
if (!all(map_compartment(unique(g125_meta$Type)) %in% c("epithelial", "fibroblast_stromal", "myeloid", "lymphoid", "endothelial", "other"))) halt("GSE125449 unmapped label")

# GSE178341: assert the registered specimen/patient propagation before loading X.
f178_meta <- file.path(raw, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz")
f178_cl <- file.path(raw, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz")
if (!all(file.exists(c(f178_meta, f178_cl)))) halt("GSE178341 metadata files missing")
m178 <- read.csv(gzfile(f178_meta), check.names = FALSE)
c178 <- read.csv(gzfile(f178_cl), check.names = FALSE)
if (nrow(m178) != 370115L || nrow(c178) != 370115L) halt("GSE178341 cell count mismatch")
if (!all(c("cellID", "SPECIMEN_TYPE", "PatientTypeID") %in% names(m178))) halt("GSE178341 patient fields missing")
if (!all(c("sampleID", "clTopLevel", "clMidwayPr", "cl295v11SubFull") %in% names(c178))) halt("GSE178341 cluster fields missing")
keep178 <- trimws(m178$SPECIMEN_TYPE) == "T"
if (sum(keep178) != 129L || length(unique(m178$PatientTypeID[keep178])) != 62L)
  halt("GSE178341 expected 129 tumour GSM channels and 62 patients not met")
if (anyDuplicated(m178$cellID) || anyDuplicated(c178$sampleID)) halt("GSE178341 identifiers are not unique")

# Peng: verify the declared Besca mirror and the celltype1 field.  Category
# decoding is deferred until all four annotation tables are available.
fp <- file.path(raw, "StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad")
if (!file.exists(fp)) halt("Peng Besca mirror missing")
obs_names <- rhdf5::h5ls(fp, recursive = TRUE)$name
if (!all(c("obs/celltype1", "obs/Patient", "obs/Type", "obs/CONDITION") %in% obs_names))
  halt("Peng required obs fields missing")
if (length(rhdf5::h5read(fp, "obs/celltype1")) != 57423L ||
    length(unique(rhdf5::h5read(fp, "obs/Patient"))) != 35L)
  halt("Peng expected 57,423 cells and 35 patients not met")

# GSE183904: fail closed.  The committed archive has expression matrices only;
# no annotation column/table is present from which the A.c map can be built.
f183 <- file.path(raw, "GSE183904_RAW.tar")
members183 <- system2("tar", c("-tf", shQuote(f183)), stdout = TRUE)
if (length(members183) != 40L || any(!grepl("\\.csv\\.gz$", members183)))
  halt("GSE183904 archive member structure changed")
halt("GSE183904 has 40 expression matrices but no committed cell-level annotation table; A.c cannot be completed without an amended/provenanced annotation source")

# The code below is intentionally unreachable until the GSE183904 halt is
# resolved by an authorised amendment.  It documents the locked A.d--A.g
# estimand without silently substituting data or labels.
pseudobulk_raw <- function(X, cells, genes) {
  comps <- c("epithelial", "fibroblast_stromal", "myeloid", "lymphoid", "endothelial", "other")
  S <- matrix(0, nrow = length(genes), ncol = length(comps), dimnames = list(genes, comps))
  for (cc in comps) { ix <- which(cells$compartment == cc); if (length(ix)) S[, cc] <- Matrix::rowSums(X[genes, ix, drop = FALSE]) }
  S
}
purity_grid <- seq(0.30, 0.70, by = 0.01)
stopifnot(length(purity_grid) == 41L)
