# 01_download.R ---------------------------------------------------------------
# Data retrieval for: compartment decomposition of STAT3 activity scores in GI cancers
#
# Run once. Everything is cached to disk as .rds so later scripts never re-download.
# Expect 1-3 hours and ~15-25 GB for the TCGA step depending on connection.
#
# NOT TESTED - written without an R environment available. Run section by section
# the first time and fix as you go.
# -----------------------------------------------------------------------------

# ---- setup ------------------------------------------------------------------

dirs <- c("data/raw", "data/tcga", "data/geo", "data/manual", "data/cache")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

needed_cran <- c("readxl", "dplyr", "readr", "Matrix", "R.utils")
needed_bioc <- c("TCGAbiolinks", "SummarizedExperiment", "GEOquery")

for (p in needed_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in needed_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)   # tribble() is used in section 3; dplyr re-exports it but be explicit

options(timeout = 3600)  # GEO/GDC downloads exceed the 60s default and fail silently otherwise

# FIX (2026-07-31): R's default internal url() method fails behind this machine's
# HTTP proxy - api.gdc.cancer.gov returns "HTTP response code said error" and
# TCGAbiolinks misreports it as "GDC server down". libcurl honours the proxy env
# vars correctly. Must be set BEFORE any GDCquery call.
options(download.file.method = "libcurl", url.method = "libcurl")


# ---- 1. TCGA expression -----------------------------------------------------
# STAR-Counts gives raw counts, TPM, and FPKM in one object. Keep the whole
# SummarizedExperiment: colData carries the barcodes you need for the CDR merge.

projects <- c("TCGA-COAD", "TCGA-READ", "TCGA-STAD", "TCGA-ESCA",
              "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL")

# FIX (2026-07-31): the GDC API through this machine's proxy fails intermittently
# (roughly 1 in 5 calls returns a connection error / 15s timeout). Unretried, a
# multi-hour download aborts on a transient blip. Retry with backoff.
with_retry <- function(expr, tries = 5, wait = 5, label = "") {
  for (i in seq_len(tries)) {
    r <- try(force(expr), silent = TRUE)
    if (!inherits(r, "try-error")) return(r)
    msg <- conditionMessage(attr(r, "condition"))
    message(sprintf("  [retry %d/%d] %s: %s", i, tries, label, substr(msg, 1, 90)))
    if (i < tries) Sys.sleep(wait * i)
  }
  stop("all ", tries, " attempts failed for ", label)
}

fetch_tcga <- function(project) {
  out <- file.path("data/tcga", paste0(project, "_se.rds"))
  if (file.exists(out)) {
    message("cached: ", project); return(invisible(out))
  }
  message("downloading: ", project)
  q <- with_retry(GDCquery(
    project      = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  ), label = paste("GDCquery", project))
  with_retry(GDCdownload(q, directory = "data/raw"),
             label = paste("GDCdownload", project))
  # GDCprepare's final step enriches colData from api.gdc.cancer.gov/cases and is
  # the most proxy-fragile call in the pipeline (TCGAbiolinks' own internal retry
  # does not recover from it). Raw files are already cached at this point, so
  # re-attempts are cheap: only the clinical lookup is repeated. Retry harder.
  se <- with_retry(GDCprepare(q, directory = "data/raw"),
                   tries = 12, wait = 20,
                   label = paste("GDCprepare", project))
  saveRDS(se, out)
  rm(se); gc()
  invisible(out)
}

for (p in projects) {
  tryCatch(fetch_tcga(p),
           error = function(e) message("FAILED ", p, ": ", conditionMessage(e)))
}

# Sanity check - assay names should include unstranded, tpm_unstrand, fpkm_unstrand
# FIX (2026-07-31): guarded. Original read TCGA-CHOL_se.rds unconditionally, so if
# that one cohort failed the script died here instead of reporting which succeeded.
chk <- file.path("data/tcga", "TCGA-CHOL_se.rds")
if (file.exists(chk)) {
  se_chk <- readRDS(chk)
  print(assayNames(se_chk))
  print(dim(se_chk))
  rm(se_chk); gc()
} else {
  message("TCGA-CHOL not downloaded; skipping assay sanity check")
}
message("cohorts on disk: ",
        paste(sub("_se\\.rds$", "", list.files("data/tcga", pattern = "_se\\.rds$")),
              collapse = ", "))


# ---- 2. GEO validation cohorts ---------------------------------------------
# GSE39582  - CRC, Marisa et al., ~560 with survival
# GSE125449 - liver cancer single cell (iCCA + HCC), Ma et al., the compartment atlas
# ACRG gastric is GSE66229 (expression) / GSE62254 (the commonly cited accession);
#   check which one carries the survival annotation you need before relying on it.

library(GEOquery)

fetch_geo_series <- function(acc) {
  out <- file.path("data/geo", paste0(acc, "_eset.rds"))
  if (file.exists(out)) { message("cached: ", acc); return(invisible(out)) }
  message("downloading series matrix: ", acc)
  gse <- getGEO(acc, GSEMatrix = TRUE, destdir = "data/geo")
  saveRDS(gse, out)
  invisible(out)
}

for (acc in c("GSE39582", "GSE66229")) {
  tryCatch(fetch_geo_series(acc),
           error = function(e) message("FAILED ", acc, ": ", conditionMessage(e)))
}

# Single cell comes as supplementary files, not a series matrix.
# GSE125449 ships 10x-style triplets (matrix.mtx / barcodes.tsv / genes.tsv) plus a
# per-cell annotation table with the author labels you need:
#   Malignant cell, HPC-like, CAF, TEC, TAM, T cell, B cell
if (!length(list.files("data/geo/GSE125449", recursive = TRUE))) {
  getGEOSuppFiles("GSE125449", baseDir = "data/geo", makeDirectory = TRUE)
  # untar/gunzip whatever arrives; the layout has changed before, so inspect first
  print(list.files("data/geo/GSE125449", recursive = TRUE))
}


# ---- 3. Annotation tables that need a browser -------------------------------
# No stable API. Download by hand into data/manual/ with these exact filenames,
# then re-run this section to verify.

manual <- tribble(
  ~file,                  ~what,                              ~where,
  "TCGA-CDR.xlsx",        "Curated survival endpoints (Liu 2018)",
                          "GDC PanCanAtlas publications page - TCGA-CDR supplemental table S1",
  "aran_purity.xlsx",     "Consensus purity estimate (CPE)",
                          "Aran 2015 Nat Commun ncomms9971 - Supplementary Data 1",
  "dong_icca_tableS1.xlsx","FU-iCCA sample annotation + normalised data",
                          "Dong 2022 Cancer Cell PMID 34971568 - Table S1",
  "NODE_OEP001105/",      "FU-iCCA bulk RNA-seq / proteome / phosphoproteome",
                          "biosino.org/node/project/detail/OEP001105 (registration required)"
)

missing <- manual$file[!file.exists(file.path("data/manual", manual$file))]
if (length(missing)) {
  message("\n--- STILL NEEDED (download by hand) ---")
  print(as.data.frame(manual[manual$file %in% missing, ]), right = FALSE)
} else {
  message("all manual downloads present")
}

# Once TCGA-CDR.xlsx is in place:
if (file.exists("data/manual/TCGA-CDR.xlsx")) {
  cdr <- readxl::read_excel("data/manual/TCGA-CDR.xlsx", sheet = 1)
  cdr <- cdr %>% filter(type %in% gsub("TCGA-", "", projects))
  saveRDS(cdr, "data/cache/cdr.rds")
  message("CDR rows retained: ", nrow(cdr))
  # Liu et al. give per-cancer endpoint recommendations - honour them.
  # OS is underpowered in COAD/READ; use PFI there.
}


# ---- 4. provenance ----------------------------------------------------------
# Preregistration will need exact versions. Write them down now, not later.

writeLines(
  c(paste("retrieved:", Sys.Date()),
    paste("R:", R.version.string),
    paste("TCGAbiolinks:", as.character(packageVersion("TCGAbiolinks"))),
    paste("GEOquery:", as.character(packageVersion("GEOquery"))),
    "",
    "GDC data release: check GDCquery output above and record it here"),
  "data/cache/provenance.txt"
)

sessionInfo()
