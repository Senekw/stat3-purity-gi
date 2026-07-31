# 02_panel.R ------------------------------------------------------------------
# Implements panel_definition.md. Run once, commit the output, do not re-run
# with different parameters after the sweep.
#
# NOT TESTED - no R available where this was written. Expect to debug.
# -----------------------------------------------------------------------------

library(dplyr)

dir.create("data/panel", recursive = TRUE, showWarnings = FALSE)

# ---- criterion A: MSigDB hallmark -------------------------------------------

if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
library(msigdbr)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING") %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

# msigdbr's argument names changed across versions (category= vs collection=).
# If the above errors, check msigdbr_collections() and adjust.

message("criterion A genes: ", length(hallmark))
msigdb_version <- as.character(packageVersion("msigdbr"))


# ---- criterion B: direct targets -------------------------------------------
# Implements panel_definition.md §2B AS AMENDED 2026-07-31 (Amendment 1):
#   >= 1 human ChIP-seq derived source AND independent curation in TRRUST v2.
#
# ChEA3's public API is enrich-only (no per-TF target endpoint), so the underlying
# GMT libraries are pulled directly and the STAT3 terms extracted.
#
# Independence handling, per §2B ("two sources querying the same underlying
# ChIP-seq experiment count as one"):
#   - ARCHS4_Coexpression, GTEx_Coexpression: EXCLUDED. Coexpression is not
#     binding evidence.
#   - Enrichr_Queries: EXCLUDED. User-submitted gene lists, not primary evidence.
#   - Literature_ChIP-seq: STAT3_1855785_CHIPSEQ_MESC_MOUSE dropped as a
#     truncated-PMID duplicate of STAT3_18555785_CHIPSEQ_MESC_MOUSE.
#   - Literature_ChIP-seq split by species. Only the human experiment
#     (PMID 23295773, U87) counts; mouse terms are recorded but do not count.

chea3_base <- "https://maayanlab.cloud/chea3/assets/tflibs"
chea3_libs_chip <- c("ENCODE_ChIP-seq", "ReMap_ChIP-seq", "Literature_ChIP-seq")
chea3_libs_excluded <- c("ARCHS4_Coexpression" = "coexpression, not binding",
                         "GTEx_Coexpression"   = "coexpression, not binding",
                         "Enrichr_Queries"     = "user query lists, not primary")

read_gmt_terms <- function(path, tf = "STAT3") {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (!length(p)) next
    term <- p[1]
    if (toupper(sub("_.*$", "", term)) != tf) next
    genes <- sub(",.*$", "", trimws(p[-(1:2)]))
    out[[term]] <- genes[nzchar(genes)]
  }
  out
}

chip_terms <- list()
for (L in chea3_libs_chip) {
  f <- file.path("data/panel", paste0(L, ".gmt"))
  if (!file.exists(f)) {
    download.file(file.path(chea3_base, paste0(L, ".gmt")), f, quiet = TRUE)
  }
  chip_terms[[L]] <- read_gmt_terms(f)
  message(L, ": ", length(chip_terms[[L]]), " STAT3 terms, ",
          length(unique(unlist(chip_terms[[L]]))), " genes")
}

# duplicate-PMID collapse
chip_terms$`Literature_ChIP-seq`$STAT3_1855785_CHIPSEQ_MESC_MOUSE <- NULL
lit_names  <- names(chip_terms$`Literature_ChIP-seq`)
lit_human  <- lit_names[grepl("_HUMAN$", lit_names)]
lit_mouse  <- lit_names[grepl("_MOUSE$", lit_names)]

chip_groups <- list(
  ENCODE_ChIP_seq           = unique(unlist(chip_terms$`ENCODE_ChIP-seq`)),
  ReMap_ChIP_seq            = unique(unlist(chip_terms$`ReMap_ChIP-seq`)),
  Literature_ChIP_seq_human = unique(unlist(chip_terms$`Literature_ChIP-seq`[lit_human]))
)
chip_mouse <- unique(unlist(chip_terms$`Literature_ChIP-seq`[lit_mouse]))

# TRRUST v2 human
trrust_file <- "data/panel/trrust_rawdata.human.tsv"
if (!file.exists(trrust_file)) {
  download.file("https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
                trrust_file, quiet = TRUE)
}
trrust <- read.delim(trrust_file, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("TF", "target", "mode", "pmids"))
trrust_stat3 <- sort(unique(trrust$target[trrust$TF == "STAT3"]))
message("TRRUST v2: ", nrow(trrust), " edges; STAT3 targets: ", length(trrust_stat3))

chip_any_human <- unique(unlist(chip_groups))
criterion_B <- intersect(chip_any_human, trrust_stat3)   # Amendment 1: AND, not any-two
message("criterion B genes (ChIP-seq AND TRRUST): ", length(criterion_B))

chea3_query_date <- Sys.Date()


# ---- assemble ---------------------------------------------------------------

origin_six <- c("SOCS3", "BCL2", "MYC", "MMP9", "HGF", "IL6")

panel <- tibble(gene = union(hallmark, criterion_B)) %>%
  mutate(
    via_A        = gene %in% hallmark,
    via_B        = gene %in% criterion_B,
    route        = case_when(via_A & via_B ~ "both", via_A ~ "A_only", TRUE ~ "B_only"),
    origin_score = gene %in% origin_six,
    # per-source evidence, so every gene's inclusion route is auditable (§7)
    ev_ENCODE_ChIPseq      = gene %in% chip_groups$ENCODE_ChIP_seq,
    ev_ReMap_ChIPseq       = gene %in% chip_groups$ReMap_ChIP_seq,
    ev_Literature_human    = gene %in% chip_groups$Literature_ChIP_seq_human,
    ev_Literature_mouse    = gene %in% chip_mouse,      # recorded, does not count
    ev_TRRUST_v2           = gene %in% trrust_stat3,
    n_chipseq_human        = ev_ENCODE_ChIPseq + ev_ReMap_ChIPseq + ev_Literature_human
  ) %>%
  arrange(desc(origin_score), route, gene)

# The original six, with full evidence, whether or not they qualified.
origin_report <- tibble(gene = origin_six) %>%
  mutate(
    via_A               = gene %in% hallmark,
    ev_ENCODE_ChIPseq   = gene %in% chip_groups$ENCODE_ChIP_seq,
    ev_ReMap_ChIPseq    = gene %in% chip_groups$ReMap_ChIP_seq,
    ev_Literature_human = gene %in% chip_groups$Literature_ChIP_seq_human,
    ev_Literature_mouse = gene %in% chip_mouse,
    ev_TRRUST_v2        = gene %in% trrust_stat3,
    n_chipseq_human     = ev_ENCODE_ChIPseq + ev_ReMap_ChIPseq + ev_Literature_human,
    via_B               = gene %in% criterion_B,
    in_panel            = via_A | via_B
  )
readr::write_csv(origin_report, "data/panel/origin_six_evidence.csv")
print(as.data.frame(origin_report))

# Any of the original six that did NOT qualify - report this, do not add them back
dropped_origin <- setdiff(origin_six, panel$gene)
if (length(dropped_origin)) {
  message("ORIGINAL SCORE GENES FAILING PANEL CRITERIA: ",
          paste(dropped_origin, collapse = ", "))
}


# ---- exclusions (sections 3.1-3.3) ------------------------------------------
# Run after 01_download.R. Each exclusion is recorded with its reason so the
# final table explains every gene that left the panel.

apply_exclusions <- function(panel, tcga_gene_sets, sc_detection, sex_genes) {
  # tcga_gene_sets : named list, one character vector of symbols per cohort
  # sc_detection   : data.frame(gene, compartment, frac_cells_detected)
  # sex_genes      : character vector of X/Y-linked symbols

  in_all_tcga <- Reduce(intersect, tcga_gene_sets)

  detected <- sc_detection %>%
    group_by(gene) %>%
    summarise(max_frac = max(frac_cells_detected), .groups = "drop") %>%
    filter(max_frac >= 0.01) %>%
    pull(gene)

  panel %>%
    mutate(excluded_reason = case_when(
      !gene %in% in_all_tcga ~ "absent_from_a_tcga_cohort",
      !gene %in% detected    ~ "undetected_in_single_cell",
      gene %in% sex_genes    ~ "sex_chromosome",
      TRUE                   ~ NA_character_
    ))
}

# panel <- apply_exclusions(panel, tcga_gene_sets, sc_detection, sex_genes)
# final <- filter(panel, is.na(excluded_reason))


# ---- lock -------------------------------------------------------------------

readr::write_csv(panel, "data/panel/panel_locked.csv")

writeLines(
  c(paste("locked:", Sys.Date()),
    paste("R:", R.version.string),
    paste("msigdbr:", msigdb_version),
    "criterion B rule: panel_definition.md Amendment 1 (2026-07-31) -",
    "  >=1 human ChIP-seq source (ChEA3, default thresholds) AND TRRUST v2 curation",
    paste("criterion A n =", length(hallmark)),
    paste("criterion B n =", length(criterion_B)),
    paste("panel n =", nrow(panel)),
    paste("  route both   =", sum(panel$route == "both")),
    paste("  route A_only =", sum(panel$route == "A_only")),
    paste("  route B_only =", sum(panel$route == "B_only")),
    paste("original-six dropped:",
          if (length(dropped_origin)) paste(dropped_origin, collapse = ",") else "none"),
    "",
    paste("ChEA3 GMT libraries queried:", chea3_query_date),
    "  source: https://maayanlab.cloud/chea3/assets/tflibs/ (API is enrich-only;",
    "  ChEA3 does not version its GMTs - file sizes recorded as fingerprint)",
    paste("   ENCODE_ChIP-seq.gmt bytes    =", file.size("data/panel/ENCODE_ChIP-seq.gmt")),
    paste("   ReMap_ChIP-seq.gmt bytes     =", file.size("data/panel/ReMap_ChIP-seq.gmt")),
    paste("   Literature_ChIP-seq.gmt bytes=", file.size("data/panel/Literature_ChIP-seq.gmt")),
    paste("   ENCODE genes    =", length(chip_groups$ENCODE_ChIP_seq)),
    paste("   ReMap genes     =", length(chip_groups$ReMap_ChIP_seq)),
    paste("   Literature human genes =", length(chip_groups$Literature_ChIP_seq_human),
          "(terms:", paste(lit_human, collapse = ";"), ")"),
    paste("   Literature mouse genes =", length(chip_mouse), "(recorded, not counted)"),
    "  collapsed: STAT3_1855785_CHIPSEQ_MESC_MOUSE dropped as truncated-PMID",
    "    duplicate of STAT3_18555785_CHIPSEQ_MESC_MOUSE",
    "  excluded libraries (not direct-binding evidence):",
    paste("   ", names(chea3_libs_excluded), "-", chea3_libs_excluded, collapse = "\n"),
    "",
    paste("TRRUST v2 human queried:", chea3_query_date),
    "  source: https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    paste("   bytes =", file.size(trrust_file), "| total edges =", nrow(trrust),
          "| STAT3 targets =", length(trrust_stat3))),
  "data/panel/panel_provenance.txt"
)

print(panel, n = Inf)
