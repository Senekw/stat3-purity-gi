
## Lock the final gene list. Derives everything from committed inputs; asserts
## every count against the Amendment 14 / Part A record rather than trusting it.
suppressPackageStartupMessages(library(SummarizedExperiment))
halt <- function(...) stop(paste0("HALT: ", paste0(c(...), collapse="")), call.=FALSE)
OUT <- "data/panel"

panel <- read.csv("data/panel/panel_locked.csv", stringsAsFactors=FALSE)
if (nrow(panel) != 152L) halt("locked panel is not 152 rows: ", nrow(panel))
g152 <- panel$gene

## ---- rule 2: <1% of cells in EVERY annotated compartment, EVERY atlas -------
det <- read.csv("output/detection_rate_by_compartment.csv", stringsAsFactors=FALSE)
d <- det[det$gene %in% g152 & det$n_cells > 0, ]        # annotated compartments only
mx <- tapply(d$pct_cells_detected, d$gene, function(x)
        if (all(is.na(x))) NA_real_ else max(x, na.rm=TRUE))
if (any(is.na(mx))) halt("panel gene(s) with no detection value anywhere: ",
                         paste(names(mx)[is.na(mx)], collapse=", "))
ex2 <- sort(names(mx)[mx < 0.01])
best <- sapply(ex2, function(g){ s <- d[d$gene==g & !is.na(d$pct_cells_detected),]
  b <- s[which.max(s$pct_cells_detected),]; sprintf("%.4f%% (%s/%s)", 100*b$pct_cells_detected, b$atlas, b$compartment) })

## ---- rules 3.1 and 3.3: read GENE ANNOTATION ONLY from each cohort ----------
coh <- sub("_se\\.rds$","",basename(Sys.glob("data/tcga/*_se.rds")))
if (length(coh) != 7L) halt("expected 7 TCGA cohorts, found ", length(coh))
pres <- matrix(NA, length(g152), length(coh), dimnames=list(g152, coh)); chrom <- NULL
for (cc in coh) {
  se  <- readRDS(sprintf("data/tcga/%s_se.rds", cc))
  sym <- as.character(rowData(se)$gene_name)          # NOT colData, NOT assays
  pres[, cc] <- g152 %in% sym
  ch  <- setNames(as.character(GenomicRanges::seqnames(rowRanges(se))), sym)
  cur <- sapply(g152, function(g){v<-unique(ch[names(ch)==g]); if(!length(v)) NA_character_ else paste(v,collapse="/")})
  if (is.null(chrom)) chrom <- cur
  else if (!identical(chrom, cur))
    halt("chromosome assignment disagrees between cohorts (", cc, "); resolve before locking")
  rm(se); gc(verbose=FALSE)
}
ex31 <- sort(g152[!apply(pres, 1, all)])
if (any(is.na(chrom))) halt("panel gene(s) without a chromosome call: ",
                            paste(names(chrom)[is.na(chrom)], collapse=", "))
## Rule 3.3 reads "Sex-chromosome genes, since cohorts differ in sex composition."
## Four panel genes sit in the PSEUDOAUTOSOMAL REGION and are annotated on BOTH
## chrX and chrY (GENCODE emits a second _PAR_Y row): CRLF2, CSF2RA, IL3RA, IL9R.
## Whether rule 3.3 catches them is NOT settled by the registered text, and the
## two readings give different scoring sets (143 vs 140 genes):
##   NARROW - exclude genes annotated ONLY on X or Y. PAR genes are retained,
##            on the rule's stated REASON: they escape X-inactivation and are
##            present in two copies in both sexes, so cohort sex composition
##            does not bias them.
##   BROAD  - exclude anything with any row on a sex chromosome, on the rule's
##            stated TEXT.
## Amendment 14 registers the final list at 143, which is the NARROW reading, so
## that is what is implemented here. The choice is recorded in
## final_gene_list_exclusions.csv and NOTES_FOR_REVIEW section 17; it was
## surfaced AFTER Amendment 14 was written and is not settled by it.
## k is 43 under both readings, so no reported Part A quantity depends on it.
## The three genes it moves are CSF2RA, IL3RA and IL9R, none of them dominant.
PAR_GENES  <- c("CRLF2", "CSF2RA", "IL3RA", "IL9R")
sex_any    <- sort(g152[grepl("chr[XY]", chrom)])
sex_only   <- sort(g152[grepl("^(chr)?[XY]$", chrom)])
par_seen   <- sort(setdiff(sex_any, sex_only))
if (!identical(par_seen, sort(PAR_GENES)))
  halt("the set of pseudoautosomal panel genes changed: observed ",
       paste(par_seen, collapse=", "), " vs registered ", paste(sort(PAR_GENES), collapse=", "),
       ". The rule 3.3 interpretation was settled against the old set; re-decide it.")
ex33 <- sex_only                      # NARROW reading, per Amendment 14

## ---- assemble, in the registered order 1 -> 2 -> 3 --------------------------
excl <- unique(c(ex31, ex2, ex33))
final <- setdiff(g152, excl)
reason <- function(g) paste(c(if(g %in% ex31)"3.1_absent_from_a_TCGA_cohort_annotation",
                              if(g %in% ex2) "2_undetectable_under_1pct_in_every_compartment",
                              if(g %in% ex33)"3.3_sex_chromosome"), collapse="+")

## ---- assertions against the Part A / Amendment 14 record --------------------
stopifnot(length(ex31)==0L, length(ex2)==6L, length(ex33)==3L, length(final)==143L)
if (length(intersect(excl, c("SOCS3","MYC","IL6")))) halt("an origin-six panel gene was excluded")
dm <- read.csv("output/compartment_dominance_matrix.csv", stringsAsFactors=FALSE)
p  <- dm[dm$in_panel, ]
nd <- tapply(p$dominant, p$gene, function(x) sum(x, na.rm=TRUE))
k152 <- sum(nd >= 2); k143 <- sum(nd[names(nd) %in% final] >= 2)
a152 <- sum(nd == 3);  a143 <- sum(nd[names(nd) %in% final] == 3)
if (k152 != 46L || k143 != 43L || a152 != 24L || a143 != 23L)
  halt("Amendment 14 records k=46/43 and k_all3=24/23; observed ",
       k152,"/",k143," and ",a152,"/",a143)

## ---- write --------------------------------------------------------------
fl <- data.frame(gene=final, in_locked_panel_152=TRUE,
                 origin_six=final %in% c("SOCS3","BCL2","MYC","MMP9","HGF","IL6"),
                 chromosome=unname(chrom[final]),
                 max_detection_pct=round(100*unname(mx[final]),4),
                 stringsAsFactors=FALSE)
write.csv(fl, file.path(OUT,"final_gene_list.csv"), row.names=FALSE)
ex <- data.frame(gene=excl, excluded_by_rule=sapply(excl, reason),
                 chromosome=unname(chrom[excl]),
                 max_detection_pct=round(100*unname(mx[excl]),4),
                 route_into_panel=panel$route[match(excl, panel$gene)],
                 dominant_in_n_atlases=unname(nd[excl]), stringsAsFactors=FALSE)
ex <- ex[order(ex$excluded_by_rule, ex$gene), ]
write.csv(ex, file.path(OUT,"final_gene_list_exclusions.csv"), row.names=FALSE)
cat("locked:", length(final), "genes |", nrow(ex), "excluded\n")
print(ex, row.names=FALSE)
