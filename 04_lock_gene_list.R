
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
## R2-NAINDEX: which(), not a bare logical.  `NA > 0` is NA, and R's `[` on an NA
## row index returns a row of all NAs rather than dropping it -- real rows would
## be lost (lowering a gene's max, possibly manufacturing an exclusion) while
## phantom NA rows vanished silently in tapply.
if (any(is.na(det$n_cells))) halt("detection table has NA n_cells")
d <- det[which(det$gene %in% g152 & det$n_cells > 0), ]  # annotated compartments only
mx <- tapply(d$pct_cells_detected, d$gene, function(x)
        if (all(is.na(x))) NA_real_ else max(x, na.rm=TRUE))
## R2-COVERAGE: tapply only creates groups for observed levels, so a panel gene
## with NO row at all is ABSENT from mx rather than NA -- it would be retained
## without ever being tested by rule 2.  setequal catches that; the range guard
## pins the units (proportion, not percent).
if (!setequal(names(mx), g152))
  halt("detection table does not cover the panel: missing ",
       paste(setdiff(g152, names(mx)), collapse=", "))
if (any(is.na(mx))) halt("panel gene(s) with no detection value anywhere: ",
                         paste(names(mx)[is.na(mx)], collapse=", "))
if (!all(mx >= 0 & mx <= 1)) halt("detection values outside [0,1]; units are not proportions")
## R2-A12: Amendment 12 drops six panel genes from GSE125449, so their max is
## taken over TWO atlases while every other gene gets three.  Two of them (CRLF2,
## LEP) are rule-2 exclusions, and a gene held to the threshold over fewer
## atlases is likelier to fail it.  Recorded per gene, not left implicit.
n_atl <- tapply(d$atlas[!is.na(d$pct_cells_detected)], d$gene[!is.na(d$pct_cells_detected)],
                function(a) length(unique(a)))
ex2 <- sort(names(mx)[mx < 0.01])
best <- sapply(g152, function(g){ s <- d[d$gene==g & !is.na(d$pct_cells_detected),]
  if (!nrow(s)) return(NA_character_)
  b <- s[which.max(s$pct_cells_detected),]; sprintf("%s/%s", b$atlas, b$compartment) })

## ---- rules 3.1 and 3.3: read GENE ANNOTATION ONLY from each cohort ----------
## COHORT-ID: assert cohort IDENTITY, not arity.  Rule 3.1 is defined against
## "every TCGA cohort in the study" -- the seven named in Amendments 8 and 10.  A
## stray or substituted *_se.rds would keep a count check true while evaluating
## the rule against the wrong annotation set.
COHORTS <- c("TCGA-CHOL","TCGA-COAD","TCGA-ESCA","TCGA-LIHC","TCGA-PAAD","TCGA-READ","TCGA-STAD")
coh <- sub("_se\\.rds$","",basename(Sys.glob("data/tcga/*_se.rds")))
if (!setequal(coh, COHORTS))
  halt("data/tcga does not hold exactly the seven registered cohorts. found: ",
       paste(sort(coh), collapse=", "), " | registered: ", paste(COHORTS, collapse=", "))
coh <- COHORTS                                   # fixed order, not glob order

pres  <- matrix(NA, length(g152), length(coh), dimnames=list(g152, coh))
chrom <- NULL; nrow_annot <- NULL
for (cc in coh) {
  se  <- readRDS(sprintf("data/tcga/%s_se.rds", cc))
  sym <- as.character(rowData(se)$gene_name)          # NOT colData, NOT assays
  pres[, cc] <- g152 %in% sym
  ch  <- setNames(as.character(GenomicRanges::seqnames(rowRanges(se))), sym)
  ## SYMBOL-MULTI: a symbol matching several annotation rows is collapsed to a
  ## slash-joined string.  That is right for the PAR case (chrX + the
  ## _PAR_Y row) but would silently pass for a symbol on two autosomes, and Part B
  ## would have no defined rule for which row to quantify.  sort() so the string
  ## is order-stable across cohorts; record multiplicity and assert it below.
  cur <- sapply(g152, function(g){v<-sort(unique(ch[names(ch)==g])); if(!length(v)) NA_character_ else paste(v,collapse="/")})
  nr  <- sapply(g152, function(g) sum(names(ch)==g))
  if (is.null(chrom)) { chrom <- cur; nrow_annot <- nr }
  else {
    ## R31-DEAD: compare only over genes PRESENT in both cohorts.  Previously an
    ## absent gene made `cur` differ from `chrom` and the script halted with
    ## "chromosome assignment disagrees" before ex31 was ever computed -- so rule
    ## 3.1 could only halt, never exclude, and its zero count was not an
    ## observation.  Verified: a gene missing from one cohort halted on the
    ## chromosome check in both orderings.
    both <- !is.na(chrom) & !is.na(cur)
    if (!identical(chrom[both], cur[both]))
      halt("chromosome assignment disagrees between cohorts (", cc, "): ",
           paste(names(chrom)[both][chrom[both] != cur[both]], collapse=", "),
           ". Resolve before locking.")
    chrom[is.na(chrom)] <- cur[is.na(chrom)]     # fill from any cohort that has it
    nrow_annot <- pmax(nrow_annot, nr)
  }
  rm(se); gc(verbose=FALSE)
}
ex31 <- sort(g152[!apply(pres, 1, all)])         # now reachable

## Chromosome is needed only for genes that survive rule 3.1.
if (any(is.na(chrom[setdiff(g152, ex31)])))
  halt("panel gene(s) present in every cohort but without a chromosome call: ",
       paste(setdiff(names(chrom)[is.na(chrom)], ex31), collapse=", "))
## Rule 3.3 reads "Sex-chromosome genes, since cohorts differ in sex composition."
## Four panel genes sit in the PSEUDOAUTOSOMAL REGION and are annotated on BOTH
## chrX and chrY (the annotation emits a second `<ENSG>_PAR_Y` row for each):
## CRLF2, CSF2RA, IL3RA, IL9R.  The deposit records `metadata(se)$data_release =
## "Data Release 45.0 - December 04, 2025"` and carries no GENCODE version field;
## an earlier revision of this comment cited "GENCODE v36", which was inferred
## from the _PAR_Y convention rather than read from the object, and is withdrawn.
## Whether rule 3.3 catches them is NOT settled by the registered text, and the
## two readings give different scoring sets (143 vs 140 genes):
##   NARROW - exclude genes annotated ONLY on X or Y (PAR genes retained). 143.
##   BROAD  - exclude anything with any row on a sex chromosome, on the rule's
##            stated TEXT. 140.
##
## AN EARLIER REVISION OF THIS COMMENT JUSTIFIED NARROW ON THE BIOLOGY, ARGUING
## THAT PAR GENES "escape X-inactivation and are present in two copies in both
## sexes, so cohort sex composition does not bias them". THAT ARGUMENT IS WRONG
## AT THE OPERATIVE STEP and has been withdrawn. Escaping X-inactivation is
## precisely the mechanism that PRODUCES sex-differential expression: an escapee
## is transcribed from both X copies in XX individuals, so it tends to be
## expressed HIGHER in females. Rule 3.3's stated reason therefore argues FOR
## excluding PAR genes, not for retaining them. The claim was asserted from
## memory and not checked; it should not have been written into a lock script.
##
## The question is consequently OPEN. This script implements NARROW only because
## it is what Amendment 14's registered "143-gene final list" describes -- and
## Amendment 14 was written before the PAR question was identified, so it records
## the number without deciding the question. That is not an argument, and the
## `PAR_UNDECIDED` flag below marks it as provisional.
##
## DIRECTION OF BIAS: narrow retains CSF2RA, IL3RA and IL9R in the scoring set.
## None is epithelial-dominant, so each adds a non-epithelial gene to the score
## -- the direction that favours this study's "substantially stromal" thesis.
## Narrow is therefore the hypothesis-friendly reading, and it is the one
## implemented; that is a reason for a human to decide it explicitly, not for the
## script to keep choosing silently. k is 43 under both readings, so no Part A
## quantity depends on it; the Amendment 3 branch is BUILT under 152, 143 and 140
## alike (30.3%, 30.1%, 30.7% of panel, all above the 20% floor).
PAR_UNDECIDED <- TRUE   # cleared only by an amendment deciding rule 3.3's scope
PAR_GENES  <- c("CRLF2", "CSF2RA", "IL3RA", "IL9R")
sex_any    <- sort(g152[grepl("chr[XY]", chrom)])
sex_only   <- sort(g152[grepl("^(chr)?[XY]$", chrom)])
par_seen   <- sort(setdiff(sex_any, sex_only))
if (!identical(par_seen, sort(PAR_GENES)))
  halt("the set of pseudoautosomal panel genes changed: observed ",
       paste(par_seen, collapse=", "), " vs the four this lock was built against (",
       paste(sort(PAR_GENES), collapse=", "), "). NOTE this halt has TWO possible ",
       "causes: the annotation changed, OR rule 3.3's scope was never decided on ",
       "the merits (see PAR_UNDECIDED). Re-decide the interpretation in a dated ",
       "amendment before relocking.")
## PAR-5: pin the chromosome string too, so a PAR gene reannotated to a single
## chromosome cannot pass the set check while changing what it means.
if (!all(chrom[PAR_GENES] %in% c("chrX/chrY", "chrY/chrX")))
  halt("a pseudoautosomal gene is no longer annotated on both sex chromosomes: ",
       paste(sprintf("%s=%s", PAR_GENES, chrom[PAR_GENES]), collapse=", "))
## SYMBOL-MULTI: every panel gene must map to exactly one locus, except the four
## PAR genes which carry the _PAR_Y second row.  Part B's expression
## lookup needs this contract to be explicit.
expect_rows <- ifelse(g152 %in% PAR_GENES, 2L, 1L)
if (!identical(unname(nrow_annot), expect_rows))
  halt("panel gene(s) with unexpected annotation-row multiplicity: ",
       paste(sprintf("%s=%d", g152[nrow_annot != expect_rows], nrow_annot[nrow_annot != expect_rows]),
             collapse=", "), ". A symbol matching several loci has no defined ",
       "quantification rule in Part B.")
ex33 <- sex_only                      # NARROW -- provisional, see PAR_UNDECIDED

## ---- assemble, in the registered order 1 -> 2 -> 3 --------------------------
excl <- unique(c(ex31, ex2, ex33))
final <- setdiff(g152, excl)
reason <- function(g) paste(c(if(g %in% ex31)"3.1_absent_from_a_TCGA_cohort_annotation",
                              if(g %in% ex2) "2_undetectable_under_1pct_in_every_compartment",
                              if(g %in% ex33)"3.3_sex_chromosome"), collapse="+")

## ---- assertions against the Part A / Amendment 14 record --------------------
stopifnot(length(ex31)==0L, length(ex2)==6L, length(ex33)==3L, length(final)==143L)
## ORIGIN6: assert which origin genes are IN the panel rather than hardcoding the
## conclusion.  Amendment 2 removed BCL2, MMP9 and HGF, so only three can appear
## in `excl` -- but if the panel were rebuilt and BCL2 requalified, a guard naming
## only three would not notice its exclusion.
ORIGIN6 <- c("SOCS3","BCL2","MYC","MMP9","HGF","IL6")
if (!setequal(intersect(ORIGIN6, g152), c("SOCS3","MYC","IL6")))
  halt("origin-six panel membership changed: ", paste(sort(intersect(ORIGIN6, g152)), collapse=", "),
       " (expected SOCS3, MYC, IL6; BCL2/MMP9/HGF are non-qualifying per Amendment 2)")
if (length(intersect(excl, ORIGIN6)))
  halt("an origin-six gene was excluded: ", paste(sort(intersect(excl, ORIGIN6)), collapse=", "))

dm <- read.csv("output/compartment_dominance_matrix.csv", stringsAsFactors=FALSE)
p  <- dm[dm$in_panel, ]
## K-DENOM: the k assertion re-derives from a separate committed Part A input, so
## it is not circular -- but it says nothing about the DENOMINATOR.  Pin it: the
## dominance matrix must cover exactly the panel, over exactly three atlases, with
## no NA dominance call standing in for a missing one.
n_atlas <- length(unique(p$atlas))
if (n_atlas != 3L) halt("dominance matrix covers ", n_atlas, " atlases, expected 3")
nd <- tapply(p$dominant, p$gene, function(x) sum(x, na.rm=TRUE))
if (!setequal(names(nd), g152))
  halt("dominance matrix in_panel rows do not cover the locked panel")
k152 <- sum(nd >= 2); k143 <- sum(nd[names(nd) %in% final] >= 2)
a152 <- sum(nd == n_atlas);  a143 <- sum(nd[names(nd) %in% final] == n_atlas)
if (k152 != 46L || k143 != 43L || a152 != 24L || a143 != 23L)
  halt("Amendment 14 records k=46/43 and k_all3=24/23; observed ",
       k152,"/",k143," and ",a152,"/",a143)

## ---- write --------------------------------------------------------------
## PAR-4: the retained PAR genes are the trace of an undecided interpretation, so
## they are FLAGGED in the list itself rather than only in prose -- previously the
## only evidence in either file was CRLF2's chromosome string, which a reader
## would have had to reverse-engineer.
fl <- data.frame(gene=final, in_locked_panel_152=TRUE,
                 origin_six=final %in% ORIGIN6,
                 chromosome=unname(chrom[final]),
                 pseudoautosomal=final %in% PAR_GENES,
                 max_detection_pct=round(100*unname(mx[final]),4),
                 n_atlases_evaluable=unname(n_atl[final]),
                 n_annotation_rows=unname(nrow_annot[final]),
                 stringsAsFactors=FALSE)
write.csv(fl, file.path(OUT,"final_gene_list.csv"), row.names=FALSE)
ex <- data.frame(gene=excl, excluded_by_rule=sapply(excl, reason),
                 chromosome=unname(chrom[excl]),
                 pseudoautosomal=excl %in% PAR_GENES,
                 max_detection_pct=round(100*unname(mx[excl]),4),
                 max_at=unname(best[excl]),                  # R2-A12: where the max was
                 n_atlases_evaluable=unname(n_atl[excl]),    # R2-A12: Amendment 12 exposure
                 route_into_panel=panel$route[match(excl, panel$gene)],
                 dominant_in_n_atlases=unname(nd[excl]), stringsAsFactors=FALSE)
ex <- ex[order(ex$excluded_by_rule, ex$gene), ]
write.csv(ex, file.path(OUT,"final_gene_list_exclusions.csv"), row.names=FALSE)
cat("locked:", length(final), "genes |", nrow(ex), "excluded\n")
print(ex, row.names=FALSE)
