
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
## R2-A12: rule 2's evidence base is NOT equal across genes, and TWO mechanisms
## compound.  Amendment 12 drops six panel genes from GSE125449 outright, and the
## per-atlas evidence floor (>=20 summed counts) then removes further atlases from
## individual genes.  Measured on this run: 142 panel genes are evaluated against
## all three atlases, 7 against two, and 3 against ONE.
##   An earlier revision of this comment said the Amendment 12 genes "get two
##   atlases while every other gene gets three".  Both halves are wrong. CRLF2 and
##   LEP have ONE evaluable atlas, not two -- Amendment 12 removes GSE125449 and
##   the floor removes Peng -- and three genes with only two atlases (INHBE,
##   OPRM1, PAX3) are not Amendment 12 genes at all.  The single-atlas set is
##   exactly CRLF2, DNTT and LEP; DNTT reaches n=1 through the evidence floor
##   alone, since it is not an Amendment 12 gene.
## This matters because a gene held to the 1% threshold over fewer atlases has
## fewer chances to clear it. Recorded per gene in n_atlases_evaluable rather than
## left implicit; asserted below so the disparity cannot drift unnoticed.
n_atl <- tapply(d$atlas[!is.na(d$pct_cells_detected)], d$gene[!is.na(d$pct_cells_detected)],
                function(a) length(unique(a)))
if (!setequal(names(n_atl), g152))
  halt("atlas-evaluability count does not cover the panel")
## Pin the observed disparity: 142 genes on three atlases, 7 on two, 3 on one.
## A change here means the deposits or the evidence floor have moved, and rule 2
## would then be applied with different power than the locked list was built on.
atl_tab <- table(factor(n_atl[g152], levels = 1:3))
if (!identical(as.integer(atl_tab), c(3L, 7L, 142L)))
  halt("rule-2 atlas coverage changed: ", paste(sprintf("%s atlas=%d", names(atl_tab), atl_tab), collapse=", "),
       " (locked against 1 atlas=3, 2 atlases=7, 3 atlases=142). Re-check which genes ",
       "lost coverage before relocking; rule 2's power is not uniform across genes.")
if (!setequal(names(n_atl)[n_atl == 1L], c("CRLF2", "DNTT", "LEP")))
  halt("the single-atlas gene set changed: ", paste(sort(names(n_atl)[n_atl == 1L]), collapse=", "),
       " (locked against CRLF2, DNTT, LEP)")
ex2 <- sort(names(mx)[mx < 0.01])
best <- sapply(g152, function(g){ s <- d[d$gene==g & !is.na(d$pct_cells_detected),]
  if (!nrow(s)) return(NA_character_)
  b <- s[which.max(s$pct_cells_detected),]; sprintf("%s/%s", b$atlas, b$compartment) })

## ---- rules 3.1 and 3.3: read GENE ANNOTATION ONLY from each cohort ----------
## COHORT-ID: assert cohort IDENTITY, not arity.  Rule 3.1 is defined against
## "every TCGA cohort in the study".  Those seven are registered in
## analysis_plan.md -- its endpoint table lists all seven and its section on PFI
## says "in six of the seven cohorts (COAD, READ, ESCA, PAAD, LIHC, CHOL)".
## NOTE panel_definition.md names six of the seven (CHOL, COAD, ESCA, LIHC, READ,
## STAD -- Amendment 7 turns on COAD and LIHC explicitly) but never mentions PAAD,
## so the plan, not the panel document, is the authority for the full cohort set.
## A stray or substituted *_se.rds would keep a
## count check true while evaluating rule 3.1 against the wrong annotation set.
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
  if (is.null(chrom)) {
    chrom <- cur; nrow_annot <- nr
    ## PAR region assignment needs coordinates, not just chromosome names.  Taken
    ## from the first cohort; the cross-cohort identity check below covers the
    ## chromosome, and all seven derive from the same annotation release.
    gr0        <- rowRanges(se)
    gene_chr   <- as.character(GenomicRanges::seqnames(gr0))
    gene_start <- setNames(GenomicRanges::start(gr0), sym)
    gene_end   <- setNames(GenomicRanges::end(gr0),   sym)
  }
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
## AMENDMENT 15, 2026-08-02 -- DECIDED: the BROAD reading. Rule 3.3 is applied on
## its literal text: any panel gene whose annotation places it on chrX or chrY is
## excluded, PAR genes included. PRIMARY scoring set = 140. The 143-gene narrow
## list is retained as a PRESPECIFIED SENSITIVITY set, written alongside.
##
##   BROAD  (PRIMARY, Amendment 15) - any row on chrX or chrY. 7 genes qualify;
##            CRLF2 is already a rule-2 exclusion, so 12 excluded, 140 retained.
##   NARROW (SENSITIVITY)           - only genes annotated solely on X or Y. 143.
##
## Amendment 15's ground is that the registered text is categorical and names no
## exception, and that the biological argument the alternative needs does not hold
## uniformly: PAR1 genes have a functional Y homolog and are broadly
## dosage-balanced, while PAR2 genes are not uniformly XCI-escaping.
##   That PAR1/PAR2 split belongs to Amendment 15's NEGATIVE argument -- its reason
##   for declining the biological route -- not to the decision itself. The decision
##   rests on the registered text, which is categorical and excludes all four genes
##   whether they span one region or two. The split is recomputed below because a
##   claim in the amendment record should be checkable, NOT because the outcome
##   depends on it: were all four PAR1, the broad reading would be unchanged.
##
## TWO JUSTIFICATIONS WERE OFFERED AND BOTH WITHDRAWN; Amendment 15 relies on
## neither, and neither is used here.
##   (1) Mine, for retaining PAR genes: "they escape X-inactivation and are
##       present in two copies in both sexes." Asserted from memory, not checked,
##       and written into this script -- the error that made the audit necessary.
##   (2) The audit's replacement: "XCI escape is the mechanism that PRODUCES
##       female-biased expression." True of non-PAR X-linked escapees, but not of
##       PAR1, which has a functional Y homolog.
## Amendment 15 decides the question on the registered TEXT instead, which is
## categorical. That is the ground implemented below.
##
## DIRECTION OF BIAS (Amendment 15): CONSERVATIVE. The three additional exclusions
## -- CSF2RA, IL3RA, IL9R -- are all non-epithelial, so retaining them would make
## the score a MORE stromal readout, the direction consistent with this study's
## hypothesis. Excluding them makes the hypothesis harder to support. This is the
## opposite of the earlier narrow reading, which was the hypothesis-friendly one.
##
## k is untouched either way (Amendment 14: k is computed over the locked 152),
## and the Amendment 3 branch is BUILT under 152, 143 and 140 alike.
PAR_GENES  <- c("CRLF2", "CSF2RA", "IL3RA", "IL9R")
sex_any    <- sort(g152[grepl("chr[XY]", chrom)])
sex_only   <- sort(g152[grepl("^(chr)?[XY]$", chrom)])
par_seen   <- sort(setdiff(sex_any, sex_only))
if (!identical(par_seen, sort(PAR_GENES)))
  halt("the set of pseudoautosomal panel genes changed: observed ",
       paste(par_seen, collapse=", "), " vs the four this lock was built against (",
       paste(sort(PAR_GENES), collapse=", "), "). Rule 3.3's SCOPE is settled -- ",
       "Amendment 15 applies the literal reading, PAR genes excluded -- so this ",
       "halt means the ANNOTATION changed, not that the interpretation is open. ",
       "Re-derive the 140/143 counts and record them in a dated amendment before ",
       "relocking; do not simply widen this constant.")
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
## PAR_REGION: assign each PAR gene to PAR1 or PAR2 from its own chrX coordinate,
## recomputed at lock time.  The PAR boundaries are a property of the ASSEMBLY,
## not of any gene, so they are stated here as GRCh38 intervals and the assignment
## is derived -- never carried as a remembered list.  Amendment 15's reasoning
## turns on these four genes spanning BOTH regions, so the claim is verified
## rather than trusted.
GRCh38_PAR1_X <- c(10001L,     2781479L)
GRCh38_PAR2_X <- c(155701383L, 156030895L)
## ASSEMBLY GUARD.  The PAR intervals above are GRCh38-specific, and these objects
## carry NO usable assembly field -- seqlengths() is all NA and there is no genome
## string, so nothing in the file states the build.  Applying GRCh38 intervals to a
## GRCh37 annotation would misassign every PAR gene silently.  Pin it with a
## coordinate that differs between builds: MYC is chr8:127,735,434-127,742,951 in
## GRCh38 and chr8:128,748,315-128,753,680 in GRCh37, ~1 Mb apart.
local({
  i <- which(names(gene_start) == "MYC" & gene_chr == "chr8")
  if (length(i) != 1L) halt("cannot pin the assembly: MYC is not uniquely annotated on chr8")
  if (!(gene_start[i] == 127735434L && gene_end[i] == 127742951L))
    halt("annotation is not GRCh38: MYC observed at chr8:", gene_start[i], "-", gene_end[i],
         " but GRCh38 places it at chr8:127735434-127742951 (GRCh37: 128748315-128753680). ",
         "The PAR intervals below are GRCh38-specific and would misassign PAR1/PAR2 ",
         "on any other build. Do not lock against a different assembly.")
})
par_region <- vapply(PAR_GENES, function(g) {
  i <- which(names(gene_start) == g & gene_chr == "chrX")
  if (length(i) != 1L) return(NA_character_)
  s <- gene_start[i]; e <- gene_end[i]
  if (s >= GRCh38_PAR1_X[1] && e <= GRCh38_PAR1_X[2]) "PAR1"
  else if (s >= GRCh38_PAR2_X[1] && e <= GRCh38_PAR2_X[2]) "PAR2"
  else NA_character_
}, character(1))
if (any(is.na(par_region)))
  halt("a pseudoautosomal gene's chrX coordinates fall outside both GRCh38 PAR ",
       "intervals: ", paste(names(par_region)[is.na(par_region)], collapse=", "),
       ". The assembly may not be GRCh38; do not lock against a different one.")
## Amendment 15 states three PAR1 and one PAR2. Asserted, not assumed.
if (!identical(sort(names(par_region)[par_region == "PAR1"]), sort(c("CRLF2","CSF2RA","IL3RA"))) ||
    !identical(names(par_region)[par_region == "PAR2"], "IL9R"))
  halt("PAR region assignment disagrees with Amendment 15 (which records CRLF2, ",
       "CSF2RA, IL3RA as PAR1 and IL9R as PAR2): observed ",
       paste(sprintf("%s=%s", names(par_region), par_region), collapse=", "))

ex33 <- sex_any                       # BROAD -- Amendment 15, the literal reading
ex33_narrow <- sex_only               # retained for the sensitivity list

## ---- assemble, in the registered order 1 -> 2 -> 3 --------------------------
excl <- unique(c(ex31, ex2, ex33))
final <- setdiff(g152, excl)
reason <- function(g) paste(c(if(g %in% ex31)"3.1_absent_from_a_TCGA_cohort_annotation",
                              if(g %in% ex2) "2_undetectable_under_1pct_in_every_compartment",
                              if(g %in% ex33) if (g %in% PAR_GENES)
                                 "3.3_sex_chromosome_pseudoautosomal" else "3.3_sex_chromosome"),
                            collapse="+")

## ---- assertions against the Part A / Amendment 14 record --------------------
## Amendment 15: 7 genes qualify under rule 3.3 (3 chrX-only + 4 PAR); CRLF2 is
## ALSO a rule-2 exclusion, so the union is 12, not 13, and the primary set is 140.
stopifnot(length(ex31)==0L, length(ex2)==6L, length(ex33)==7L,
          length(excl)==12L, length(final)==140L)
if (!identical(sort(intersect(ex2, ex33)), "CRLF2"))
  halt("Amendment 15 records CRLF2 as the single gene caught by both rule 2 and ",
       "rule 3.3; observed overlap: ",
       if (length(intersect(ex2, ex33))) paste(sort(intersect(ex2, ex33)), collapse=", ") else "(none)")
## The sensitivity set is the narrow reading: 143.
final_narrow <- setdiff(g152, unique(c(ex31, ex2, ex33_narrow)))
stopifnot(length(final_narrow)==143L, all(final %in% final_narrow),
          setequal(setdiff(final_narrow, final), c("CSF2RA","IL3RA","IL9R")))
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
## K-DENOM-NA: which(), for the same reason as R2-NAINDEX above -- a bare logical
## with an NA injects an all-NA row rather than dropping it.  It happened to be
## caught downstream by the n_atlas guard, but with a message pointing at the
## wrong cause.
p  <- dm[which(dm$in_panel), ]
## NA dominance is EXPECTED and correct: a gene below the >=20-count evidence
## floor has no fraction, so A.d gives it no dominance call. What must hold is the
## CORRESPONDENCE -- NA exactly where the gene is not evaluable, never where it
## is. An NA on an evaluable gene would be a lost call that sum(na.rm=TRUE) would
## silently score as non-dominant, deflating k. (An earlier revision of this guard
## forbade NA outright and halted on the real data; that was my error, not the
## data's.)
if (any(is.na(p$dominant) & p$evaluable))
  halt("dominance matrix has NA dominance calls on EVALUABLE genes: ",
       paste(utils::head(p$gene[is.na(p$dominant) & p$evaluable], 5), collapse=", "),
       ". A missing call on an evaluable gene would be summed as non-dominant.")
if (any(!is.na(p$dominant) & !p$evaluable))
  halt("dominance matrix has a dominance call on a NON-evaluable gene: ",
       paste(utils::head(p$gene[!is.na(p$dominant) & !p$evaluable], 5), collapse=", "))
if (!all(table(p$gene) == 3L))
  halt("dominance matrix is not exactly one row per gene per atlas")
## K-DENOM: the k assertion re-derives from a separate committed Part A input, so
## it is not circular -- but it says nothing about the DENOMINATOR.  Pin it: the
## dominance matrix must cover exactly the panel, over exactly three atlases, with
## no NA dominance call standing in for a missing one.
n_atlas <- length(unique(p$atlas))
if (n_atlas != 3L) halt("dominance matrix covers ", n_atlas, " atlases, expected 3")
nd <- tapply(p$dominant, p$gene, function(x) sum(x, na.rm=TRUE))
if (!setequal(names(nd), g152))
  halt("dominance matrix in_panel rows do not cover the locked panel")
## K-DENOM-140: `final` is the 140-gene PRIMARY list, so a variable named k143
## computed over it was a mislabel -- it was being asserted against 43, the value
## Amendment 14 registers for the 143-gene list. The assertion passed only because
## the two coincide (the three extra genes are non-dominant). Both are now
## computed over their own denominators and their coincidence is asserted rather
## than relied on.
k152 <- sum(nd >= 2);       a152 <- sum(nd == n_atlas)
k140 <- sum(nd[names(nd) %in% final] >= 2)
a140 <- sum(nd[names(nd) %in% final] == n_atlas)
k143 <- sum(nd[names(nd) %in% final_narrow] >= 2)
a143 <- sum(nd[names(nd) %in% final_narrow] == n_atlas)
if (k140 != k143 || a140 != a143)
  halt("k differs between the primary 140 and sensitivity 143 lists (k140=", k140,
       ", k143=", k143, "; k_all3 ", a140, " vs ", a143, "). Amendment 15 records ",
       "the three PAR genes it removes as non-dominant, so the counts must agree; ",
       "if they do not, that claim is false and the amendment needs revisiting.")
if (k152 != 46L || k143 != 43L || k140 != 43L || a152 != 24L || a143 != 23L || a140 != 23L)
  halt("Amendment 14 records k=46 over the panel and 43 over the final list, ",
       "k_all3=24 and 23; Amendment 15 leaves both unchanged at 140. observed ",
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
## OUT-PRIMARY-NOCOL / LIST-IDENTITY: the two lists are built by ONE function so
## they are column-identical by construction.  Previously the primary lacked
## in_primary_140, and in R `d$in_primary_140` on a frame without that column
## returns NULL while `d[d$in_primary_140, ]` returns ZERO ROWS with no error --
## a Part B script filtering the wrong file would silently score nothing.
## Each frame also carries list_id and n_genes, so a reader can assert WHICH list
## it holds from the data rather than trusting a filename that has already changed
## meaning once (final_gene_list.csv held 143 before Amendment 15 and holds 140 now).
make_list <- function(genes, id) data.frame(
  list_id             = id,
  n_genes             = length(genes),
  gene                = genes,
  in_locked_panel_152 = TRUE,
  origin_six          = genes %in% ORIGIN6,
  chromosome          = unname(chrom[genes]),
  pseudoautosomal     = genes %in% PAR_GENES,
  par_region          = unname(ifelse(genes %in% names(par_region), par_region[genes], NA_character_)),
  in_primary_140      = genes %in% final,
  max_detection_pct   = round(100*unname(mx[genes]),4),
  n_atlases_evaluable = unname(n_atl[genes]),
  n_annotation_rows   = unname(nrow_annot[genes]),
  stringsAsFactors    = FALSE)

fl    <- make_list(final,        "primary_140")
fl143 <- make_list(final_narrow, "sensitivity_143")
stopifnot(identical(names(fl), names(fl143)),
          sum(fl$in_primary_140) == 140L, sum(fl143$in_primary_140) == 140L)
write.csv(fl,    file.path(OUT,"final_gene_list_140.csv"), row.names=FALSE)  # PRIMARY
write.csv(fl,    file.path(OUT,"final_gene_list.csv"),     row.names=FALSE)  # alias -> primary
write.csv(fl143, file.path(OUT,"final_gene_list_143.csv"), row.names=FALSE)  # SENSITIVITY
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
