# STAT3 target panel — prespecification

**Status:** LOCKED.
**Locked on:** 2026-07-31  **By:** Sean GP Lee
**Commit hash at lock:** ac9c5e0 (commit adding Amendments 1–4)

Once locked, this file is not edited. Changes go in an amendments section at the
bottom, dated, with a reason. That record is what makes the panel defensible.

---

## 1. Why this document exists

The paper's claim is that a clinically deployed STAT3 score is substantially a
non-tumour readout. If genes were chosen after seeing which ones came out
epithelial, the claim would be circular in exactly the way the paper accuses
others of being. So the panel is defined by external criteria, in advance, and
compartment fraction is never an inclusion or exclusion criterion.

## 2. Inclusion criteria

A gene enters the panel if it meets **A**, or **B**, or both.

**A — Pathway membership.** Member of MSigDB
`HALLMARK_IL6_JAK_STAT3_SIGNALING`. Record MSigDB release version.

**B — Direct transcriptional target.** Supported as a direct STAT3 target by at
least **two independent** lines of evidence from the following, at default
thresholds:
  - ChIP-seq derived TF-target databases (e.g. ChEA3, ENCODE TF targets)
  - Curated TF-target databases (e.g. TRRUST v2)
  - A promoter-level STAT3 binding site reported in primary literature

Two sources querying the same underlying ChIP-seq experiment count as one.
Record the exact resource versions and query date.

## 3. Exclusion criteria

Applied in this order, and none of them involve compartment fraction:

1. Not present in the gene annotation of every TCGA cohort in the study.
2. Not detected in the single-cell atlases at a minimum of 1% of cells in at
   least one annotated compartment (undetectable genes carry no compartment
   information either way).
3. Sex-chromosome genes, since cohorts differ in sex composition.

Genes excluded under rule 2 are **reported**, not silently dropped — a canonical
STAT3 target that is undetectable in tumour tissue is itself informative.

## 4. The original six

`SOCS3, BCL2, MYC, MMP9, HGF, IL6` are flagged `origin_score = TRUE` and
analysed as a highlighted subset. They receive no special treatment in panel
construction and are not guaranteed inclusion — a gene from the original score
that fails section 2 is dropped, and that is reported.

## 5. Prespecified decision rule after the compartment sweep

Let *k* = number of panel genes that are epithelial-dominant (abundance-weighted
epithelial fraction > 0.50) across the full 30–70% purity band, in **at least
two** tissue-matched atlases.

| *k*   | What happens to the decomposition |
|-------|-----------------------------------|
| ≥ 8   | Epithelial subscore is built; subscore survival models run as a **secondary** analysis |
| 3–7   | Subscore reported as **exploratory** only; no formal survival inference on it |
| < 3   | Decomposition is **descriptive only**; reported as a finding, no subscore |

This rule is fixed now so the result cannot be argued into whichever branch
looks better later. All three outcomes are publishable.

## 6. The primary endpoint does not depend on any of this

Primary analysis is the purity-adjusted Cox model on the **intact** score, per
cohort and meta-analysed. It requires no compartment fractions and cannot
degenerate. The compartment work explains whatever attenuation is observed; it
does not carry the headline claim.

Practical consequence: cohorts with no usable single-cell atlas (ESCA, possibly
CHOL) still contribute fully to the primary result.

## 7. Provenance to capture at lock time

- MSigDB release version
- ChEA3 / TRRUST / other resource versions and query date
- Final gene count, and count entering via A only, B only, both
- `sessionInfo()` and `renv.lock`
- Every gene's inclusion route, in the output table — one row per gene

## 8. Registration timing

Preregister after the panel is locked and after the compartment sweep is
specified, but **before** any survival model is fitted. Disclose in the
registration that HPA and GSE125449 compartment fractions for the original six
genes were inspected beforehand — that inspection motivated the study and hiding
it would be worse than declaring it.

---

## Amendments

### Amendment 1 — 2026-07-31

Section 2, criterion B, is amended from "at least two independent lines of
evidence" to:

  B — Direct transcriptional target. Gene must satisfy BOTH:
      (i)  reported as a STAT3 target in a ChIP-seq-derived resource
           (ChEA3, at default thresholds), AND
      (ii) independently curated as a STAT3 target in TRRUST v2.

Reason: as originally written, criterion B admitted 1,387 genes (1,448 in union
with criterion A). ChIP-seq occupancy is far broader than direct transcriptional
regulation, and an any-two-sources rule counts multiple re-slices of the same
occupancy data as independent. The resulting panel would no longer be the class
of object under study — a compact clinically deployed activity score — and would
render the §5 thresholds vacuous.

Made BLIND to all compartment and survival results. No compartment fraction, no
purity sweep, and no outcome data had been examined at the time of this
amendment. Panel size at amendment: not yet computed.

### Amendment 2 — 2026-07-31

Criterion B(i) is clarified to require ChIP-seq evidence from a HUMAN
experiment. Mouse ChIP-seq does not satisfy B(i).

Reason: cross-species conservation of promoter occupancy is unreliable, and the
study's cohorts are human. Applied at implementation of Amendment 1; recorded
here as a formal amendment because it was not in the adopted text.

Disclosure: this decision was NOT blind. Pilot compartment results in
feasibility_assessment.md were already known, and MMP9 and HGF — the two genes
excluded by this clarification — were the two most stroma/myeloid-dominant genes
in that pilot. The restriction removes the clearest illustrations of the paper's
own thesis and is therefore conservative with respect to it. No compartment
fraction or survival result for the amended 152-gene panel had been computed at
the time of this amendment.

MMP9, HGF and BCL2 are reported as a labelled non-qualifying subset per §4.

### Amendment 3 — 2026-07-31

Section 5's thresholds are restated as proportions of the final panel. They were
written anticipating a 30–50 gene panel; the panel is 152, which would make
k >= 8 a 5% bar rather than the ~20% intended.

Amended rule, where k = panel genes that are epithelial-dominant
(abundance-weighted malignant-epithelial fraction > 0.50) across the full
30–70% purity band in at least two tissue-matched atlases:

  k >= 20% of final panel AND k >= 8   -> epithelial subscore built; subscore
                                          survival models run as SECONDARY
  5% <= k < 20%                        -> subscore reported as EXPLORATORY only
  k < 5%                               -> decomposition DESCRIPTIVE only

Made blind: no compartment fraction for the 152-gene panel had been computed.

### Amendment 4 — 2026-07-31

The compartment estimand is changed from malignant-epithelial fraction to
EPITHELIAL fraction computed on tumour-channel samples only.

Primary definition: abundance-weighted contribution to pseudobulk from the
lineage-level EPITHELIAL compartment, computed using only samples designated by
the source atlas as tumour (excluding adjacent-normal, peritoneal and
normal-donor samples). Non-epithelial compartments are grouped as
fibroblast/stromal, myeloid, lymphoid, endothelial and other, per each atlas's
own level-1 labels.

Sensitivity analysis: where the source atlas supplies CNV-based per-cell
malignancy calls (GSE125449, GSE183904), the fraction is recomputed restricted
to malignant cells. Peng (PDAC) is reported separately because its malignancy
call is a cluster-level ductal-subtype proxy rather than a per-cell threshold.
Agreement between primary and sensitivity definitions is reported.

Reason: the five compartment atlases infer malignancy by three incompatible
methods (inferCNV per-cell thresholds in Ma; inferCNV with unstated cutoff in
Kumar; a trained classifier in Pelka, where ~11% of malignant-called cells show
no substantial CNV difference; a cluster-level ductal proxy in Peng; and no
malignancy inference at all in Steele — see
output/atlas_malignant_annotation_audit.csv). A cross-tissue comparison whose
numerator is defined differently in each tissue would not measure the same
quantity. Lineage-level epithelial labels are available and comparably defined
in all five. Re-deriving malignancy uniformly was rejected because it would make
every malignancy call this project's own construction on data of differing
platform and depth.

This supersedes the phrase "abundance-weighted malignant-epithelial fraction" in
Amendment 3; k is now defined on the primary estimand above.

Direction of bias: residual normal epithelium in tumour samples inflates the
epithelial fraction, biasing AGAINST the paper's thesis that the score is
substantially non-epithelial. The change is therefore conservative.

Consequence: GSE155698 (Steele) becomes usable at lineage level, giving PDAC two
tissue-matched atlases as Amendment 3 requires.

Made blind: no compartment fraction for the 152-gene panel had been computed.

### Amendment 5 — 2026-07-31

Amendment 3's replication requirement, "in at least two tissue-matched atlases",
is amended to: "in at least two of the five compartment atlases, of any GI
tissue."

Reason: only pancreatic has two atlases (GSE155698, Peng). Liver, colorectal and
gastric have one each and oesophageal none, so under a tissue-matched reading k
could only ever be driven by pancreatic data. The requirement's purpose is
robustness against a single dataset's dissociation protocol or annotation scheme
determining the call — a threat that is atlas-specific, not tissue-specific. The
subscore, if built, is applied across seven TCGA cohorts, so cross-tissue
consistency is the relevant property.

Reported alongside k: the per-atlas dominance matrix, a per-tissue breakdown, and
a strict count of genes dominant in ALL five atlases.

Rejected alternative: adding the Sci Data 2026 integrated atlas to supply second
atlases for gastric and colorectal. It is an integration of constituent GEO
series, and if GSE183904 or GSE178341 are among them the two atlases would share
cells, making the replication requirement appear satisfied while providing no
independent evidence.

Direction of bias: this is looser than the tissue-matched reading, yielding a
larger k and making the epithelial-subscore branch more reachable — i.e. biasing
AGAINST the paper's thesis that the score is substantially non-epithelial.
Conservative.

Made blind: no compartment fraction for the 152-gene panel had been computed.

### Amendment 6 — 2026-07-31

Disclosure and sensitivity analysis for the "epithelial-dominant across the full
30–70% band" criterion. The primary rule is UNCHANGED.

Because f(pi) = pi*A / (pi*A + (1-pi)*B), where A is epithelial intensity per
cell and B the abundance-weighted non-epithelial intensity, dominance at pi is
equivalent to A/B > (1-pi)/pi. Dominance across the full band is therefore
evaluated entirely at the lower boundary and requires A/B > 2.33, whereas
dominance at typical tumour purity (pi = 0.50) requires only A/B > 1.

The full-band bar is thus substantially stricter, yields a smaller k, and makes
the descriptive-only branch more likely — the branch most consistent with the
paper's thesis. This is the one criterion in the design that leans toward the
hypothesis, and it is disclosed as such.

Additional prespecified reporting: k_50, defined identically but evaluated at
pi = 0.50 only, reported alongside primary k with its bootstrap interval. The
branch decision uses primary k. If primary k and k_50 fall in different branches,
that fact is stated in the paper.

Made blind: no compartment fraction for the 152-gene panel had been computed.

### Amendment 7 — 2026-07-31

Endpoint designation is restated as a single deterministic rule, applied
identically to all seven cohorts:

  Primary endpoint = overall survival (OS) wherever the TCGA-CDR (Liu et al.
  2018, Table 3) marks OS usable without caution. Where OS carries a caution
  mark, the primary endpoint is the CDR-preferred alternative for that cancer
  type. The non-primary endpoint is reported as a prespecified sensitivity
  analysis in every cohort.

Reason: the v1.0 provisional table assigned COAD to PFI on the higher event
count while assigning LIHC to OS despite PFI having more events — two
contradictory tiebreaks, i.e. discretion rather than a rule. The v1.0 COAD
justification ("CDR notes OS underpowered in colon") was also factually wrong;
Table 3 marks COAD OS as usable without caution. Correcting the fact while
retaining the conclusion under a substituted justification is not defensible,
so the rule is restated and re-derived from scratch.

Rationale for defaulting to OS: it is objective mortality, not subject to
assessment or ascertainment variation, and departing from it only where the CDR
itself flags a problem removes all analyst discretion from the choice.

Made blind: no score has been computed, no model fitted, and no patient-level
outcome data merged with expression. Only the Table 3 recommendation matrix was
consulted.

### Amendment 8 — 2026-08-01

Prespecified model-degradation rule for small cohorts.

Events per variable (EPV) is computed per cohort as (events on that cohort's
primary endpoint) / (number of estimated parameters in the model being fitted).
Applied per model, not per cohort, so a cohort may support M2 but not M4:

  EPV >= 10  -> model fitted and entered into the meta-analysis as specified
  5 <= EPV < 10 -> model fitted, entered into the meta-analysis, and flagged in
                   the results table as low-EPV; a leave-one-out meta-analysis
                   omitting every low-EPV cohort is reported as a prespecified
                   sensitivity
  EPV < 5    -> model NOT fitted for that cohort; the cohort is reported
                descriptively and excluded from the meta-analysis for that model
                only

Cohort-level inclusion, superseding the ad hoc treatment of CHOL in B.i: a cohort
whose TCGA-CDR Table 3 explanation states the sample size is too small for all
endpoints is reported descriptively and excluded from the meta-analysis as a
weighted stratum, regardless of EPV. This currently applies to TCGA-CHOL and is
stated as a rule so it does not read as a post hoc judgement about one cohort.

Reason: TCGA-READ has 26 OS and 39 PFI events against four covariates in M2 and
six in M4, so it may not support the full model sequence. Deciding that after
seeing the fits would reintroduce exactly the discretion Amendment 7 removed.

Made blind: no score computed, no model fitted, no outcome data merged with
expression. Event counts are the CDR's published 2018 snapshot, read from
Table 3.

### Amendment 9 — 2026-08-01

GSE155698 (Steele et al., pancreatic) is removed from the compartment atlas set.
The set is reduced from five to four: GSE125449 (liver/biliary), GSE178341
(colorectal), GSE183904 (gastric), and Peng (pancreatic).

Reason: Amendment 4 states that lineage-level epithelial labels "are available
and comparably defined in all five" atlases and that GSE155698 "becomes usable at
lineage level." That is factually incorrect for this atlas. GSE155698_RAW.tar
contains 41 per-sample archives holding only CellRanger output (barcodes,
features, matrix), with no cell-type annotation table, no metadata file, and no
annotation column in any sample; the four samples carrying extra files carry
only format variants of the same matrices. The GEO series has no other
supplementary file. The annotations exist in Steele et al.'s publication but
were not deposited.

Deriving cell-type labels de novo was rejected: it would make pancreatic the only
tissue whose compartment labels are this project's own construction, on data of
different platform and depth from the other atlases, reintroducing precisely the
cross-atlas method asymmetry Amendment 4 was written to eliminate. Substituting
a different pancreatic atlas was also rejected, since selecting a replacement
dataset after registration is a researcher degree of freedom the amendment
record cannot neutralise.

Consequent restatements:
- Amendment 5's replication requirement becomes "in at least two of the four
  compartment atlases, of any GI tissue."
- k_all5 becomes k_all4 (dominance in all four atlases).
- k_eval4 is replaced by TWO reported quantities, because no single threshold on
  four atlases is the exact analogue of "at least four of five" (80%):
  k_eval3 = k restricted to genes evaluable in at least three of four atlases;
  k_evalall = k restricted to genes evaluable in all four atlases. Both are
  reported alongside k. k_evalall is the stricter analogue and preserves the
  original's intent; k_eval3 is retained because k_evalall may prove small
  enough to be uninformative. Neither gates the branch decision, which uses
  primary k. Reporting both rather than choosing one avoids a threshold selected
  after the atlas set changed.
- Amendment 4's sentence "GSE155698 (Steele) becomes usable at lineage level,
  giving PDAC two tissue-matched atlases as Amendment 3 requires" is withdrawn
  as factually incorrect.

Each of the four tissues now has exactly one atlas. The two-atlas replication
requirement is therefore necessarily satisfied across two different tissues,
which is a stronger generality claim than two atlases of one tissue would have
been.

DIRECTION OF BIAS — this amendment is NOT conservative. Removing an atlas can
only reduce k or leave it unchanged: a gene dominant in exactly two atlases, one
of which was GSE155698, now falls to one and is no longer counted. A smaller k
makes the descriptive-only branch more likely, and that is the branch most
consistent with this study's hypothesis. This is the one amendment in the record
that favours the hypothesis. It is unavoidable — no handling of a missing
annotation is neutral — and it is disclosed rather than absorbed. The full
per-atlas dominance matrix is reported so readers can see which genes were
affected.

Made blind: no compartment fraction for the 152-gene panel had been computed.
The only atlas content inspected was the file manifest of GSE155698_RAW.tar,
Peng's cell-type vocabulary, and Peng's tumour/normal sample counts. No
expression values entered any computation.

### Amendment 10 — 2026-08-01

GSE183904 (Kumar et al., gastric) is removed from the compartment atlas set. The
set is reduced to three: GSE125449 (liver/biliary), GSE178341 (colorectal), and
Peng (pancreatic). No gastric or oesophageal atlas remains.

Reason: direct inspection of GSE183904_RAW.tar shows 40 flat .csv.gz members, one
per GSM, each a gene x cell count matrix (26,572 gene rows, barcode column
headers). No annotation row, no metadata table, no cluster file, in any sample.
The GEO series carries no supplementary file beyond the RAW tar and its 40
GSM-level equivalents. As with GSE155698, the annotations exist in the source
publication but were not deposited.

This is the second atlas to fail the same premise. Amendment 4's statement that
lineage-level labels are "available and comparably defined in all five" atlases
was derived from the source publications' methods sections rather than from the
deposits themselves, and was wrong for two of five. The error and its detection
are reported in the paper rather than silently corrected: annotation availability
in public single-cell deposits is materially worse than the corresponding
publications imply, which is itself a finding relevant to any study of this design.

Consequent restatements:
- Amendment 5's replication requirement becomes "in at least two of the three
  compartment atlases, of any GI tissue."
- k_all4 becomes k_all3 (dominance in all three atlases).
- k_eval3 and k_evalall collapse to a single quantity at three atlases and are
  replaced by k_evalall = k restricted to genes evaluable in all three atlases.
- Compartment coverage is now liver/biliary, colorectal and pancreatic only.
  TCGA-STAD, TCGA-ESCA and TCGA-READ have no tissue-matched atlas. This does not
  affect Part B, which does not depend on the decomposition, and those cohorts
  contribute fully to the primary analysis.

Rejected alternative: substituting the Sci Data 2026 integrated gastrointestinal
atlas for gastric coverage. Amendment 5's two grounds for rejecting it no longer
apply — its lack of malignant-cell annotation ceased to matter under Amendment 4's
lineage-level estimand, and the risk of shared cells with GSE183904 is moot now
that GSE183904 is excluded. It was nonetheless rejected because selecting a
replacement dataset after registration is a researcher degree of freedom the
amendment record cannot neutralise, on the same reasoning applied to Steele in
Amendment 9. The loss of gastric coverage is reported as a limitation.

DIRECTION OF BIAS — this amendment is NOT conservative, on two counts. First,
removing an atlas can only reduce k or leave it unchanged. Second, the
replication requirement rises from two of five (40%) through two of four (50%) to
two of three (67%) without ever being deliberately restated as a proportion; the
bar a gene must clear is now substantially higher than at registration. Both
effects shrink k and make the descriptive-only branch more likely, which is the
branch most consistent with this study's hypothesis. Retaining "at least two" is
the only coherent replication requirement at three atlases, so the escalation is
unavoidable, but it is disclosed rather than absorbed. The full per-atlas
dominance matrix is reported so readers can see exactly which genes each removal
affected.

Made blind: no compartment fraction for the 152-gene panel had been computed. The
only atlas content inspected was file manifests, matrix dimensions, gene-row
names, and cell-type vocabularies. No expression values entered any computation.

Also corrected in this amendment, from direct inspection rather than metadata:
GSE178341's expected tumour channel count is 129 GSMs in GEO but 128 in the
deposited matrices. Channel C144_T_1_1_12_c1_v2 (GSM5388094, CD45pCD3nCD19nMACS
fraction) contributes zero cells, consistent with a QC failure during processing;
patient C144 retains two other tumour channels. Both counts are asserted. The
expected patient count remains 62 (34 MMRd, 28 MMRp), confirmed identical in the
metatables and GEO title parsing, matching Pelka's published cohort.

### Amendment 11 — 2026-08-01

For the Peng pancreatic atlas, compartment mapping uses the celltype1 annotation
level rather than celltype0.

Reason: Amendment 4 specifies six compartments (epithelial, fibroblast/stromal,
myeloid, lymphoid, endothelial, other) mapped from "each atlas's own level-1
labels." Peng's celltype0 collapses myeloid and lymphoid into a single
"hematopoietic" category and therefore cannot express Amendment 4's compartment
scheme. celltype1 separates myeloid leukocyte, T lineage and B lineage, and maps
onto the six compartments directly.

DIRECTION OF BIAS: none, and this is provable rather than argued. The estimand is
the epithelial contribution to pseudobulk, i.e. epithelial counts divided by total
counts. Subdividing the non-epithelial mass changes neither numerator nor
denominator, so f(pi) is numerically identical under celltype0 and celltype1 for
every gene at every grid point, and k is unaffected. Only the reported breakdown
of non-epithelial compartments differs.

Made blind: no compartment fraction for the 152-gene panel had been computed.

### Amendment 12 — 2026-08-01

For GSE125449, the two deposited sets (Set1, 20,124 gene rows; Set2, 19,572) are
combined on the INTERSECTION of their gene universes. Genes present in only one
set are excluded from that atlas and reported.

Reason: analysis_plan.md v1.5 specifies combining the sets but does not state how
to reconcile their differing gene universes, and load_GSE125449's guard correctly
halts rather than cbind-ing mismatched row spaces. Of the three available
reconciliations, only intersection avoids inventing a value. Union with zero-fill
would record a gene as measured at zero when it was not measured, violating A.d's
rule that absence is NA and never 0. Union with NA-fill would give Set1 and Set2
different denominators, so their pseudobulk fractions would not be comparable
across the very sets being combined. Intersection keeps every retained value a
measured one.

Affected panel genes, excluded from GSE125449 only: CCL7, CRLF2, CSF2, IL9R,
ITGB3, LEP. Panel coverage in this atlas becomes 143 of 152. No origin-six gene is
affected. These genes remain fully evaluable in GSE178341 and Peng, so each can
still reach k on those two atlases under Amendment 5's two-of-three requirement.

DIRECTION OF BIAS: this amendment is NOT conservative. Six genes lose one of their
three possible atlases, so each must now be dominant in both remaining atlases
rather than any two of three. That can only reduce k or leave it unchanged, making
the descriptive-only branch more likely — the branch most consistent with this
study's hypothesis. This is the third such effect, after Amendments 9 and 10, and
is disclosed on the same footing. The evaluability distribution and the per-atlas
dominance matrix are reported so readers can see exactly which genes each
reduction affected, and the six names above are reported explicitly.

Made blind: no compartment fraction for the 152-gene panel had been computed. The
only content inspected was gene row names and set dimensions.
