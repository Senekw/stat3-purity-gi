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