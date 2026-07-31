# STAT3 target panel — prespecification

**Status:** DRAFT — lock before running any compartment analysis.
**Locked on:** ____________  **By:** ____________
**Commit hash at lock:** ____________

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