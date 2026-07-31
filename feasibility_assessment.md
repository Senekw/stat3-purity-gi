# Feasibility and Novelty Assessment
## Proposed study: compartment decomposition of a six-gene pSTAT3 activity score in GI cancers

**Prepared:** 30 July 2026
**Score under study:** SOCS3, BCL2, MYC, MMP9, HGF, IL6 (mean log2(expr+1)), from the
ESMO Asia 2026 abstract "A 6-Gene pSTAT3 Transcriptomic Score Identifies an
Immunosuppressive, Chemotherapy-Resistant Phenotype and Predicts Poor Survival in
Biliary Tract Cancer."
**Related prior work by the same author:** CCRC-D-26-00363 (napabucasin
complete-enumeration meta-analysis; pooled OS-HR 1.028, CO.23 pSTAT3 non-replication
P = 0.006).

---

## 1. Verdict

**Plausible and worth doing, with one required reframing.** The mechanics are sound
and the data are all public and open-tier. The novelty is real but narrower than
the framing in the planning notes implies: "prognostic signal in bulk tumour
transcriptomes is substantially stromal" has been established since 2015 and cannot
be the headline. What is genuinely unoccupied is the same question asked of *this
score*, whose clinical relevance is anchored by a failed Phase 3 programme.

Confidence that the analysis is executable as specified: **high**.
Confidence that it is publishable in a specialty journal: **moderate-to-high**,
conditional on the reframing in §4.

---

## 2. The novelty check

Systematic search of Europe PMC (30 July 2026), eight query families across
signature/purity/single-cell/STAT3 vocabulary.

### 2.1 What is already established — and constrains you

| Prior work | What it established | Consequence |
|---|---|---|
| Isella 2015 (~849 cites); Calon 2015 (~495) | Stromal/CAF programs carry much of the prognostic signal in bulk CRC | "Signal is stromal" is settled. Not your discovery. |
| Guinney 2015 (~4,020 cites) | CMS4 mesenchymal = stroma-rich, poor prognosis | Reviewers will ask whether your stromal subscore is CMS4 relabelled. Test it. |
| Venet 2011 (~470 cites) | Most random signatures associate with breast cancer outcome | Your null-signature analysis *is* this method. Cite it; do not claim it. |
| Aran 2015 (~962 cites) | Pan-cancer purity; source of CPE | Purity-as-confounder is generic. Your claim must be score-specific. |
| Stroma admixture in serous ovarian cancer, CEBP 2020 | The same design, different disease | Design is publishable; "their paper in GI" is the criticism to pre-empt. |
| STAT3 prognostic meta-analyses (Oncotarget 2016; PLoS One 2015) | pSTAT3 IHC associates with poor prognosis in digestive cancers | IHC is cell-resolved. A bulk-stromal finding does **not** contradict them — say so before a reviewer does. |

### 2.2 What is genuinely open

Two searches came back essentially empty:

- **STAT3 + tumour purity**: 11 records in Europe PMC, none asking whether a STAT3
  activity score is a purity artefact.
- **Prespecified epithelial/stromal decomposition of a named clinical score into
  subscores**: 1 record, unrelated.

Nobody has asked whether a STAT3 transcriptional activity score — the class of
biomarker that selected patients in the napabucasin programme — is measuring tumour
cells or measuring stroma. That is a real gap, and it is the right size for the
paper you describe.

---

## 3. Empirical feasibility check on the lead figure

Rather than assert that step 1 will work, it was run. Human Protein Atlas
single-cell consensus data (cell-type–enriched nCPM) were pulled for all six score
genes and collapsed onto epithelial / fibroblast-stromal / myeloid / lymphoid
compartments.

| Gene | Epithelial | Fibroblast/stromal | Myeloid | Lymphoid | Other |
|---|---|---|---|---|---|
| *SOCS3* | 0.17 | **0.62** | 0.14 | 0.00 | 0.07 |
| *BCL2*  | 0.20 | 0.13 | 0.00 | 0.17 | 0.51 |
| *MYC*   | **0.53** | 0.40 | 0.00 | 0.00 | 0.06 |
| *MMP9*  | 0.00 | 0.00 | **1.00** | 0.00 | 0.00 |
| *HGF*   | 0.00 | **0.71** | 0.13 | 0.00 | 0.16 |
| *IL6*   | 0.38 | **0.59** | 0.00 | 0.00 | 0.04 |

Under a prespecified 50% epithelial threshold, **one of six genes (*MYC*) is
epithelial-restricted.** *MMP9* is 100% myeloid (neutrophils and neutrophil
progenitors). *HGF* is 71% stromal, driven by hepatic stellate cells and pericytes.
*SOCS3* and *IL6* are majority stromal.

Three consequences you should absorb before committing:

1. **The lead figure exists and is striking.** This is a genuine result, not a
   hoped-for one.
2. **Your epithelial subscore will be one gene.** A one-gene "epithelial subscore"
   that is *MYC* is not a STAT3 score — it is a MYC score, and MYC is prognostic in
   nearly everything. Your prespecified threshold, applied honestly, degenerates the
   primary comparison. **Fix this before preregistering**, not after: either use a
   continuous purity-weighted decomposition instead of a hard split, or widen to a
   canonical STAT3 target panel (~30-50 genes) so the epithelial arm has enough
   genes to be a score. The wider panel is the better fix and the planning notes
   already gesture at it.
3. **These are normal-tissue atlas fractions, not tumour fractions.** HPA is a
   defensible screen and a legitimate sensitivity analysis, but the paper needs
   tumour-derived compartment fractions from GI tumour atlases. Treat the numbers
   above as a feasibility preview, not as Figure 1.

---

## 4. The reframing that makes this publishable

Do not write "STAT3 scores are stromal." Write:

> *A biomarker that has already been used to select patients in a Phase 3 trial
> programme is, on inspection of its constituent genes, largely a readout of
> non-tumour compartments — which offers a concrete, testable explanation for why
> the CO.23 pSTAT3 signal did not replicate.*

That framing does three things the generic version does not. It makes the paper a
sequel to your own meta-analysis rather than a re-run of Isella 2015. It gives the
finding a clinical consequence instead of a methodological observation. And it
explains a specific published non-replication, which is a far more defensible
contribution than another demonstration that stroma is prognostic.

It also correctly bounds the claim. Your own meta-analysis found pSTAT3 behaving
prognostically but not predictively across three trials; a stromal-dominant score is
exactly what a prognostic-not-predictive marker looks like. The two papers
interlock.

**A caution on self-citation.** Your ESMO abstract reports this score as a validated
prognostic biomarker. This paper argues it is substantially a stromal readout. That
is legitimate scientific self-correction and reviewers respect it — but it must be
explicit in the discussion, not left for someone else to notice. Frame it as
characterising the score you built, which is a strength.

---

## 5. Data availability — verified

| Resource | Status |
|---|---|
| TCGA COAD/READ/STAD/ESCA/PAAD/LIHC/CHOL | Open tier, no application. CHOL n≈36: descriptive only. |
| Survival endpoints | Use the curated TCGA Clinical Data Resource (Liu, Cell 2018). PFI for COAD/READ. |
| Purity | ESTIMATE + Aran CPE consensus table. Two sources = sensitivity analysis. |
| GSE39582 (CRC), ACRG (gastric), ICGC PACA | Open. |
| **Dong iCCA proteogenomics (n=262)** | Confirmed: Cancer Cell 2022, PMID 34971568. Raw WES/RNA-seq/proteome/phosphoproteome in the **BioSino NODE database, accession OEP001105**; sample annotation and processed/normalised data in **Table S1** of the paper. |

The iCCA answer to your open question: processed data are in the supplementary
tables, raw data in NODE OEP001105 — not a CPTAC portal deposit. Note the paper is
paywalled (no PMC copy), so budget for access to Table S1.

One caution on the iCCA validation. A measured STAT3 pY705 phosphosite in 262
patients is a genuinely strong anchor, but it is measured in **bulk tissue**, so it
carries the same compartment ambiguity as the mRNA score. It validates that your
score tracks STAT3 phosphorylation; it does not by itself establish which cells the
phosphorylation is in. State that explicitly.

---

## 6. Recommended changes to the analysis plan

1. **Widen the gene panel.** Run the compartment analysis on a canonical STAT3
   target set (~30-50 genes), reporting your six as a highlighted subset. This
   rescues the epithelial subscore from being one gene and turns the paper from
   "my score" into "this class of score."
2. **Replace the hard 50% threshold** with a continuous compartment-weighted
   decomposition, retaining the 50% split as a prespecified sensitivity analysis.
   Prespecify both.
3. **Add a CMS orthogonality analysis** in the CRC cohorts. One panel; pre-empts the
   most likely reviewer objection.
4. **Add pSTAT3 IHC reconciliation** to the discussion: bulk-stromal does not
   contradict cell-resolved IHC positivity.
5. **Preregister after step 1, before step 3.** You need the compartment fractions to
   set a sensible threshold, but you must fix the survival models before touching
   survival data. Register the analysis plan with the fractions already in hand and
   say so — that is honest and defensible, and it avoids repeating the post-hoc
   registration timing you disclosed in CCRC-D-26-00363.
6. **Expect the mushy outcome.** Attenuation-but-not-abolition with cohort
   heterogeneity remains the single most likely result. Under the reframing in §4 it
   is still a publishable paper, because the explanatory target is the
   non-replication, not a clean verdict.

---

## 7. Realistic assessment of venue

The analysis as scoped is a specialty-journal paper, not a high-impact one — the
methodological core is established and the contribution is a well-executed
application to a clinically anchored target. That is consistent with what you asked
for. The napabucasin link is what lifts it above a routine bioinformatics signature
paper, so lead with it.

Timeline of two to three months is realistic for someone already fluent in the
tooling, assuming tumour single-cell atlas selection does not stall.

---

# Addendum — estimand correction (30 July 2026)

## The critique was right about the metric

HPA cell-type–enriched nCPM is an **unweighted cluster share**: it treats every
cell type as equally abundant. Bulk contribution is expression × abundance. That
metric systematically inflates rare compartments and deflates epithelium, so the
original Figure was computing the wrong quantity.

Corrected as recommended: abundance-weighted pseudobulk (summed raw UMI counts per
compartment) computed in **GSE125449** (Ma et al., 19 liver-cancer tumours; 9,946
annotated cells; 10 iCCA and 9 HCC), using the authors' own labels — Malignant cell,
HPC-like, CAF, TEC, TAM, T cell, B cell.

## The numbers moved. The conclusion did not.

Epithelial share of each gene's bulk contribution, four estimands:

| Gene | HPA cluster share (normal) | GSE125449 HCC | GSE125449 iCCA | iCCA reweighted to 65% purity |
|---|---|---|---|---|
| *SOCS3* | 0.17 | 0.03 | 0.15 | 0.14 |
| *BCL2*  | 0.20 | 0.32 | 0.15 | 0.31 |
| *MYC*   | **0.53** | 0.11 | **0.67** | **0.69** |
| *MMP9*  | 0.00 | 0.03 | 0.12 | 0.26 |
| *HGF*   | 0.00 | 0.17 | 0.03 | 0.03 |
| *IL6*   | 0.38 | 0.01 | 0.12 | 0.23 |

**No gene passes the 50% epithelial threshold under all four estimands. Five of six
pass under none.** *MYC* is the only gene that clears 50% in tumour data at
realistic purity, which is what the original assessment concluded — so the
degeneracy warning stands, but it now rests on the correct quantity.

Specifically on the predicted threshold crossings: the correction moved *SOCS3*
(0.17 → 0.15), *IL6* (0.38 → 0.12), and *MMP9* (0.00 → 0.12) — but **downward or
flat**, not across the 50% line. The abundance weighting does inflate the epithelial
share relative to a pure cluster share (MMP9 goes from 0% to 12-26%), but not nearly
enough to rescue any gene. The prediction that abundance weighting would move the
53/59/62% genes across the line was directionally reasonable and empirically wrong.

*MMP9* remains substantially myeloid (69% of iCCA pseudobulk), as anticipated.

## Purity sensitivity, computed explicitly

Because compartment fraction depends on assumed purity, the analysis was run as a
sweep rather than at a single point. Compartment-specific per-cell intensities were
computed CP10K-normalised (composition-free), then reweighted to bulk compositions
from 20% to 95% purity.

Purity at which each gene's epithelial share first reaches 50%:

| Gene | Crossing purity |
|---|---|
| *MYC*   | 46% |
| *BCL2*  | 81% |
| *MMP9*  | 85% |
| *IL6*   | 86% |
| *SOCS3* | 92% |
| *HGF*   | never within 20–95% |

Typical GI tumour purity is roughly 30–70%. Only *MYC* crosses inside that band.
Every other gene requires a purity higher than GI tumours plausibly reach.

**This sweep should be a figure in the paper.** It converts "the threshold is
sensitive to the estimand" from a reviewer's objection into a reported result, and
it is the analysis that makes the threshold choice defensible.

## Two specific checks requested

**BCL2's 51% "other" in the HPA figure was a collapse artifact.** In tumour data
with proper labels there is no "other" mass: BCL2 is **65% lymphoid** in iCCA
pseudobulk. The critique's suspicion was correct and the cause was the normal-tissue
cell-type vocabulary (pituitary stem cells, melanocytes, adipocytes) not mapping onto
tumour compartments.

**SOCS3 is confirmed non-epithelial in tumours, and the biology is more interesting
than "stromal."** In iCCA it is 48% lymphoid, 20% myeloid, 15% fibroblast, 15%
malignant. As the most canonical direct STAT3 negative-feedback target in the set,
its being immune-dominated rather than tumour-dominated is a substantive finding.
Note this is *immune*, not *fibroblast* — the HPA figure's "62% stromal" was wrong
about the compartment as well as the magnitude.

## A caveat neither of us named

The iCCA and HCC columns disagree sharply (*MYC* 0.67 vs 0.11; *IL6* 0.12 vs 0.01).
Some of that is real biology, but some is **dissociation bias**: malignant cells are
34.5% of captured iCCA cells and 25.9% of captured HCC cells, and neither reflects
true tissue composition, because dissociation protocols recover immune cells far more
efficiently than epithelium.

This means abundance-weighted pseudobulk from a dissociated atlas is *also* not the
bulk estimand — it is a different biased estimate. The fix is the purity-sweep
approach above: use single-cell data for compartment-specific **intensities** (which
are composition-independent) and supply the composition from bulk purity estimates
(ESTIMATE, CPE) rather than from cell counts. That is what the sweep does, and it is
what should be preregistered.

Per-patient variability is also large — *MYC*'s epithelial fraction ranges 0.12 to
0.98 across 10 iCCA patients (median 0.62). Any threshold rule must be applied to a
cohort-level summary with per-patient variation reported, not to a single pooled
number.

## Revised recommendation

1. Use **compartment intensity × bulk-purity composition**, not raw pseudobulk shares
   and not HPA cluster shares. Report the purity sweep as a figure.
2. Widen to 30–50 canonical STAT3 targets — for **scope** (characterising the class
   of score), which was the critique's correct reason, not to rescue degeneracy.
3. Set the threshold **after** running the widened panel through the sweep, and
   preregister the sweep procedure rather than a fixed cutoff.
4. Report per-patient variation, and report iCCA and HCC separately — the divergence
   is itself a finding about tissue context.
5. GSE125449 pairs with FU-iCCA bulk (NODE OEP001105) as suggested; that pairing is
   the right validation axis for the biliary arm.
