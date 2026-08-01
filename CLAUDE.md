# Project brief — compartment decomposition of STAT3 activity scores in GI cancers

Handoff document. Read fully before writing code.

Companion files that belong in this repo:
- `panel_definition.md` — the prespecification. **Binding.** Do not edit after lock.
- `feasibility_assessment.md` — prior-art review and pilot compartment results. **Authoritative for all pilot numbers.** Do not restate figures from memory; read them from this file.
- `01_download.R`, `02_panel.R` — drafted, untested.

---

## 1. Who and what came before

Sole/first author is a high-school researcher (Cinco Ranch HS, Katy TX) working
with collaborators at Baylor College of Medicine and MD Anderson. Two prior
works, both already submitted or presented — **neither is to be revised**:

**(a) CCRC-D-26-00363**, under review at *Clinical Colorectal Cancer*. A
complete-enumeration meta-analysis of all four registrational Phase 3
napabucasin trials in GI cancer (CO.23, BRIGHTER, CanStem303C, CanStem111P;
n=3,383). Pooled general/ITT OS-HR 1.028, PFS-HR 1.027, both at the null. Key
result for this project: the CO.23 pSTAT3-positive OS-HR of 0.41 (n=55) did not
replicate in CanStem303C (0.969); the two differ significantly (z=2.77,
P=0.006). Control-arm contrasts in BRIGHTER (1.32) and CanStem303C (1.518)
showed pSTAT3 is **prognostic, not predictive**.

**(b) ESMO Asia 2026 abstract / Cholangiocarcinoma Foundation poster.** A
6-gene transcriptomic pSTAT3 score (`SOCS3, BCL2, MYC, MMP9, HGF, IL6`) in
biliary tract cancer. High score → shorter PFS and OS, Treg/Th2-skewed
immunosuppressive TME, stromal enrichment, reduced *SLC29A1*.

Known weaknesses in (b) that this project exists partly to address: the six
genes were selected in TCGA-CHOL for tracking RPPA STAT3_pY705 and then
"validated" against the same data (in-sample); genes were averaged as
`mean(log2(expr+1))` without z-scoring; the cutoff derivation is unstated.

---

## 2. The thesis

A biomarker used to select patients in a Phase 3 programme is largely a readout
of **non-tumour compartments** — which is a concrete, testable explanation for
why the CO.23 pSTAT3 signal did not replicate.

This is deliberately *not* framed as "STAT3 scores are stromal." That general
claim is already established (see §4). The contribution is the clinically
anchored application plus the link to a documented biomarker non-replication.

The two papers interlock: prognostic-but-not-predictive is exactly what a
stroma-dominant score looks like.

---

## 3. Design — read this before proposing any analysis

**The primary endpoint does not depend on the compartment decomposition.**

- **Primary:** purity-adjusted Cox model on the **intact** score, per cohort,
  meta-analysed. Requires no compartment fractions. Cannot degenerate. The
  estimand is the *attenuation* in log-HR on adding purity and stromal score,
  not a binary significant/not verdict.
- **Secondary / explanatory:** compartment decomposition, gated by the decision
  rule in `panel_definition.md` §5.

This inversion is deliberate. An earlier design made the primary analysis
depend on splitting the score into epithelial and stromal subscores; pilot work
showed the epithelial arm may not exist. Do not revert to that design.

Practical consequence: cohorts with no usable single-cell atlas (ESCA, likely
CHOL) still contribute fully to the primary result.

### The three-branch contingency

Pilot data suggest very few canonical STAT3 targets are epithelial-dominant.
`panel_definition.md` §5 commits in advance to what happens at each possible
count. **All three branches are publishable. Do not argue toward a branch.**

---

## 4. Prior art that constrains the framing

These must be cited and explicitly responded to. Do not present any of them as
novel findings of this work.

| Work | What it established | Required response |
|---|---|---|
| Isella et al. 2015 (Nat Genet) | Stromal contribution to the CRC transcriptome | Cite; our claim is not this |
| Calon et al. 2015 (Nat Genet) | Stromal gene expression defines poor-prognosis CRC | Same |
| Guinney et al. 2015 (Nat Med) | CMS; CMS4 is mesenchymal/stromal | Show the STAT3 score is not just CMS4 — orthogonality analysis required |
| Venet et al. 2011 (PLoS Comput Biol) | Most random signatures associate with outcome | This is a **required control, not a contribution** |
| Ovarian analogue (CEBP 2020) | Same study design, different tumour | Cite as precedent for the design |
| Aran et al. 2015 (Nat Commun ncomms9971) | Purity confounds bulk expression; CPE | Method source |
| Liu et al. 2018 (Cell) | TCGA-CDR curated endpoints | Use their per-cancer endpoint recommendations |

Novelty gap identified by literature sweep: no work asks whether a STAT3
activity score is a purity artefact, and no prespecified epithelial/stromal
decomposition of a named clinical score exists.

---

## 5. Pilot results already in hand

See `feasibility_assessment.md` for exact numbers. Summary of what is settled:

- HPA single-cell consensus was used as a first pass, then **correctly
  discarded as a primary source** — its metric sums per-cell-type expression
  and ignores cell abundance, so it inflates rare compartments.
- The correct estimand is **abundance-weighted contribution to pseudobulk**,
  computed in tumour single-cell data.
- Recomputing on GSE125449 with abundance weighting did **not** rescue the
  epithelial arm. Only *MYC* exceeded 50% epithelial, and its per-patient
  fraction varied enormously across ~10 iCCA patients.
- A one-gene "epithelial subscore" that is *MYC* is a MYC score, not a STAT3
  score. Hence the panel widening and the contingency rule.

Open caveat: GSE125449 is 19 tumours. Compartment estimates from it are
fragile. The CRC atlas (62 patients) will be far more stable — report both.

---

## 6. Data inventory

### Tier 1 — bulk discovery (required)
| Item | Access |
|---|---|
| TCGA COAD, READ, STAD, ESCA, PAAD, LIHC, CHOL (STAR-Counts) | `TCGAbiolinks`; https://portal.gdc.cancer.gov/ |
| TCGA-CDR survival (Liu 2018) | https://gdc.cancer.gov/about-data/publications/pancanatlas |
| Aran CPE purity | https://www.nature.com/articles/ncomms9971 → Suppl. Data 1 |

ESTIMATE is computed, not downloaded.

### Tier 2 — compartment atlases
| Tissue | Accession | Link |
|---|---|---|
| Liver / biliary | GSE125449 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449 |
| Colorectal | GSE178341 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341 |
| Colorectal (pre-annotated mirror) | SCP1162 | https://singlecell.broadinstitute.org/single_cell/study/SCP1162/human-colon-cancer-atlas-c295 |
| Pancreatic | GSE155698 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155698 |
| Pancreatic (annotated, 10 cell types) | Peng et al. via Zenodo | https://doi.org/10.5281/zenodo.3969339 |
| Gastric (weak — premalignant/early skew) | GSE134520 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134520 |

**Unresolved:** gastric needs a better atlas. Two candidates, accessions to be
pulled from their data-availability statements:
- Kumar et al., Cancer Discov 2022 — https://aacrjournals.org/cancerdiscovery/article/12/3/670/681898/Single-Cell-Atlas-of-Lineage-States-Tumor
- Integrated GI atlas, Sci Data 2026 — https://www.nature.com/articles/s41597-026-07108-3 (574,532 stomach cells / 479,629 colorectal cells, harmonised annotations — **check this first**, harmonised vocabulary across two tissues would remove a real comparability problem)

No usable ESCA atlas is expected. Handle per §3.

### Tier 3 — external validation
| Item | Access |
|---|---|
| GSE39582 (CRC, ~560 w/ survival) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582 |
| GSE66229 (ACRG gastric) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66229 |
| ICGC pancreatic | https://dcc.icgc.org/ (portal has migrated before; verify) |
| FU-iCCA n=255 bulk + phosphoproteome | https://www.biosino.org/node/project/detail/OEP001105 |
| Dong et al. Table S1 | https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00659-0 |

**Manual only** — NODE requires account registration, and the Dong supplement
is a paper download. Do not attempt to automate account creation.

Note on the FU-iCCA phosphoproteome: it is measured in **bulk tissue**, so it
carries the same compartment ambiguity as the mRNA score. It validates that the
score tracks STAT3 phosphorylation, **not** which cells the phosphorylation is
in. Do not overclaim it.

### Tier 4 — annotation
- MSigDB `HALLMARK_IL6_JAK_STAT3_SIGNALING` — https://www.gsea-msigdb.org/gsea/msigdb
- ChEA3 / TRRUST v2 for direct-target evidence (see `panel_definition.md` §2B)
- `CMScaller` or equivalent for the CMS orthogonality analysis

**Avoid:** CRC raw data at https://www.ega-archive.org/studies/phs002407 is
controlled-access and needs a DAC application. Use GEO or Broad instead.

---

## 7. Script plan

| Script | Purpose |
|---|---|
| `01_download.R` | Retrieval, caching to `.rds`. Drafted, untested. |
| `02_panel.R` | Implements `panel_definition.md`. Criterion B left unfilled deliberately. |
| `03_compartments.R` | Abundance-weighted pseudobulk fractions per atlas; purity sweep 30–70%; bootstrap over patients |
| `04_clinical.R` | TCGA-CDR merge; per-cancer endpoint selection; primary-tumour filter; dedup |
| `05_purity.R` | ESTIMATE + Aran CPE; cross-method agreement as sensitivity |
| `06_score.R` | Z-score within cohort, average. Functions, not scripts. |
| `07_survival.R` | Cox per cohort, unadjusted then purity-adjusted. `cox.zph`. |
| `08_pool.R` | `metafor` meta-analysis of per-cohort log-HRs |
| `09_null.R` | ~1,000 random gene sets matched on mean expression and stromal correlation |
| `10_validation.R` | GSE39582, ACRG, FU-iCCA — reuses 06/07 as functions |
| `11_figures.R` | Compartment bars, forest, null histogram, KM |

### Compartment computation — the one thing not to get wrong

```r
# counts: genes x cells, RAW. meta$compartment: factor.
pseudo <- t(rowsum(t(as.matrix(counts)), group = meta$compartment))
frac   <- sweep(pseudo, 1, rowSums(pseudo), "/")
```

Summing **raw** counts by compartment is abundance-weighted by construction.
Normalising per cell type first reproduces the HPA error. Do not normalise
first.

---

## 8. Conventions and gotchas

- `options(timeout = 3600)` before any GDC/GEO download — the 60s default fails
  as apparent file corruption.
- Z-score **within cohort**, never across. Batch.
- Scale score per SD so HRs are comparable across cohorts.
- Harmonise stage coding; it differs by TCGA project.
- Set a seed in `09_null.R`.
- Use `renv`; commit `renv.lock`. Version pinning matters more than usual given
  preregistration.
- **Write 06 and 07 as pure functions taking a matrix**, not scripts on globals.
  They run four-plus times on different cohorts; every copy-paste variant is a
  silent inconsistency.
- GSE125449's supplementary file layout has changed across GEO revisions —
  inspect before writing parsers.

---

## 9. Hard constraints

1. **Do not fit any survival model before the panel is locked** and
   `panel_locked.csv` is committed. Tooling speed makes accidental
   out-of-order execution the main threat to the paper's central claim.
2. **Compartment fraction is never an inclusion or exclusion criterion** for the
   panel. That is the outcome, not a filter.
3. Preregister on OSF after panel lock and sweep specification, **before** any
   survival fit. Disclose that pilot compartment fractions for the original six
   were inspected beforehand.
4. The prior work in §1 is submitted. This project may *reinterpret* it but must
   not be described as correcting it. The discussion must explicitly address
   that the ESMO abstract reports the score as a validated prognostic biomarker
   while this paper argues it is substantially stromal — that self-correction
   is legitimate and must be stated openly, not buried.
5. Reviewer findings must be queried explicitly. `host.findings()` is called at
   the end of every work item and ALL findings reported regardless of severity.
   Warn-level findings are never injected into the working conversation — during
   Part A, ten warns sat unaddressed while work continued on top of them,
   including two bearing on a numerical claim underpinning Amendment 13. Never
   rely on injection. Report the full list, with disposition for each.

---

## 10. Open decisions

- Gastric atlas selection (§6, Tier 2).
- Whether ESCA stays in the primary and out of the decomposition, or is dropped.
- Criterion B evidence sources and versions (`panel_definition.md` §2B).
- Whether Baylor / MD Anderson collaborators join as co-authors — relevant
  because this partly reframes their shared prior work.
- Target venue. Working assumption: specialty journal (*British Journal of
  Cancer*, *npj Precision Oncology*), and AACR Annual Meeting for the abstract.
  Lead with the napabucasin link — it is what lifts this above a routine
  signature paper.

---

## 11. Immediate next step

Run `02_panel.R` after filling criterion B. Report the gene count and the
route breakdown (A only / B only / both), and whether any of the original six
failed. That number determines whether the decomposition is a subscore analysis
or descriptive only, and therefore what `03_compartments.R` needs to do.