# STAT3 activity, tumour purity and outcome in GI cancers — study summary

Standalone record of Parts A and B, written for the paper's methods and
discussion. Every figure here is read from a committed output file; the
provenance column in each table names it.

**Preregistration:** `panel_definition.md` (panel, exclusion criteria, Amendments
1–15) and `analysis_plan.md` v1.5 (sections A.a–A.g, B.h–B.o). The panel was
locked at commit `ac9c5e0` before any compartment fraction was computed.

**Scripts:** `01_download.R`, `02_panel.R`, `03_compartments.R`,
`04_lock_gene_list.R` (Part A); `05_clinical.R`, `06_purity.R`, `07_score.R`,
`08_survival.R`, `09_null.R` (Part B). Script 10 (external validation) is not run.

---

# PART A — panel, atlases, compartment estimand

## A.1 Panel construction

A gene enters the panel if it meets **A**, **B**, or both.

- **A — Pathway membership:** member of MSigDB `HALLMARK_IL6_JAK_STAT3_SIGNALING`.
- **B — Direct transcriptional target:** at least two independent lines of
  evidence at default thresholds, from ChIP-seq-derived TF-target databases,
  curated TF-target databases (TRRUST v2), or a promoter-level STAT3 binding site
  in primary literature. Two sources querying the same underlying ChIP-seq
  experiment count as one.

**Result: 152 genes** (`data/panel/panel_locked.csv`).

| Route | n |
|---|---|
| A only (pathway membership) | 77 |
| B only (direct target evidence) | 65 |
| Both | 10 |
| **Total** | **152** |

### The origin six

The six genes that motivated the study were assessed against the same criteria.
**Three qualified, three did not** (`data/panel/origin_six_evidence.csv`):

| Gene | Route | Human ChIP-seq sources | In panel | Reason if excluded |
|---|---|---|---|---|
| SOCS3 | A + B | 2 | **yes** | — |
| MYC | B | 2 | **yes** | — |
| IL6 | A + B | 2 | **yes** | — |
| BCL2 | — | 0 | no | no ChIP-seq evidence in any group; fails criterion B |
| MMP9 | — | 0 | no | mouse ChIP-seq only; excluded by Amendment 2 (human-only) |
| HGF | — | 0 | no | mouse ChIP-seq only; excluded by Amendment 2 (human-only) |

The three non-qualifying genes are **reported throughout as a labelled
non-qualifying subset** (prespecification §4, Amendment 2). They enter no k
variant, no dominance matrix, and no evaluability count, but their compartment
fractions are computed and reported so readers can see whether excluding them was
consequential. This is why several outputs carry 155 rows (152 + 3) while every
inferential quantity is computed over 152.

## A.2 Atlas selection and attrition — five to three

The registration named five compartment atlases. **Two were removed after direct
inspection of the deposits**, both for the same reason.

| Atlas | Tissue | Status | What direct inspection found |
|---|---|---|---|
| GSE125449 | liver/biliary | **retained** | Two deposited sets with different gene universes (20,124 vs 19,572 rows) — see Amendment 12 |
| GSE178341 | colorectal | **retained** | 43,113 × 370,115; too large to materialise in 16 GB, read by a chunked subset reader |
| Peng (PRJCA001063) | pancreatic | **retained** | `raw/X` contains no count layer; it is log1p-CP10K — see Amendment 13 |
| GSE155698 | pancreatic | **removed (Amendment 9)** | `GSE155698_RAW.tar` holds 41 per-sample CellRanger archives (barcodes, features, matrix) with **no cell-type annotation table, no metadata file, no annotation column in any sample**. The GEO series has no other supplementary file. The annotations exist in Steele et al. but were not deposited. |
| GSE183904 | gastric | **removed (Amendment 10)** | 40 flat `.csv.gz` members, one per GSM, each a gene × cell count matrix (26,572 rows, barcode headers). **No annotation row, no metadata table, no cluster file, in any sample.** Annotations exist in Kumar et al. but were not deposited. |

Both amendments record why the obvious workarounds were rejected. Deriving
cell-type labels de novo would have made those tissues the only ones whose
compartment labels are this project's own construction, on data of different
platform and depth, reintroducing the cross-atlas asymmetry Amendment 4 was
written to eliminate. Substituting a different atlas would have been selecting a
replacement dataset after registration.

**Consequence:** no gastric or oesophageal atlas remains. TCGA-STAD and TCGA-ESCA
are scored and modelled in Part B with no tissue-matched compartment atlas, which
Amendment 10 records as a stated limitation rather than a resolved problem.

**This is a candidate finding for the discussion** (`NOTES_FOR_REVIEW.md` §13):
three of the five originally selected atlases carried deposit defects invisible
from their publications — two with no deposited annotation at all, one with no
count layer. Amendment 4's claim that lineage-level labels were "available and
comparably defined in all five" atlases was derived from the publications'
methods sections rather than from the deposits, and was wrong for two of five.

## A.3 The compartment estimand, and how it changed

**Amendment 4** changed the estimand from *malignant-epithelial fraction* to
**epithelial fraction on tumour-channel samples only**. Malignancy was inferred by
three different methods across the atlases (`NOTES_FOR_REVIEW.md` §7), so a
malignant-cell estimand would not have been comparably defined; the epithelial
lineage label is.

**Amendment 11** changed the Peng atlas's compartment mapping from `celltype0` to
`celltype1`. Verified before adoption, not assumed: `celltype1` is a strict
refinement of `celltype0`, and **epithelial membership is bitwise identical across
levels** (24,156 cells either way). The amendment's no-bias claim therefore holds
by construction — the change affects only how non-epithelial mass is subdivided,
which alters neither the numerator nor the denominator of f(π).

## A.4 Part A results

Dominance rule (Amendments 5, 10): a gene is epithelial-dominant if it is dominant
**in at least two of the three atlases, of any GI tissue**.

| Quantity | Estimate | 95% CI | Provenance |
|---|---|---|---|
| **k** (dominant, full 30–70% band) | **46** | 38–57 | `output/k_estimates.csv` |
| k_50 (dominant at π = 0.50 only) | 96 | 86–103 | same |
| k_all3 (dominant in all three atlases) | 24 | 15–31 | same |
| k_evalall (dominant, evaluable in all three) | 43 | 35–55 | same |

Intervals are patient-level bootstrap percentile CIs, B = 2,000, seed 20260731.

**Amendment 3 branch: `epithelial subscore BUILT; subscore survival models run as
SECONDARY`.** k = 46 of 152 is 30.3%, above the 20% threshold. The Amendment 6
disclosures (`output/amendment6_disclosures.txt`) record that the **entire CI
selects the same branch** (branch at 38 = branch at 57 = branch at 46), so no
band-straddle disclosure is required, and k_50 selects the same branch as k.

**Per-atlas dominance** (`output/per_tissue_dominance.csv`): GSE178341 colorectal
67, GSE125449 liver/biliary 43, Peng pancreatic 38.

**Evaluability** (`output/evaluability_distribution.csv`), over the 152 panel
genes: evaluable in 3 atlases 140, in 2 atlases 9, in 1 atlas 3, in 0 atlases 0.

**The origin six in the compartment data** (`output/origin_six_compartment.csv`):
of the three qualifying genes, MYC is epithelial-dominant in 1 of 3 atlases, SOCS3
in 0 of 3, IL6 in 0 of 3. The three non-qualifying genes are reported with their
fractions and are dominant in none.

## A.5 The final gene list

k is computed over the **locked 152-gene panel** (Amendment 14), not over the
scoring list. The scoring list applies prespecification exclusion rules §3.1–§3.3
to the 152.

- **Primary scoring list: 140 genes** (`data/panel/final_gene_list_140.csv`,
  `list_id = "primary_140"`, md5 of sorted symbols
  `eb167e8c7a33b4202bd609a17defa629`)
- **Sensitivity list: 143 genes** (`final_gene_list_143.csv`,
  `list_id = "sensitivity_143"`, md5 `8a54835eedcfe3b23dcb56ce74805a87`)

The 143 list is the narrow reading of rule 3.3 (pseudoautosomal genes retained);
Amendment 15 settled the broad reading, making 140 primary. The 140 is a strict
subset of the 143, asserted at every read.

### All 12 exclusions

| Gene | Rule | Max detection | Atlases evaluable |
|---|---|---|---|
| DNTT | 2 — undetectable (<1% in every compartment) | 0.027% | 1 |
| GFAP | 2 | 0.502% | 3 |
| LEP | 2 | 0.303% | 1 |
| OPRM1 | 2 | 0.872% | 2 |
| PAX3 | 2 | 0.136% | 2 |
| CRLF2 | **2 + 3.3** (undetectable *and* pseudoautosomal) | 0.139% | 1 |
| IL13RA1 | 3.3 — chrX | 44.92% | 3 |
| IL2RG | 3.3 — chrX | 56.72% | 3 |
| TIMP1 | 3.3 — chrX | 98.38% | 3 |
| CSF2RA | 3.3 — pseudoautosomal (chrX/chrY) | 44.10% | 3 |
| IL3RA | 3.3 — pseudoautosomal | 61.10% | 3 |
| IL9R | 3.3 — pseudoautosomal | 4.17% | 2 |

152 − 12 = 140; the sets are disjoint and partition the panel (CRLF2 is caught by
both rules, which is why the union is 12 and not 13).

**No origin-six gene is excluded** under either list.

Rule §3.1 (absent from any TCGA cohort annotation) excluded nothing: all 152
genes are present in all seven cohorts. That path can only halt, never exclude —
an earlier consistency check fires first if a gene is missing.

**Two findings recorded against this list.** First, half the rule-2 exclusions
rest on a single atlas: DNTT, LEP and CRLF2 were judged on one atlas each while
138 of the 140 retained genes were held to all three (`NOTES_FOR_REVIEW.md` §19).
Second, the assembly is GRCh38, established by matching MYC's coordinates
(`chr8:127,735,434–127,742,951`) rather than read from a version field — the
deposit objects record "Data Release 45.0" and **no GENCODE version anywhere**
(§25). An earlier claim of "GENCODE v36" was inferred from a `_PAR_Y` naming
convention and has been withdrawn.

---

# PART B — clinical, purity, score, survival, null benchmark

Seven TCGA cohorts: COAD, READ, STAD, ESCA, PAAD, LIHC, CHOL. Six enter the
meta-analysis; **CHOL is descriptive-only** (Amendment 8: a cohort whose TCGA-CDR
Table 3 explanation states the sample size is too small for all endpoints is
excluded as a weighted stratum regardless of EPV).

## B.1 Clinical assembly (05)

Endpoints follow Amendment 7, a deterministic rule applied identically to all
seven: primary = OS wherever CDR Table 3 marks OS usable without caution. **READ
is the only cohort whose OS mark is cautioned, so READ is PFI and the other six
are OS.** Event counts play no part in the designation. Administrative censoring:
10 years OS, 5 years PFI, applied identically.

| Cohort | Endpoint | n | Events (primary) | Events (sensitivity) | Stage missing | Redacted | use_sex | Meta |
|---|---|---|---|---|---|---|---|---|
| COAD | OS | 455 | 102 | 117 | 11 | 0 | yes | yes |
| READ | **PFI** | 165 | 36 | 25 | 9 | 0 | yes | yes |
| STAD | OS | 406 | 157 | 130 | 24 | 0 | yes | yes |
| ESCA | OS | 184 | 77 | 87 | 23 | 0 | yes | yes |
| PAAD | OS | 178 | 93 | 104 | 3 | 0 | yes | yes |
| LIHC | OS | 370 | 130 | 178 | 24 | 1 | yes | yes |
| CHOL | OS | 35 | 18 | 20 | 0 | 0 | yes | **no** |

**n across the six meta-eligible cohorts = 1,758** (1,755 after complete-case
attrition in 08). Source: `output/clinical_cohort_summary.csv`.

**Filter flow** (`output/clinical_filter_flow.csv`), reconciling exactly per
cohort: 1,998 samples in → 171 not primary tumour → 23 duplicate aliquots → 3 no
CDR match → 8 no usable endpoint → **1,793 final**.

Aliquot tie-break (registered decision, 2026-08-02): earliest vial letter, then
lowest portion, then lowest plate — deliberately **not** a full-barcode sort,
which leads with plate and centre codes and so orders on processing batch and
sequencing site. Only COAD had multi-aliquot patients: 3 with two aliquots, 10
with three.

Redaction: retained in the primary, counted per cohort, with a prespecified
exclusion sensitivity. One redacted sample in the entire set (LIHC, 0.3%).

### Realised EPV and dispositions

EPV = realised events ÷ registered parameter count, **per model** (Amendment 8).
Parameter counts are fixed by registration at M1=1, M2=5, M3=6, M4=7, with
`stage_group` contributing two parameters for its three levels.

| Cohort | M1 | M2 | M3 | M4 | Disposition | Provisional | Differs |
|---|---|---|---|---|---|---|---|
| COAD | 102.0 | 20.4 | 17.0 | 14.6 | fitted | fitted | no |
| READ | 36.0 | **7.2** | **6.0** | **5.1** | low-EPV | low_EPV | no |
| STAD | 157.0 | 31.4 | 26.2 | 22.4 | fitted | fitted | no |
| ESCA | 77.0 | 15.4 | 12.8 | 11.0 | fitted | fitted | no |
| PAAD | 93.0 | 18.6 | 15.5 | 13.3 | fitted | fitted | no |
| LIHC | 130.0 | 26.0 | 21.7 | 18.6 | fitted | fitted | no |
| CHOL | 18.0 | **3.6** | **3.0** | **2.6** | cohort excluded | cohort_excluded | no |

**The realised disposition matches the provisional table in all seven cohorts.**
READ is flagged low-EPV on M2–M4 and **included** in the primary pool (Amendment 8
flags, it does not exclude), with a drop-all-low-EPV pool reported separately.
CHOL's M2–M4 are **not fitted** (EPV < 5).

## B.2 Purity (06)

`data/manual/aran_purity.xlsx` was absent at the start of Part B. The plan (§3.4)
anticipated this and prespecifies a switch rule, but letting a *missing download*
flow through that rule would have inverted B.j's registered CPE-primary ordering
for all seven cohorts on an acquisition gap rather than on the data. The file was
obtained (Aran et al. 2015, *Nat Commun* 6:8971, Supplementary Data 1; md5
`c459e6a965789b96860fc77bd346c681`, 9,364 rows, 21 cancer types) and the rule
applied to genuine coverage. **06 now halts if it is absent.**

| Cohort | CPE coverage | Primary source | Purity mean (SD) | Agreement n | Pearson | Spearman |
|---|---|---|---|---|---|---|
| COAD | 99.8% | **CPE** | 0.770 (0.103) | 455 | 0.698 | 0.654 |
| READ | 100% | **CPE** | 0.780 (0.081) | 165 | 0.659 | 0.590 |
| LIHC | 100% | **CPE** | 0.777 (0.127) | 370 | 0.736 | 0.693 |
| STAD | 0% | ESTIMATE | 0.670 (0.157) | 0 | — | — |
| ESCA | 0% | ESTIMATE | 0.789 (0.130) | 0 | — | — |
| PAAD | 0% | ESTIMATE | 0.605 (0.162) | 0 | — | — |
| CHOL | 0% | ESTIMATE | 0.791 (0.125) | 0 | — | — |

The four zeroes are **not a join failure**: STAD, ESCA, PAAD and CHOL are absent
from the Aran 2015 freeze entirely, verified by listing the table's own 21 cancer
types. The switch rule (CPE primary at ≥80% coverage) therefore applies exactly as
registered.

**Consequence for the reader** (`NOTES_FOR_REVIEW.md` §22): the purity covariate
is not one quantity across cohorts. Three use a consensus of four orthogonal
methods; four use an ESTIMATE-derived value whose conversion
`cos(0.6049872018 + 0.0001467884 · S)` was calibrated on **Affymetrix arrays** and
is applied here to Illumina RNA-seq TPM. `purity_calibration` in
`output/purity_summary.csv` records which is which per cohort. The registered
cross-method agreement is computable only in the three CPE cohorts.

All 1,793 purity values are bounded [0,1]; ESTIMATE scored an identical 9,911-gene
set in every cohort.

## B.3 Score construction (07)

Registered pipeline (B.h): primary tumour, one aliquot per patient, log2(TPM+1),
z-score **each gene within cohort**, mean across the final gene list, scale the
score to unit SD within cohort. No dichotomisation, no tertiles, no data-derived
cutpoint.

**Within-cohort SD is exactly 1 in all seven cohorts**, mean 0 (per-gene z-scoring
makes the column mean exactly zero before scaling). The 140- and 143-gene scores
correlate at r = 0.9992–0.9998. **No zero-variance genes** in any cohort, so all
140/143 genes survive the global exclusion rule (registered decision D: excluded
globally, never per cohort).

### The stromal collinearity result

| Cohort | r(score, panel subscore) | Exceeds 0.9 |
|---|---|---|
| COAD | 0.982 | yes |
| READ | 0.981 | yes |
| STAD | 0.965 | yes |
| ESCA | 0.962 | yes |
| PAAD | 0.978 | yes |
| LIHC | 0.986 | yes |
| CHOL | 0.957 | yes |

**The registered ESTIMATE stromal fallback triggered in all seven cohorts.** B.h
prespecifies this: "If Part A returns the descriptive branch, the stromal subscore
is defined on all non-dominant panel genes, which is nearly the whole panel — in
that case model 4 is collinear with model 1 by construction, and the plan
substitutes a purity-orthogonal stromal index (ESTIMATE stromal score) instead.
Stated now so the choice is not made after seeing k." `stromal_score` in M4 is
therefore the standardised ESTIMATE StromalScore throughout; the panel-derived
subscore is retained as `stromal_score_subscore`.

The covariate actually used correlates with the main score far more weakly than
the one it replaced (`output/stromal_correlations.csv`):

| Cohort | r(score, **ESTIMATE stromal**) — used | r(score, panel subscore) — replaced | r(score, purity) |
|---|---|---|---|
| COAD | 0.666 | 0.982 | −0.484 |
| READ | 0.622 | 0.981 | −0.456 |
| STAD | 0.413 | 0.965 | −0.427 |
| ESCA | 0.434 | 0.962 | −0.398 |
| PAAD | 0.634 | 0.978 | −0.603 |
| LIHC | 0.399 | 0.986 | −0.433 |
| CHOL | 0.600 | 0.957 | −0.439 |

## B.4 Survival — the primary result (08)

Four nested Cox models per cohort on an **identical complete-case set** (complete
cases across M4's covariates, so attenuation is not confounded by changing n):

```
M1: Surv(time, event) ~ score                                        1 parameter
M2: + age + sex + stage_group                                        5
M3: + purity                                                         6
M4: + stromal_score                                                  7
```

### Per-cohort score log-HR

| Cohort | n | Events | Model | β | SE | HR (95% CI) | p |
|---|---|---|---|---|---|---|---|
| COAD | 455 | 102 | M1 | 0.0052 | 0.0943 | 1.005 (0.836–1.209) | 0.956 |
| | | | M2 | 0.0046 | 0.0981 | 1.005 (0.829–1.218) | 0.963 |
| | | | M3 | −0.0197 | 0.1188 | 0.981 (0.777–1.238) | 0.868 |
| | | | M4 | −0.0713 | 0.1378 | 0.931 (0.711–1.220) | 0.605 |
| READ | 165 | 36 | M1 | −0.0351 | 0.1669 | 0.966 (0.696–1.339) | 0.833 |
| | | | M2 | −0.0457 | 0.1577 | 0.955 (0.701–1.301) | 0.772 |
| | | | M3 | 0.0165 | 0.1877 | 1.017 (0.704–1.469) | 0.930 |
| | | | M4 | −0.0221 | 0.2218 | 0.978 (0.633–1.511) | 0.921 |
| STAD | 403 | 157 | M1 | 0.1131 | 0.0835 | 1.120 (0.951–1.319) | 0.175 |
| | | | M2 | 0.0755 | 0.0855 | 1.079 (0.912–1.275) | 0.377 |
| | | | M3 | 0.0122 | 0.0920 | 1.012 (0.845–1.212) | 0.894 |
| | | | M4 | 0.0084 | 0.0920 | 1.008 (0.842–1.208) | 0.927 |
| ESCA | 184 | 77 | M1 | 0.0542 | 0.1157 | 1.056 (0.842–1.324) | 0.639 |
| | | | M2 | 0.1465 | 0.1182 | 1.158 (0.918–1.460) | 0.215 |
| | | | M3 | 0.2116 | 0.1370 | 1.236 (0.945–1.616) | 0.122 |
| | | | M4 | 0.2078 | 0.1392 | 1.231 (0.937–1.617) | 0.136 |
| PAAD | 178 | 93 | M1 | 0.2769 | 0.1179 | 1.319 (1.047–1.662) | 0.019 |
| | | | M2 | 0.2754 | 0.1211 | 1.317 (1.039–1.670) | 0.023 |
| | | | M3 | 0.4578 | 0.1692 | 1.581 (1.135–2.202) | 0.007 |
| | | | M4 | 0.4111 | 0.1733 | 1.508 (1.074–2.119) | 0.018 |
| LIHC | 370 | 130 | M1 | 0.2331 | 0.0885 | 1.263 (1.062–1.502) | 0.008 |
| | | | M2 | 0.2339 | 0.0904 | 1.264 (1.058–1.509) | 0.010 |
| | | | M3 | 0.3215 | 0.1048 | 1.379 (1.123–1.694) | 0.002 |
| | | | M4 | 0.3626 | 0.1068 | 1.437 (1.166–1.772) | 0.0007 |

**CHOL, descriptively** (n=35, 18 events, never pooled): M1 β = −0.1606,
SE 0.2168, HR 0.852 (0.557–1.303), p = 0.459. M2–M4 not fitted (EPV rule).

### Attenuation per cohort

Primary estimand `attenuation_total = β(M2) − β(M4)`. Interval by paired
nonparametric bootstrap over patients within cohort, **B = 2,000, seed 20260731**,
refitting all four models in each resample so β₂ and β₄ carry their correlation.
All 2,000 resamples succeeded in every cohort.

| Cohort | β(M2) | β(M4) | attenuation_total | 95% CI (percentile) | SE |
|---|---|---|---|---|---|
| COAD | 0.0046 | −0.0713 | **0.0759** | −0.1519 to 0.2967 | 0.1122 |
| READ | −0.0457 | −0.0221 | **−0.0236** | −0.4363 to 0.3453 | 0.1934 |
| STAD | 0.0755 | 0.0084 | **0.0671** | −0.0089 to 0.1653 | 0.0434 |
| ESCA | 0.1465 | 0.2078 | **−0.0614** | −0.2327 to 0.0631 | 0.0752 |
| PAAD | 0.2754 | 0.4111 | **−0.1357** | −0.5065 to 0.1065 | 0.1552 |
| LIHC | 0.2339 | 0.3626 | **−0.1286** | −0.2775 to −0.0367 | 0.0619 |

Components — `attenuation_purity` (β₂−β₃): COAD 0.0243, READ −0.0622, STAD 0.0633,
ESCA −0.0652, PAAD −0.1824, LIHC −0.0876. `attenuation_stroma` (β₃−β₄): COAD
0.0516, READ 0.0385, STAD 0.0038, ESCA 0.0038, PAAD 0.0468, LIHC −0.0410.

`prop_attenuated` (secondary): READ 0.517, STAD 0.889, ESCA −0.419, PAAD −0.493,
LIHC −0.550. **NA for COAD and CHOL** — β₂ and β₄ have opposite signs, which B.k
states makes it undefined.

### Meta-analysis (HKSJ primary, k = 6, n = 1,755)

| Model | est | 95% CI | p | HR (95% CI) | τ² | I² | Q | Q p |
|---|---|---|---|---|---|---|---|---|
| M1 | 0.1212 | 0.0002 to 0.2421 | 0.0497 | 1.129 (1.000–1.274) | 0.0022 | 16.8% | 6.08 | 0.299 |
| M2 | 0.1229 | 0.0001 to 0.2457 | 0.0499 | 1.131 (1.000–1.279) | 0.0022 | 16.1% | 6.04 | 0.303 |
| M3 | 0.1572 | −0.0414 to 0.3558 | 0.0975 | 1.170 (0.960–1.427) | 0.0186 | 54.0% | 10.97 | 0.052 |
| M4 | 0.1522 | −0.0620 to 0.3665 | 0.1273 | 1.164 (0.940–1.443) | 0.0237 | 57.5% | 11.95 | 0.035 |

**PRIMARY RESULT: pooled attenuation_total = −0.0248, 95% CI −0.1253 to 0.0757,
p = 0.554.** τ² = 0.0052, I² = 44.4%, Q = 8.54 (df 5, p = 0.129), prediction
interval −0.2356 to 0.1860. Wald CI −0.1163 to 0.0668; fixed-effect −0.0081.

Components pooled: `attenuation_purity` = −0.0227 (−0.1080, 0.0626), p = 0.524,
I² 37.9%. `attenuation_stroma` = 0.0011 (−0.0252, 0.0275), p = 0.917, τ² = 0,
I² = 0%. **Each component uses its own bootstrap SE**, not the total's — an audit
caught them sharing one, which would have given wrong τ², CI, Q and I² for both.

The registered direction is `attenuation_total > 0` (adjustment reduces the
association). B.k states in advance that a negative attenuation "would be reported
as evidence against the thesis, not reframed."

**Contributing cohort sets are identical for M2 and M4** (all six), so Amendment 8
item 4's matched re-pool is not required in the primary. It **is** required inside
the alternative-endpoint sensitivity, where READ on OS gives 25 events and EPV 3.57
at M4 and drops out: matched over the five-cohort intersection (n = 1,590),
M2 = 0.0826 (−0.0066, 0.1718), M4 = 0.1013 (−0.0337, 0.2364).

### Proportional hazards — 8 of 25 model-cohort pairs violated

`p_score` < 0.05: STAD M4 (0.045), LIHC M1 (0.049). `p_global` < 0.05: COAD M2
(0.042) and M4 (0.0073), PAAD M3 (0.050), LIHC M1–M4 (M4 p = 9.2e-05).

The registered `score × log(time)` sensitivity is reported **alongside** the
primary, never replacing it. The interaction is significant only in LIHC M1
(β_tt = −0.133, p = 0.044); LIHC M2–M4 give p_tt = 0.102 / 0.100 / 0.094, and
COAD M2/M4, STAD M4 and PAAD M3 all p_tt > 0.26.

### VIF (M4, Fox–Monette GVIF)

| Cohort | score | purity | stromal_score |
|---|---|---|---|
| COAD | 2.02 | 2.00 | 2.61 |
| READ | 2.10 | 2.07 | 2.78 |
| STAD | 1.17 | 5.21 | 5.09 |
| ESCA | 1.51 | 7.96 | 8.03 |
| PAAD | 1.80 | 7.24 | 7.23 |
| LIHC | 1.39 | 1.77 | 1.68 |

Age, sex and stage are all below 1.2 in every cohort.

### B.4a Multiplicity (B.n)

Per-cohort estimates are secondary and descriptive. B.n adjusts them by
Benjamini–Hochberg **within model**, across the **six meta-analysed cohorts**;
CHOL is not in the family and receives `NA`. **No correction is applied across
M1–M4** — B.n states they are "a prespecified nested sequence addressing a single
question, not four independent hypotheses". The pooled estimate is untouched: it
is "one estimand, one test".

Of 24 per-cohort tests, **8 reach p < 0.05 raw and 3 survive BH at q = 0.05**:

| Cohort | Model | β | p raw | p BH |
|---|---|---|---|---|
| LIHC | M4 | 0.3626 | 0.00068 | **0.0041** |
| LIHC | M3 | 0.3215 | 0.00215 | **0.0129** |
| PAAD | M3 | 0.4578 | 0.00681 | **0.0204** |
| LIHC | M1 | 0.2331 | 0.00841 | 0.0505 |
| PAAD | M1 | 0.2769 | 0.01886 | 0.0566 |
| PAAD | M4 | 0.4111 | 0.01767 | 0.0530 |
| LIHC | M2 | 0.2339 | 0.00967 | 0.0580 |
| PAAD | M2 | 0.2754 | 0.02299 | 0.0690 |

Verified by independent recomputation of `p.adjust(..., method = "BH")` within
each model, and by asserting that no adjusted value falls below its raw value and
that no descriptive cohort received one. Adding B.n left **every other committed
number bit-identical** — 48 output files compared; only `survival_per_cohort.csv`
and its CHOL subset changed, and only by gaining the two new columns.

### Registered sensitivities

| Sensitivity | M2 pooled | M4 pooled | k |
|---|---|---|---|
| Alternative endpoint | 0.0814 (0.0091, 0.1536) | 0.1013 (−0.0337, 0.2364) | 6 / **5** |
| ESTIMATE purity in place of CPE | 0.1229 (0.0001, 0.2457) | 0.1519 (−0.0692, 0.3729) | 6 |
| 143-gene score | 0.1221 (0.0042, 0.2400) | 0.1534 (−0.0595, 0.3663) | 6 |
| Redaction-excluded | 0.1224 (0.0002, 0.2446) | 0.1517 (−0.0618, 0.3653) | 6 |
| Drop all low-EPV (READ) | 0.1382 (0.0019, 0.2745) | 0.1723 (−0.0872, 0.4317) | 5 |

Drop-low-EPV attenuation_total = −0.0254 (−0.1501, 0.0992). Leave-one-cohort-out
on attenuation_total ranges from −0.0763 (minus STAD) to +0.0180 (minus LIHC).
Fixed-effect cross-check is reported for every pool.

## B.5 Null-signature benchmark (09)

Registered B.m: N = **10,000 per cohort**, set size **152**, each panel gene
matched to a random gene from the same **decile of mean log2 expression** and the
same **decile of expression variance within that cohort**, sampled without
replacement within a set, seed `withr::with_seed(20260731 + 1000 + cohort_index)`.

Three configurations, each with a role declared **before** the run so that two
full-N p-values could not arrive with equal standing:

| Config | Role | N |
|---|---|---|
| 152 genes, expression × variance | **primary_registered** | 10,000 |
| 140 genes, expression × variance | size_sensitivity | 10,000 |
| 140 genes, + stromal correlation | **EXPLORATORY_POSTHOC** | 1,000 |

### Empirical p-values

Registered construction: `p = (1 + #{|β_null| ≥ |β_obs|}) / (1 + N)`, two-sided,
add-one. Granularity 9.999e-05 at N = 10,000.

| Config | Statistic | Observed | p two-sided (PRIMARY) | hits | p one-sided | tail | hits |
|---|---|---|---|---|---|---|---|
| **primary** | p_crude_M1 | 0.1212 | **0.3953** | 3,952 | 0.3953 | upper | 3,952 |
| | p_adjusted_M2 | 0.1229 | **0.4731** | 4,730 | 0.4731 | upper | 4,730 |
| | p_atten | −0.0226 | **0.1157** | 1,156 | 9.999e-05 | lower | **0** |
| size sens. | p_crude_M1 | 0.1212 | 0.2652 | 2,651 | 0.2652 | upper | 2,651 |
| | p_adjusted_M2 | 0.1229 | 0.3607 | 3,606 | 0.3607 | upper | 3,606 |
| | p_atten | −0.0226 | 0.1266 | 1,265 | 9.999e-05 | lower | 0 |
| *exploratory* | p_crude_M1 | 0.1212 | 0.7213 | 721 | 0.7213 | upper | 721 |
| | p_adjusted_M2 | 0.1229 | 0.8551 | 855 | 0.8551 | upper | 855 |
| | p_atten | −0.0226 | 0.0230 | 22 | 0.0230 | lower | 22 |

The two-sided and one-sided counts for `p_atten` differ because the observed value
sits **below all 10,000 nulls** (lower tail 0) while 1,156 nulls exceed it in
absolute value on the positive side. Both are correct and measure different things.

### Null distributions (registered primary)

| Statistic | min | q25 | median | mean | q75 | max | Observed | Percentile |
|---|---|---|---|---|---|---|---|---|
| pooled M1 log-HR | 0.0650 | 0.1079 | 0.1175 | 0.1175 | 0.1270 | 0.1681 | 0.1212 | **60.5** |
| pooled M2 log-HR | 0.0693 | 0.1122 | 0.1217 | 0.1219 | 0.1317 | 0.1780 | 0.1229 | **52.7** |
| pooled attenuation | −0.0080 | 0.0112 | 0.0153 | 0.0152 | 0.0195 | 0.0369 | −0.0226 | **0.0** |

**Proportion of null signatures whose pooled M1 CI excludes 1:** 28.1% (primary),
17.5% (140-gene), 97.2% (exploratory).

### Per-cohort percentile of the real panel (registered primary, M1)

| Cohort | Observed β | Null median | Null 2.5–97.5% | Percentile |
|---|---|---|---|---|
| COAD | 0.0052 | 0.0335 | −0.0287 to 0.0979 | 18.3 |
| READ | −0.0351 | 0.0397 | −0.0489 to 0.1385 | 4.9 |
| STAD | 0.1131 | 0.0994 | 0.0436 to 0.1561 | 67.9 |
| ESCA | 0.0542 | −0.0250 | −0.0798 to 0.0324 | **99.8** |
| PAAD | 0.2769 | 0.2150 | 0.1296 to 0.2952 | 93.2 |
| LIHC | 0.2331 | 0.2847 | 0.2127 to 0.3553 | 7.7 |
| CHOL | −0.1606 | −0.0717 | −0.1463 to 0.0032 | 1.0 |

Matching held in every cohort: realised set 152/152 everywhere, median set-level
difference 0.05–0.19 SD, paired per-gene 0.29–0.36 SD, **zero cell-widening**.

## B.6 CMS orthogonality (09, B.o)

Classifier: **CMScaller 2.0.1** (Eide et al. 2017), `templates.CMS`, single-sample
NTP, nPerm 1,000, seed 20261731, on log2(TPM+1) with Entrez rownames mapped from
the templates' own symbol column. Unclassified samples are retained as a level,
never dropped.

| Cohort | CMS1 | CMS2 | CMS3 | CMS4 | Unclassified |
|---|---|---|---|---|---|
| COAD | 18.2% (83) | 34.5% (157) | 16.0% (73) | 27.5% (125) | 3.7% (17) |
| READ | 9.7% (16) | 31.5% (52) | 20.0% (33) | 32.7% (54) | 6.1% (10) |

**Score tertiles × CMS1–4:** COAD χ² = 122.93 (df 6, p = 3.9e-24), READ χ² = 38.80
(df 6, p = 7.8e-07); simulated p = 9.999e-05 in both (10,000 replicates, used
because READ's minimum expected count is 5.16).

**R² of score on CMS membership alone:** COAD 0.286, READ 0.238.

**r(score, CMS4 membership):** COAD 0.405, READ 0.440. Mean score in CMS4 vs
elsewhere: COAD 0.658 vs −0.249; READ 0.629 vs −0.306.

**With CMS added as a covariate:**

| Cohort | Model | β without CMS | β with CMS | Δβ | Params | EPV | Fitted |
|---|---|---|---|---|---|---|---|
| COAD | M2 | 0.0046 | −0.1119 | −0.1165 | 9 | 11.33 | yes |
| COAD | M4 | −0.0713 | −0.1371 | −0.0658 | 11 | 9.27 | yes (flagged) |
| READ | M2 | −0.0457 | — | — | 9 | 4.00 | **no** |
| READ | M4 | −0.0221 | — | — | 11 | 3.27 | **no** |

READ is not fitted because `+ cms` adds four parameters and pushes EPV below
Amendment 8's floor of 5. COAD attenuation: 0.0759 without CMS, 0.0252 with.

**Stratified by CMS4 vs non-CMS4:**

| Cohort | Stratum | Model | n | Events | EPV | β | HR | Fitted |
|---|---|---|---|---|---|---|---|---|
| COAD | CMS4 | M1 | 125 | 32 | 32.0 | −0.2318 | 0.793 | yes |
| COAD | CMS4 | M2 | 125 | 32 | 6.4 | −0.0542 | 0.947 | yes (flagged) |
| COAD | CMS4 | M4 | 125 | 32 | 4.6 | — | — | **no** |
| COAD | non-CMS4 | M1 | 330 | 70 | 70.0 | −0.0460 | 0.955 | yes |
| COAD | non-CMS4 | M2 | 330 | 70 | 14.0 | −0.0669 | 0.935 | yes |
| COAD | non-CMS4 | M4 | 330 | 70 | 10.0 | −0.0674 | 0.935 | yes |
| READ | CMS4 | M1 | 54 | 14 | 14.0 | −0.4580 | 0.633 | yes |
| READ | non-CMS4 | M1 | 111 | 22 | 22.0 | 0.0347 | 1.035 | yes |

READ's stratified M2 and M4 are not fitted in either stratum (EPV 2.0–4.4).

---

# CROSS-CUTTING

## C.1 All 15 amendments

Every amendment was supplied verbatim by the study author and inserted without
editing. Full text in `panel_definition.md`; the "direction of bias" column
records what each amendment itself states, not a later gloss.

| # | Date | What it changed | Why | Stated direction of bias |
|---|---|---|---|---|
| 1 | 2026-07-31 | Criterion B restated from "at least two independent lines of evidence" to an explicit enumeration | The original wording admitted ~1,400 genes; the enumeration is what was intended | Narrows the panel. Neutral with respect to outcome — applied before any expression data was opened |
| 2 | 2026-07-31 | Criterion B(i) requires **human** ChIP-seq; mouse does not satisfy it | Cross-species TF binding is not evidence of a human direct target | Removes BCL2, MMP9, HGF from the panel. Conservative: it excludes three of the six genes that motivated the study |
| 3 | 2026-07-31 | §5 thresholds restated as **proportions** of the final panel | Written anticipating 30–50 genes; the panel is 152, which would have made k ≥ 8 a 5% bar rather than the ~20% intended | Raises the bar. Makes the BUILT branch harder to reach |
| 4 | 2026-07-31 | Estimand changed from **malignant-epithelial** to **epithelial** fraction, on tumour-channel samples only | Malignancy is inferred by three different methods across atlases; the epithelial lineage label is comparably defined | Removes a cross-atlas method asymmetry. Direction on k not predictable a priori |
| 5 | 2026-07-31 | Replication requirement: "at least two of the five atlases, **of any GI tissue**" (was tissue-matched) | Tissue-matching would have made most genes unevaluable | Loosens replication; increases k |
| 6 | 2026-07-31 | Disclosure + sensitivity for the 30–70% band rule. **Primary rule unchanged** | The band choice is consequential and should be visible | None — disclosure only. Realised: whole CI selects one branch, k_50 agrees |
| 7 | 2026-07-31 | Endpoint designation as a single deterministic rule from CDR Table 3 | Removes any per-cohort judgement | READ becomes PFI, six become OS. Event counts play no part, so no outcome-dependent selection |
| 8 | 2026-08-01 | EPV rule per **model**: ≥10 fit+pool; 5–10 fit+pool+flag; <5 not fitted. Cohort-level exclusion for cohorts too small for all endpoints | Prevents unstable small-cohort fits driving the pool | Low-EPV cohorts are **included** in the primary — flagging, not excluding, avoids conditioning the pool on precision. CHOL excluded at cohort level |
| 9 | 2026-08-01 | **GSE155698 removed**; five atlases → four | Deposit contains only CellRanger output — no annotation table, no metadata, no annotation column in any of 41 samples | Reduces replication breadth. Amendment 4's premise that labels were "available in all five" was factually wrong for this atlas |
| 10 | 2026-08-01 | **GSE183904 removed**; four → three. No gastric or oesophageal atlas remains | 40 flat count matrices, no annotation row or cluster file in any sample | Same premise failed a second time. STAD and ESCA are now modelled with no tissue-matched atlas — a stated limitation |
| 11 | 2026-08-01 | Peng compartment mapping uses `celltype1`, not `celltype0` | celltype0's haematopoietic label does not separate myeloid from lymphoid | **No bias claimed, and verified:** epithelial membership is bitwise identical across levels (24,156 cells either way), so neither numerator nor denominator of f(π) changes |
| 12 | 2026-08-01 | GSE125449's two sets combined on the **intersection** of their gene universes | Set1 has 20,124 gene rows, Set2 19,572; the script's own guard halted on the mismatch | 18,367 genes retained; 143 of 152 panel genes survive in that atlas. Six panel genes lost. Conservative: fewer genes evaluable, not more |
| 13 | 2026-08-01 | Peng raw counts **recovered by inverting the deposited normalisation** | The deposit has no count layer; `raw/X` is log1p-CP10K, verified as `sum(expm1(x)) = 10000` per cell | Recovery verified exact: recovered sums equal the deposited `n_counts` to 0.000e+00, max deviation from integers 1.17e-02, far below the 0.5 at which rounding could err |
| 14 | 2026-08-01 | k and variants computed over the **locked 152-gene panel**, not the final scoring list | k describes the panel's compartment behaviour; the scoring list is a different object serving a different purpose | Keeps k independent of Part B's exclusion rules. **Its stated reason is partly false** — see §C.4 |
| 15 | 2026-08-02 | Exclusion rule 3.3 on its **literal** reading: any gene annotated chrX or chrY, including pseudoautosomal. Final list 140 | Settles an ambiguity the author's earlier justification had gotten backwards | Removes 3 more genes (140 vs 143). Explicitly relies on neither the assistant's withdrawn dosage argument nor the auditor's replacement objection |
| 16 | 2026-08-02 | **Specifies script 10 (external validation)**, absent from analysis_plan.md v1.5 which ends at B.o. Four cohorts fixed in advance; FU-iCCA phosphoproteomic concordance designated primary | The original six-gene score was validated against RPPA in the same dataset its genes were selected from; this removes that circularity | **Made after the Part B primary result was known — disclosed, not absorbed.** Every parameter fixed before any validation data is opened. No pooling with discovery. Stated limitation: the phosphosite is bulk, so it can show the score tracks STAT3 phosphorylation but not which cells it occurs in |

## C.2 Every audit pass

The Implementation Auditor ran before every script executed, as standing practice
from Part A onward. Six passes; every one returned findings.

| Pass | Target | Verdict | Blocking | Total | What was fixed |
|---|---|---|---|---|---|
| 1 | `03_compartments.R` | fix_before_running | 2 | 19 | Guards that could not fire (a set operation deduplicated its own input; an unreachable exclusion path); a positional file round-trip that never asserted sample identity; a cohort constant validated by count rather than name |
| 2 | `03_compartments.R` re-audit | fix_before_committing | 0 | 8 | All 10 prior repairs verified correct; the auditor **retracted** its own earlier concern about the dominance-indicator CI after confirming A.f prespecifies it |
| 3 | `04_lock_gene_list.R` | fix_before_proceeding | 0 | 8 | A dead rule-3.1 path that could only halt; the PAR biology objection (see §C.4); NA-handling in the k denominator |
| 4 | `05` + `06` | fix_before_running | 2 | 18 | The purity conversion folding back below ESTIMATEScore −4121.5; a missing download being converted into a registered analytic branch |
| 5 | `07_score.R` | fix_before_running | 2 | 8 | Cohort-dependent duplicate-symbol collapse in the primary exposure; the stromal subscore readmitting excluded genes |
| 6 | `08_survival.R` | fix_before_running | 4 | 13 | PH sensitivity silently failing on zero-time events; the `do_not_fit` band never enforced; components sharing the total's bootstrap SE; `prop_attenuated` missing its sign condition |
| 7 | `09_null.R` | fix_before_running | 2 | 17 | The null comparison not being like-for-like; two full-N configurations shipping with equal standing and no declared primary |

**What was rejected, and why:**

- **Deriving EPV denominators from the realised fit** (pass 6). The auditor
  proposed counting parameters from `sum(!is.na(coef(fit)))` so a cohort with an
  empty stage level would count one fewer. Rejected: `analysis_plan.md:1186-1187`
  fixes the counts at M1=1, M2=5, M3=6, M4=7, and deriving them would make the EPV
  band depend on realised level occupancy — a cohort could cross a band boundary
  because one stage cell happened to be empty.
- **"The matching gate cannot fail"** (pass 7). The auditor simulated an unmatched
  draw at 0.103 SD against a 0.25 SD threshold and concluded the stop condition
  was vacuous. My own measurement on COAD gives **median 1.62 SD** for unmatched
  draws (p95 1.74, max 1.84) — the gate fires correctly. The underlying point was
  still worth acting on, so a **paired per-gene** gate was added (matched 0.32 SD
  vs unmatched 1.84 SD), and it immediately showed the exploratory scheme cuts
  stromal mismatch from 1.52 to 0.25 SD — a difference the set-level statistic
  could not see.
- **The PAR sex-chromosome argument** (pass 3) was *accepted against my own prior
  position*: I had retained four pseudoautosomal genes on a dosage argument
  asserted from memory, and the auditor showed it was backwards at its operative
  step (XCI-escape genes are typically expressed *higher* in females — that is the
  mechanism, not an exception). I withdrew the justification. Amendment 15 later
  settled the question on the registered text alone, relying on neither argument.

## C.3 Deviations from run orders, and which text governed

| Deviation | Run order said | What governed | Resolution |
|---|---|---|---|
| **B.m parameters** (script 09) | N ≈ 1,000; set size 140; matching on mean expression × stromal correlation | **Registered B.m** | Author confirmed the run order had misstated three parameters from memory. N = 10,000, size 152, matching on expression × variance deciles within cohort. A first draft written to the run order was discarded and its audit cancelled mid-flight. `NOTES_FOR_REVIEW.md` §29 |
| **`list_id` label** (script 07) | assert `list_id == "final_140"` | **The committed artefact** | The locked file says `primary_140`. The instructed label was erroneous (author confirmed). The artefact was not edited to make an assertion pass; identity is enforced by an md5 of the sorted gene symbols. §24 |
| **Function sources** (script 09) | "reuse the functions from 06 and 07" | **What the code path requires** | The scoring functions are in 07, but the model functions are in 08 — 06 computes purity and defines nothing the null path needs. 06 + 07 would never reach the models. Documented in the source. §32 |
| **Amendment 11 text** | arrived as a literal `[paste]`, then "use your best judgment" | **The rule that amendment text is verbatim** | I declined to compose prespecification text and proceeded with the other items until the real text arrived |

## C.4 Findings whose reasoning is not obvious from the code

These are recorded in `NOTES_FOR_REVIEW.md` and are the items a reader of the
scripts alone would miss.

**1. The ESTIMATE purity conversion folds back (§21).** `cos(a + bS)` is monotone
in S only while the *angle* lies in [0, π/2], i.e. ESTIMATEScore ∈ [−4121.5,
+6579.6]. Below −4121.5 the cosine turns back: purity *decreases* as the tumour
gets purer, and two scores map to one purity (S = −6000 and S = −2243 both give
0.962223). A guard testing the **cosine** for [0,1] cannot catch this because the
folded values are in range; the guard must test the **angle**. Verified latent on
this data — all 1,793 samples span −3205.9 to +5371.5, none in the fold region —
and fixed so it cannot become active silently in a validation cohort.

**2. The null comparison was not like-for-like (§30).** 08's pooled attenuation
was meta-analysed with a **paired bootstrap SE**; every null signature is pooled
with `sqrt(se₂² + se₄²)`, which ignores the positive covariance between β₂ and β₄
and is **1.72× larger** on this data (0.1835 vs 0.1069). Different weights give a
different pooled *point* estimate, not merely a different interval. `p_atten` is
therefore computed against the observed value re-pooled with the nulls' own
estimator (**−0.022577**); 08's −0.024766 remains the reported estimate. The
difference is 0.0022 log-HR units — small, but this is the decisive analysis and
the two numbers are not interchangeable.

**3. The stromal fallback moots an open specification question (§23).** The
auditor objected that the panel-derived subscore was built over the locked 152
genes and so readmitted genes the exclusion rules had removed — six in total, five
of them excluded *only* under rule 3.3, including TIMP1 (chrX, yet detected in
98.4% of stromal cells). The registered text says "panel genes" (152); the
auditor read 140. **Both were computed**: 103 vs 97 genes, correlating at r ≥
0.997, and both exceed 0.9 in every cohort, so both trigger the same fallback and
neither reaches the model. The question does not need deciding for Part B — but it
would if the fallback were ever not triggered, so `stromal_score_140` is retained.

**4. Amendment 14's stated reason is partly false (§20).** It says the exclusion
rules "did not exist when k was defined" and are "derived from the compartment
output". Both clauses fail: the rules existed at the first commit, and 6 of the 9
exclusions then in force *were* derived from Part A output. The conclusion — that
k belongs to the panel, not the scoring list — holds on the different ground that
they are distinct objects serving distinct purposes.

**5. Amendment 13's stated deviation is 4× low (§16).** The registered text gives
a figure for the maximum deviation from integers under count recovery that is
about four times smaller than the full-corpus measurement (1.17e-02). The recovery
still passes its registered 0.1 gate by a wide margin. Registered text was not
edited; the discrepancy is recorded.

**6. Half the rule-2 exclusions rest on one atlas (§19).** DNTT, LEP and CRLF2
were judged undetectable on a single atlas each, while 138 of the 140 retained
genes were held to all three. Two mechanisms compound: Amendment 12 removes
GSE125449 for some genes, and the evidence floor then removes Peng.

**7. GSE178341 has no endothelial compartment (§15).** Amendment 4 maps each
atlas's own level-1 labels, and GSE178341 has none for endothelium. The code
matches the registration and the estimand is unaffected — subdividing
non-epithelial mass changes neither numerator nor denominator — but the *reported*
breakdown of non-epithelial compartments differs between atlases.

**8. Bootstrap convergence warnings are about a nuisance term (§27).** The 08 run
emitted 50+ "coefficient may be infinite" warnings, all traceable to
`stage_groupmissing` in READ (9 patients, 1 event). In 300 test resamples **zero**
had any |coefficient| > 20 and the score coefficient stayed within [−0.79, 0.90].
`boot_nuisance_unstable` now records it per cohort (READ 1/2000, PAAD 93/2000,
others 0).

**9. Warn-level reviewer findings never reach the working conversation (§14).**
Ten accumulated unaddressed during Part A while work continued on top of them,
including two bearing on a numerical claim underpinning Amendment 13. This is now
`CLAUDE.md` §9 constraint 5: `host.findings()` is queried explicitly at the end of
every work item and all findings reported regardless of severity.

**10. Corrections to my own reporting (§25, §31, §32).** Recorded because a reader
comparing the transcript to the outputs would otherwise find discrepancies: the 05
unit-test count is 23 not 21; HGF's GSE125449 f(0.30) is 0.104 not 0.011 (a
chat-only transcription error — the CSV was always correct); a "GENCODE v36" claim
was inferred from a naming convention, not read, and is withdrawn; a
three-cohort spot check was described as covering the realised data; and the
code-path identity guard in 09 was **tautological** — it compared a binding to the
object it had just been copied from, and now compares function bodies against the
freshly-sourced committed files with a negative control.

---

# WHAT REMAINS

## R.1 Script 10 — external validation (not started)

**Now specified by Amendment 16** (2026-08-02), which supplies the section
`analysis_plan.md` v1.5 lacked. Four cohorts fixed in advance: FU-iCCA (NODE
OEP001105), GSE39582, GSE66229 and ICGC PACA-AU/PACA-CA. The **primary** validation
analysis is FU-iCCA phosphoproteomic concordance — the transcriptomic score
correlated against the directly measured STAT3 pY705 phosphosite — because the
original six-gene score was validated against RPPA in the same dataset its genes
came from, and that circularity is what this removes. Survival replication in the
other three is secondary. Purity is ESTIMATE-derived throughout; the 140-gene list
is primary with 143 as sensitivity; **no pooling with the TCGA discovery cohorts**.

The amendment discloses that it was made after the Part B primary result was
known — unavoidable, since the plan contained no validation section — and fixes
every parameter before any validation data is opened.

Scripts 05–08 were written as pure functions precisely so a validation cohort
reuses the identical code path, so the machinery is ready. **The data is not**:
see `DATA_NEEDED.md`. Two of the four cohorts need credentials a human must
obtain, and the primary cohort has the hardest access path of the four.

## R.2 Registered but not yet run

| Item | Section | Status |
|---|---|---|
| **B.n multiplicity — Benjamini–Hochberg FDR** | B.n | **IMPLEMENTED 2026-08-02.** `survival_per_cohort.csv` now carries `p_adj_BH` and `multiplicity_family` beside the raw `p`. See §B.4a |
| **`CMSclassifier` concordance check** against CMScaller's calls | B.o | Deferred. The package is not installed; CMScaller 2.0.1 alone produced the reported calls |
| **Score distribution by CMS as a figure** | B.o.1 | Deferred to the figures script. The underlying cross-tabulation and χ² are computed and committed |

Neither B.o deferral gates any reported result — B.o "does not gate the primary
analysis" by its own text.

## R.3 Figures

No figure script exists. Nothing in Parts A or B depends on one; every figure is a
presentation of a committed CSV. The registered figure obligations are B.o.1
(score distribution by CMS) and the compartment sweep figure deferred from A.f
under Amendment 12's dual-grouping requirement
(`output/gse125449_compartment_counts_dual.csv` and
`gse125449_stromal_fraction_dual.csv` hold its inputs).

## R.4 Reviewer findings still open

31 unresolved at the close of 10: 1 fail/medium, 7 warn/medium, 23 warn/low.
(An earlier draft of this section gave the split as 5 medium / 18 low at a count
of 24; the total was right but the split was not traced — corrected here and
recorded in `NOTES_FOR_REVIEW.md` §25.) Every finding with substantive content is
recorded in §14, §25 and §29–37; the remainder are anchored to sub-agent frames or
predate an exhausted review budget and cannot be cleared from this session. `NOTES_FOR_REVIEW.md` is
the durable record, not the findings panel.

---

# REPRODUCTION

```
Rscript 02_panel.R          # panel -> data/panel/panel_locked.csv (152)
Rscript 03_compartments.R   # A.a-A.g, ~16 min, B=2000, seed 20260731
Rscript 04_lock_gene_list.R # -> final_gene_list_140.csv / _143.csv
Rscript 05_clinical.R       # -> clinical_analysis_set.csv (1,793)
Rscript 06_purity.R         # CPE + ESTIMATE -> purity_per_patient.csv
Rscript 07_score.R          # -> scores_per_patient.csv
Rscript 08_survival.R       # M1-M4, bootstrap, meta -> ~16 min
Rscript 09_null.R           # N=10,000 x 7 x 3 configs + CMS -> ~27 min
```

Each script halts on any assertion failure rather than warning. 03 verifies its
own point estimates against a committed per-gene dominance signature; 09 verifies
that it reuses 07's and 08's function definitions by comparing function bodies
against the freshly-sourced committed files.

**Environment:** R 4.4 via conda env `stat3-gi`. The `estimate` package (v1.0.13)
is loaded from the system library at `/opt/homebrew/lib/R/4.6/site-library` —
appended to `.libPaths()`, never prepended, because prepending shadows `readxl`
with a build linking a second OpenMP runtime and aborts the session.
`CMScaller` 2.0.1 is installed at `renv/library-local`.

**Seeds:** 20260731 throughout (bootstrap and null draws), with documented
per-cohort and per-configuration offsets.
