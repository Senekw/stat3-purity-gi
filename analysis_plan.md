# Analysis plan — compartment decomposition of STAT3 activity scores in GI cancers

**Version:** 1.0
**Written:** 2026-07-31
**Author:** Sean GP Lee
**Panel:** locked 2026-07-31 at commit `ac9c5e0`; 152 genes; see `panel_definition.md`
**Status:** this is the document to be preregistered. Nothing in Part B is run
until this plan is registered.

Companion documents, all binding:
- `panel_definition.md` — panel prespecification and Amendments 1–4. LOCKED.
- `feasibility_assessment.md` — prior-art review and pilot compartment results.
  Authoritative for all pilot numbers.
- `output/atlas_malignant_annotation_audit.csv` — how each atlas calls malignancy.
- `output/atlas_tumour_designation_audit.csv` — tumour-vs-normal designation routes.
- `NOTES_FOR_REVIEW.md` — open items not acted on.

## 0. Two specification gaps requiring a decision before Part A runs

These are recorded here rather than resolved, because resolving either one is an
amendment, not an implementation choice.

### 0.1 Amendment 3's "at least two tissue-matched atlases" is unsatisfiable for
### three of four tissues

Amendment 3 defines *k* as panel genes epithelial-dominant "in at least two
tissue-matched atlases". Atlas coverage by tissue:

| Tissue | Atlases available | Count |
|---|---|---|
| Pancreatic | GSE155698, Peng | **2** |
| Liver / biliary | GSE125449 (covers both iCCA and HCC) | 1 |
| Colorectal | GSE178341 | 1 |
| Gastric | GSE183904 | 1 |
| Oesophageal | none | 0 |

Amendment 4 notes that making Steele usable "giv[es] PDAC two tissue-matched
atlases as Amendment 3 requires" — correct, but only for PDAC. Under a strict
reading (two atlases of the *same* tissue), *k* can only ever be driven by
pancreatic data, and the other three tissues contribute nothing. Under a loose
reading (two atlases anywhere in the panel's tissue set), "tissue-matched" does no
work at all.

Three ways out, all requiring an amendment:
- **(a)** Redefine the replication unit as two atlases *of any GI tissue*, and
  report per-tissue results alongside. Loses tissue-specificity; keeps all data.
- **(b)** Keep two-same-tissue but add a second atlas per tissue. Candidates
  exist (the Sci Data 2026 integrated atlas for gastric and colorectal), but it
  has no malignant-cell annotation — which under Amendment 4's lineage-level
  primary estimand is **no longer disqualifying**. This is the option Amendment 4
  quietly opens up and it deserves explicit consideration.
- **(c)** Require dominance in ≥2 atlases *or* in the single atlas for tissues
  where only one exists, and label the latter unreplicated.

**Nothing in Part A depends on this choice except step (g).** Steps (a)–(f) can
run as specified; *k* cannot be computed until this is settled.

### 0.2 GSE178341's expected count of 129 is GSMs, not patients

`specimen_type == "T"` gives 129 GSMs, which matches the instruction. But GSMs in
this series are *channels*, not patients: the same specimen appears under multiple
`processing_type` values (`unosrted` [sic], `LiveMACS`). Pelka's cohort is 62
patients (28 MMRp + 34 MMRd). So the halt-check should assert 129 GSMs, while the
bootstrap unit (step f) must be **patient**, not GSM or cell. Sorted and unsorted
channels of the same specimen are not independent samples, and treating them as
such would understate the interval. Flagged rather than assumed.

---

# PART A — compartment sweep (`03_compartments.R`)

## A.a Sample inclusion, as explicit filters with halt-on-mismatch

Every filter below is implemented as a filter *plus* an assertion. The script
stops on any mismatch rather than proceeding with a silently wrong sample set.

```r
assert_n <- function(observed, expected, atlas, rule) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(sprintf(
      "HALT [%s]: rule '%s' kept %d samples, expected %d. Do not proceed.",
      atlas, rule, observed, expected))
  }
  message(sprintf("OK [%s]: %s -> %d samples", atlas, rule, observed))
}
```

### GSE155698 — pancreatic (Steele)

Expected surviving samples: **17**.

```r
# Series is 41 GSMs: 17 PDAC_TISSUE, 3 AdjNorm_TISSUE, 17 PDAC_PBMC, 4 Healthy_PBMC.
keep <- grepl("^PDAC_TISSUE", sample_title)
```

**A filter keying on "PDAC" is WRONG and must not be used.** `grepl("PDAC", ...)`
matches the 17 `PDAC_PBMC` samples as well. Peripheral blood contains essentially
no epithelium, so including it would drive the epithelial fraction toward zero and
manufacture this study's own conclusion rather than test it. This is the single
most consequential filter in the plan.

Assertions: 17 kept; zero kept matching `PBMC`; zero kept matching `AdjNorm`;
`source_name` of all survivors equals `PANCREAS TUMOR`.

### GSE183904 — gastric (Kumar)

Expected surviving samples: **26**.

```r
# Labels are not clean: one GSM reads "Peritonium tissue  (Tumor)" (double space).
tissue_norm <- gsub("\\s+", " ", trimws(tissue_field))
keep <- tissue_norm == "Primary Gastric Tissue (Tumor)"
```

Whitespace normalisation is mandatory and must precede matching. Excluded: 10
`Primary Gastric Tissue (Normal)`, 3 `Peritonium tissue (Tumor)`, 1
`Peritonium tissue (Normal)` — peritoneal samples are excluded per Amendment 4
even when they are tumour.

Assertions: 26 kept; 40 GSMs seen in total; every one of the 40 normalised labels
matches one of the four known strings, else halt (an unrecognised label means the
series changed).

### GSE178341 — colorectal (Pelka)

Expected: **129 GSMs** with `specimen_type == "T"` (52 are `N`).

```r
keep_gsm <- specimen_type == "T"
# Per-cell propagation must be verified on open, not assumed:
if (!any(c("specimen_type","channel","sample_type") %in% colnames(cell_meta))) {
  stop("HALT [GSE178341]: no tumour/normal column in per-cell metatable. ",
       "Join via GSM before filtering; do not proceed.")
}
```

The published per-cell file is `GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz`.
If it carries no specimen column, the GSM-level designation must be joined onto
cells through the sample identifier before filtering — and the script halts rather
than guessing. See §0.2 on the GSM-vs-patient distinction.

### GSE125449 — liver, iCCA + HCC (Ma)

**No filter. Recorded as an explicit no-op.**

All 19 samples are tumours (10 iCCA, 9 HCC). The per-cell table carries only
`Sample`, `Cell Barcode`, `Type` — no tumour/normal field — but none is needed.

```r
# no-op by design; assert the premise rather than skipping silently
assert_n(n_distinct(samples$Sample), 19, "GSE125449", "no filter (all tumour)")
```

### Peng — pancreatic

Expected surviving samples: **24 of 35** (24 primary PDAC tumours, 11 control
pancreases).

```r
keep <- sample_group == "tumour"   # exact field name to be confirmed on open
```

## A.b Peng source decision

**Decision: use the primary deposit, GSA `CRA001160` / project `PRJCA001063`
(Genome Sequence Archive, BIG Data Center).**

Reasons. It is the accession named in the paper's own data-availability
statement, so the provenance chain is the authors' to the reader with no
intermediary. A Zenodo mirror, if one exists, is a third-party re-upload whose
processing (alignment version, filtering, cell calling, any re-annotation) is not
documented in the paper and would have to be characterised and reported
separately — added work and an added reviewer question, for no analytical gain.

If GSA access proves impractical (registration or download constraints) and a
mirror is used instead, this plan requires recording, in
`data/cache/provenance.txt`: the mirror URL and DOI, its upload date and
depositor, the stated relationship to CRA001160, file-level checksums, and any
processing the mirror applied. That record is a precondition of use, not a
follow-up.

Note the brief has referred to this as "the Peng Zenodo release"; the primary
deposit is GSA, and the plan follows the primary deposit.

## A.c Compartment harmonisation map

Target compartments, per Amendment 4: `epithelial`, `fibroblast_stromal`,
`myeloid`, `lymphoid`, `endothelial`, `other`.

Complete and verified for GSE125449 (label strings read from the downloaded
per-cell table):

| Atlas | Source label | Target compartment |
|---|---|---|
| GSE125449 | `Malignant cell` | epithelial |
| GSE125449 | `HPC-like` | epithelial |
| GSE125449 | `CAF` | fibroblast_stromal |
| GSE125449 | `TEC` | endothelial |
| GSE125449 | `TAM` | myeloid |
| GSE125449 | `T cell` | lymphoid |
| GSE125449 | `B cell` | lymphoid |
| GSE125449 | `unclassified` | other |

Note a change from the pilot: `TEC` (tumour endothelial cell) mapped to
fibroblast/stromal in `feasibility_assessment.md`, because that pilot used a
five-compartment scheme with no endothelial group. Amendment 4 makes endothelial
its own compartment, so `TEC` moves. Pilot stromal fractions are therefore not
directly comparable to Part A output, and the plan will report both groupings for
GSE125449 so the pilot remains auditable.

For the remaining four atlases the mapping rule is fixed now, but the exact source
label strings must be read from each file on open (see §3.1) and the map completed
before any pseudobulk is computed. The rule:

| Marker vocabulary in source label | Target |
|---|---|
| epithelial, malignant, tumour cell, ductal, enterocyte, goblet, Paneth, gastric/pit/chief/parietal, acinar, endocrine/islet, hepatocyte, cholangiocyte, HPC-like | epithelial |
| fibroblast, CAF, myofibroblast, stellate, pericyte, smooth muscle, mesothelial | fibroblast_stromal |
| macrophage, TAM, monocyte, neutrophil, granulocyte, mast, dendritic, myeloid | myeloid |
| T cell, CD4, CD8, Treg, NK, B cell, plasma, lymphocyte | lymphoid |
| endothelial, TEC, lymphatic endothelial, vascular | endothelial |
| unclassified, doublet, ambiguous, heterogeneous, enteric glial, Schwann | other |

**No unmapped labels are permitted.** The script builds the map, then:

```r
unmapped <- setdiff(unique(cell_meta$label), names(compartment_map))
if (length(unmapped)) {
  stop("HALT: unmapped level-1 labels: ", paste(unmapped, collapse = ", "),
       ". Extend the map in analysis_plan.md as a dated amendment, then re-run.")
}
```

Two mapping calls to record explicitly because they are judgement, not lookup:
- **Acinar and endocrine/islet cells → epithelial.** They are epithelial in
  lineage. In PDAC samples they are largely non-malignant residual parenchyma, so
  this inflates the epithelial fraction — conservative under Amendment 4's stated
  bias direction.
- **Enteric glial and Schwann cells → other**, not stromal. Neural crest lineage,
  not fibroblast.

## A.d Pseudobulk computation

Stated as code, because the order of operations is the whole point.

```r
# X : genes x cells, RAW UMI counts. No normalisation applied.
# cells: data.frame(cell_id, patient_id, atlas, compartment)

pseudobulk_raw <- function(X, cells, genes) {
  comps <- c("epithelial","fibroblast_stromal","myeloid","lymphoid",
             "endothelial","other")
  S <- matrix(0, nrow = length(genes), ncol = length(comps),
              dimnames = list(genes, comps))
  n_cells <- setNames(integer(length(comps)), comps)
  for (cc in comps) {
    idx <- which(cells$compartment == cc)
    n_cells[cc] <- length(idx)
    if (!length(idx)) next
    S[, cc] <- Matrix::rowSums(X[genes, idx, drop = FALSE])   # RAW counts summed
  }
  list(counts = S, n_cells = n_cells)
}

# Row-normalise AFTER summing: each gene's compartment shares sum to 1.
compartment_fraction <- function(S) S / rowSums(S)
```

**Forbidden: normalising per cell or per cell type before summing.** No CP10K, no
log, no scaling, no per-compartment mean, no library-size correction ahead of the
`rowSums`. The quantity of interest is each compartment's contribution to a bulk
library, and that contribution is proportional to cells × transcripts per cell.
Normalising first discards the transcript-per-cell term and silently reweights
every compartment to equal RNA content — the exact error that made the
HPA-based pilot figure measure the wrong estimand
(`feasibility_assessment.md`, addendum).

Genes absent from an atlas's annotation are recorded as `NA` for that atlas, never
as zero.

## A.e Purity sweep

The pseudobulk fraction from A.d reflects the *dissociated* composition of the
atlas, which is not the composition of a bulk tumour: dissociation recovers immune
cells more efficiently than epithelium. The sweep replaces the atlas's own
composition with a specified one.

Definitions. For gene *g* and compartment *c*, let `S[g,c]` be the raw summed
counts from A.d and `n_c` the number of cells in that compartment. Mean counts
per cell:

```
I[g,c] = S[g,c] / n_c
```

`I` retains between-compartment differences in total RNA content, which is
intended — a neutrophil and a carcinoma cell do not contribute equal transcript
mass, and that asymmetry is part of the bulk contribution.

Reweighting. For target epithelial cell proportion π, weights are

```
w_epithelial = π
w_c          = (1 - π) * n_c / sum(n_c over non-epithelial c)     for c != epithelial
```

That is, non-epithelial compartments retain their observed proportions *relative
to one another* and are jointly scaled to (1 − π). The epithelial fraction at π is

```
f[g](π) = w_epithelial * I[g,epithelial] / sum_c ( w_c * I[g,c] )
```

**Grid:** π from **0.30 to 0.70 inclusive, step 0.01** (41 points). This is the
band named in Amendment 3 and stated there as typical GI tumour purity.

**Evaluation of "epithelial-dominant across the full band":** gene *g* is
epithelial-dominant in an atlas iff `f[g](π) > 0.50` for **every** π on the grid.
Note `f[g](π)` is monotone increasing in π, so this is equivalent to
`f[g](0.30) > 0.50`; the full grid is still evaluated and reported, and the
monotonicity is asserted numerically rather than assumed:

```r
if (any(diff(f_grid) < -1e-9)) stop("HALT: f(pi) not monotone; check weights.")
```

The reported sweep figure is `f[g](π)` against π with the 0.50 line and the
30–70% band marked, one panel per atlas.

Purity source for context, not for the sweep grid: ESTIMATE and the Aran CPE
consensus table, two independent estimates, agreement between them reported as a
sensitivity analysis. The sweep itself is a specified grid and does not depend on
either.

## A.f Patient-level bootstrap

**Resampling unit: patient.** Not cell, not GSM, not channel. Cells are not
independent within a patient, and in GSE178341 several GSMs share a specimen
(§0.2). Resampling cells would produce intervals far too narrow.

- **Resamples:** B = 2000, within atlas, patients drawn with replacement to the
  observed patient count.
- **Seed:** `set.seed(20260731)`, recorded here so the interval is reproducible.
  One seed for the whole sweep; per-atlas streams derived by `withr::with_seed`
  on `20260731 + atlas_index` so adding an atlas later does not perturb existing
  atlases' draws.
- **Recomputed inside each resample:** the full A.d → A.e chain (raw sums,
  per-cell means, reweighting) on the resampled patient set. The bootstrap must
  wrap the whole pipeline, not just the final ratio.
- **Reported interval:** 95% percentile bootstrap CI (2.5th, 97.5th) on `f[g](π)`
  at every grid point, and on the dominance indicator, per gene per atlas.
- Genes with fewer than 20 summed counts across all compartments in an atlas are
  reported as insufficient-evidence rather than given a fraction; this threshold
  is prespecified and matches the pilot.

Bootstrap is also used for *k*: within each resample, *k* is recomputed under the
A.g rule, giving a distribution and a 95% interval for *k* itself. The branch
decision uses the point estimate; the interval is reported alongside, and if the
interval spans a branch boundary that fact is stated in the paper.

## A.g Computing k, per Amendment 3

*k* = number of panel genes that are epithelial-dominant (A.e) in at least two
tissue-matched atlases (**pending §0.1**).

Panel size is locked at 152, so Amendment 3's proportional thresholds are fixed
arithmetic:

- 20% of 152 = 30.4 → the `k >= 20%` condition is **k ≥ 31**
- 5% of 152 = 7.6 → the `5% <= k` condition is **k ≥ 8**

| Realised *k* | Branch |
|---|---|
| **k ≥ 31** (and ≥ 8, automatically satisfied) | epithelial subscore built; subscore survival models run as SECONDARY |
| **8 ≤ k ≤ 30** | subscore reported as EXPLORATORY only; no formal survival inference on it |
| **k ≤ 7** | decomposition DESCRIPTIVE only; no subscore |

The branch is read off the table and not renegotiated. All three outcomes are
publishable, and the primary endpoint (Part B) does not depend on which obtains.

The origin-score genes are reported as a labelled subset throughout: `SOCS3`,
`MYC`, `IL6` are in the panel; `BCL2`, `MMP9`, `HGF` failed criterion B and are
reported as a labelled non-qualifying subset per §4 of the prespecification and
Amendment 2. Non-qualifying genes do not contribute to *k*.

---

# PART B — primary survival analysis (`07_score.R`, `08_survival.R`, `09_meta.R`)

Part B is the **primary** analysis. Per §6 of the prespecification it requires no
compartment fractions and cannot degenerate: cohorts with no usable single-cell
atlas (ESCA, and CHOL for practical purposes) contribute fully here.

## B.h Score construction

```r
# se : SummarizedExperiment for one cohort, STAR-Counts
# 1. primary tumour samples only, one aliquot per patient
#    (sample_type == "Primary Tumor"; TCGA-CHOL and others carry normals)
# 2. expression matrix: tpm_unstrand, log2(x + 1)
# 3. drop genes with zero variance in this cohort
# 4. z-score EACH GENE WITHIN COHORT (mean 0, sd 1 across that cohort's patients)
# 5. score = rowMeans over the 152 panel genes present in the cohort
# 6. scale the resulting score to unit SD within cohort

score_cohort <- function(se, panel_genes) {
  x   <- log2(SummarizedExperiment::assay(se, "tpm_unstrand") + 1)
  x   <- x[rownames(x) %in% panel_genes, , drop = FALSE]
  x   <- x[matrixStats::rowSds(x) > 0, , drop = FALSE]
  z   <- t(scale(t(x)))                       # per-gene z within cohort
  s   <- colMeans(z, na.rm = TRUE)            # mean across panel
  as.numeric(scale(s))                        # per-SD scaling
}
```

Three points fixed here deliberately, all of which the ESMO score got differently
and which the discussion will address:

- **z-scoring within cohort before averaging.** The original six-gene score used
  `mean(log2(expr + 1))` with no z-scoring, so genes with larger dynamic range
  dominated the average arbitrarily. Z-scoring gives each panel gene equal weight,
  which is what "mean of a panel" is normally taken to mean.
- **Per-SD scaling.** All hazard ratios are therefore **per 1 SD of score within
  cohort**, which is what makes them poolable across cohorts of different scale.
  No dichotomisation, no tertiles, no optimised cutpoint — the ESMO abstract's
  unstated cutoff derivation is a known weakness this plan removes rather than
  repeats.
- **Genes missing from a cohort** are dropped for that cohort with the count
  reported; the score is the mean over genes present. Cohorts retaining fewer than
  140 of 152 genes are flagged in the results table.

The **stromal subscore** used as a covariate in B.j model 4 is constructed
identically from panel genes that are *not* epithelial-dominant. If Part A returns
the descriptive branch (k ≤ 7), the stromal subscore is defined on all non-dominant
panel genes, which is nearly the whole panel — in that case model 4 is collinear
with model 1 by construction, and the plan substitutes a purity-orthogonal stromal
index (ESTIMATE stromal score) instead. Stated now so the choice is not made after
seeing k.

## B.i Endpoint per cohort

Endpoints follow the per-cancer-type recommendations of the TCGA Clinical Data
Resource (Liu et al., *Cell* 2018), which is the point of using the CDR rather
than raw TCGA clinical files.

**The per-cohort table below is provisional and must be filled from Table 1 /
Table S1 of the CDR file before registration** — `TCGA-CDR.xlsx` is a manual
download and is not yet on disk (see §3.2). The general rule the CDR encodes: OS
is usable where event counts suffice; PFI is preferred for cancer types with long
survival and few deaths within follow-up; DSS and DFI are secondary.

| Cohort | Primary endpoint (provisional) | Basis |
|---|---|---|
| TCGA-COAD | PFI | CDR notes OS underpowered in colon; recorded in `feasibility_assessment.md` §5 |
| TCGA-READ | PFI | as COAD; rectal event counts lower still |
| TCGA-STAD | OS | event count adequate |
| TCGA-ESCA | OS | event count adequate |
| TCGA-PAAD | OS | high event rate, short survival |
| TCGA-LIHC | OS | CDR flags PFI as also usable; OS primary, PFI sensitivity |
| TCGA-CHOL | OS, **descriptive only** | n = 44 samples; not entered into the meta-analysis as a weighted stratum, reported separately |

One endpoint per cohort is designated primary in advance. The alternative endpoint
is reported as a prespecified sensitivity analysis, never substituted for the
primary if the primary is null.

Follow-up is administratively censored at 10 years for OS and 5 years for PFI, per
common CDR practice, and the censoring rule is applied identically across cohorts.

## B.j Model sequence

Four nested Cox proportional-hazards models per cohort, fitted in this order, on
the same patient set (complete cases across all covariates in model 4, so that all
four models are fitted on an identical set and the attenuation in B.k is not
confounded by changing n):

```r
# M1  unadjusted
coxph(Surv(time, event) ~ score)

# M2  + clinical
coxph(Surv(time, event) ~ score + age_at_diagnosis + sex + stage_group)

# M3  + purity
coxph(Surv(time, event) ~ score + age_at_diagnosis + sex + stage_group + purity)

# M4  + stromal score
coxph(Surv(time, event) ~ score + age_at_diagnosis + sex + stage_group + purity
                          + stromal_score)
```

Covariate handling, fixed in advance:
- `age_at_diagnosis`: continuous, years.
- `sex`: factor. Dropped in any cohort with < 10 patients of either sex.
- `stage_group`: AJCC stage collapsed to I/II vs III/IV. Collapsed rather than
  four-level because CHOL and PAAD have sparse cells; patients with missing stage
  are retained with an explicit `missing` level rather than dropped.
- `purity`: consensus CPE (Aran et al. 2015) as primary; ESTIMATE-derived purity
  as sensitivity. Continuous, on its native 0–1 scale.
- `stromal_score`: as B.h.

Proportional-hazards assumption tested by `cox.zph` on every model; violations
reported and, where present, handled by adding a `score × log(time)` interaction as
a prespecified sensitivity rather than by silently accepting the violation.

## B.k The estimand: attenuation in log-HR

This is the quantity the paper is about, and it is stated as a number with an
interval, not as a narrative comparison of two p-values.

Let β₂, β₃, β₄ be the score log-HR in M2, M3, M4. Define:

```
attenuation_purity  = beta2 - beta3                  (absolute, log-HR units)
attenuation_stroma  = beta3 - beta4
attenuation_total   = beta2 - beta4
prop_attenuated     = 1 - beta4 / beta2              (reported only if beta2 != 0
                                                      and sign(beta2)==sign(beta4))
```

**Primary estimand: `attenuation_total`**, the absolute reduction in score log-HR
on adding purity and stromal score to the clinically adjusted model, meta-analysed
across cohorts.

`prop_attenuated` is reported as a secondary, interpretable summary but is not the
primary, because it is unstable when β₂ is near zero and undefined when the sign
flips. Both are reported; the primary is the absolute difference.

**Interval:** by nonparametric bootstrap over **patients within cohort**,
B = 2000, `set.seed(20260731)`, refitting all four models in every resample so
that β₂ and β₄ are paired within resample and their difference carries the
correlation between the estimates. A naive interval from two independent standard
errors would be wrong — the estimates come from nested models on identical data
and are strongly positively correlated, so an independence assumption
overstates the width. 95% percentile CI reported.

Direction is prespecified: the thesis predicts `attenuation_total > 0` (adjustment
reduces the association). A negative attenuation — adjustment *strengthening* the
association — is a real possible outcome and would be reported as evidence against
the thesis, not reframed.

## B.l Meta-analysis across cohorts

- **Model:** random-effects, REML estimation of τ², on the per-cohort score log-HR
  (and separately on `attenuation_total`), inverse-variance weighted.
  `metafor::rma(yi, vi, method = "REML")`.
- **Small-k interval handling:** with six cohorts entering (CHOL descriptive only,
  per B.i) the standard normal-approximation CI is anticonservative. The
  **Hartung–Knapp–Sidik–Jonkman** adjustment is used for the pooled CI
  (`test = "knha"`), prespecified as the primary interval. The Wald interval is
  reported alongside for comparability with the literature.
- **Heterogeneity:** τ², I², Q with its p-value, and a prediction interval. I² is
  reported but not used as a decision rule — with six studies its own uncertainty
  is large.
- **Fixed-effect pooling** reported as a sensitivity analysis only.
- **Influence:** leave-one-cohort-out pooled estimates reported as a table. No
  cohort is dropped on the basis of influence.
- Funnel plots and small-study-effect tests are **not** performed: with six
  cohorts they have no power and would be decorative.

## B.m Null-signature procedure

Following Venet, Dumont & Detours (*PLoS Comput Biol* 2011), cited as the method
source and not claimed as a contribution.

- **Number of random sets:** N = 10,000 per cohort.
- **Set size:** 152 genes, matching the panel exactly (or the cohort's realised
  panel size where genes are missing, matched per cohort).
- **Matching variables:** each panel gene is matched to a random gene drawn from
  the same **decile of mean log2 expression** and the same **decile of expression
  variance** within that cohort, sampled without replacement within a set.
  Matching on both is necessary: unmatched random sets are dominated by
  low-expression, low-variance genes and produce an artificially weak null.
- **Seed:** `set.seed(20260731)`; the null draw is a separate RNG stream from the
  bootstrap (`withr::with_seed(20260731 + 1000 + cohort_index)`).
- **Statistic:** each null set is scored by the identical B.h pipeline and fitted
  through M2, giving a null distribution of the score log-HR.
- **Empirical p-value:**

```
p_emp = (1 + #{ |beta_null| >= |beta_observed| }) / (1 + N)
```

  Two-sided, add-one so p is never exactly zero. The same construction is applied
  to `attenuation_total`, giving an empirical p for the attenuation itself — the
  more relevant null here, since the question is whether *this* panel attenuates
  more than an arbitrary panel would.

## B.n Multiplicity across the seven cohorts

- The **primary inference is the meta-analytic pooled estimate**: one estimand,
  one test, no multiplicity correction required or applied.
- The **per-cohort estimates are secondary and descriptive**. Where they are
  tested, p-values are adjusted across the six meta-analysed cohorts by
  Benjamini–Hochberg FDR at q = 0.05, and both raw and adjusted values are
  reported in the same table.
- CHOL is reported descriptively and is not counted in the multiplicity family.
- The four models M1–M4 within a cohort are a prespecified nested sequence
  addressing a single question, not four independent hypotheses; no correction is
  applied across them, and this is stated.
- No subgroup analyses beyond those named in this plan. Any post-hoc subgroup that
  appears in the paper is labelled post-hoc explicitly.

## B.o CMS orthogonality analysis (CRC cohorts)

Purpose: pre-empt the objection that a stroma-weighted STAT3 score is CMS4
relabelled (Guinney et al. 2015; `prior_art_matrix.csv`).

- **Cohorts:** TCGA-COAD and TCGA-READ.
- **Classifier:** `CMScaller` (Bioconductor) on log2 TPM, which is designed for
  single-sample CMS calling; `CMSclassifier` as a concordance check. Samples the
  classifier leaves unlabelled are retained as an `unclassified` level, not dropped.
- **Analyses, all prespecified:**
  1. Cross-tabulation of score tertiles against CMS1–4 with a χ² test, and the
     score's distribution by CMS as a figure.
  2. Variance decomposition: R² of the score on CMS membership alone. If CMS
     explains most of the score's variance, the score is largely a CMS proxy and
     the paper says so.
  3. The B.j model sequence refitted **with CMS added as a covariate**, and
     separately **stratified by CMS4 vs non-CMS4**. Persistence of the score
     association within non-CMS4 patients is the evidence that the score is not
     simply CMS4.
- This analysis is confirmatory of orthogonality or it is not; either result is
  reported. It does not gate the primary analysis.

---

# 3. What could not be specified without computing or downloading

Listed because each is a real gap in this plan, not a formality. Nothing below was
guessed at in the text above; each is marked provisional where it appears.

## 3.1 Exact level-1 label strings for four of five atlases (A.c)

GSE125449's labels are exact — the per-cell table is on disk. GSE183904,
GSE178341, GSE155698 and Peng have their **mapping rule** fixed in A.c, but the
literal source strings must be read from each file. The published texts give
lineage vocabulary, not the annotation column's exact values, and Kumar's 34
lineage states and Pelka's 88 subsets are not enumerated verbatim in either paper's
main text. Mitigation: the halt-on-unmapped-label check in A.c makes an incomplete
map a loud failure rather than a silent miscount. Requires downloading the atlases,
which this session was instructed not to do.

## 3.2 Per-cohort endpoint recommendations (B.i)

`TCGA-CDR.xlsx` is a manual download and is not on disk. The table in B.i is
provisional: the COAD/READ→PFI assignment carries forward from
`feasibility_assessment.md`, and the others follow the CDR's general rule, but the
authoritative per-type recommendation must be read from the file. **This must be
resolved before registration**, since designating the primary endpoint after seeing
outcome data would void the preregistration.

## 3.3 Event counts, and whether any cohort is too small to model (B.i, B.l)

Which cohorts have adequate events for a five-covariate Cox model cannot be known
without the CDR merge. CHOL is designated descriptive-only on sample count alone
(n = 44). Whether READ (177 samples) supports model M4 with events-per-variable
above the conventional 10 is unknown until events are counted. The plan commits to
reporting events-per-variable per cohort and to flagging any model below 10 rather
than to a fixed exclusion rule — and that flagging rule is itself prespecified here.

## 3.4 Purity coverage (B.j)

`aran_purity.xlsx` is not on disk, so the proportion of patients in each cohort
with a CPE value is unknown. Since M3 and M4 require purity, missing CPE reduces
the complete-case set that all four models share (B.j). If coverage is poor in a
cohort, the ESTIMATE-derived purity becomes primary for that cohort — the rule is
prespecified as: CPE primary where coverage ≥ 80% of the cohort's patients,
ESTIMATE primary otherwise, with the switch reported per cohort.

## 3.5 GSE178341 per-cell specimen propagation (A.a, §0.2)

Whether the published metatable carries a tumour/normal column per cell, or whether
the GSM-level designation must be joined on, is unresolved. Handled by the halt
check rather than an assumption. Related: the patient identifier needed for the
A.f bootstrap unit must be confirmed present in the same file.

## 3.6 GSE183904 per-cell malignancy cutoff (Amendment 4 sensitivity analysis)

Amendment 4's sensitivity analysis recomputes the fraction restricted to malignant
cells "where the source atlas supplies CNV-based per-cell malignancy calls
(GSE125449, GSE183904)". Kumar's main text reports inferCNV scores compared between
epithelial and macrophage cells but does not state a per-cell malignancy threshold.
Whether GSE183904 ships per-cell malignant calls, or only a CNV score requiring a
threshold this project would have to choose, is unknown until the file is opened.
If the latter, the GSE183904 arm of the sensitivity analysis cannot proceed on the
atlas's own definition, and that limits the sensitivity analysis to GSE125449 —
which should be recorded as an amendment rather than absorbed silently.

## 3.7 Peng field names and sample-level tumour flag (A.a, A.b)

The exact metadata field distinguishing the 24 tumours from the 11 control
pancreases is written as `sample_group == "tumour"` pending confirmation from the
GSA deposit's own annotation. Also unverified: whether a Zenodo mirror exists at
all, which A.b requires characterising *if* it is used.

## 3.8 Whether Part A's k interval spans a branch boundary (A.f, A.g)

By construction unknowable before the sweep. The plan commits in advance to
reporting the fact if it happens, and to using the point estimate for the branch
decision, so that the disclosure is prespecified rather than discretionary.

## 3.9 Not a computation gap, but unresolved: §0.1 and §0.2

The Amendment 3 tissue-matching ambiguity (§0.1) blocks step A.g and needs a
decision, not a computation. The GSM-vs-patient distinction (§0.2) is resolved
within this plan (bootstrap on patients) but the expected-count assertion is
written on GSMs, and that mismatch of units should be read and confirmed rather
than discovered later.

