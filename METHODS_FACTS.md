# METHODS_FACTS — values for the manuscript methods section

Every value below is transcribed from a committed file, with the file named. **No
value was computed for this document.** Where a fact the methods section needs is
not in any committed file, it is marked **`[NOT FOUND]`** with what is there
instead — see §15 for the consolidated list of gaps.

Assembled 2026-08-05 at commit `353cc95`.

---

## 1. Panel construction

Source: `data/panel/panel_provenance.txt`, `data/panel/criterionB_provenance.txt`.

| Item | Value |
|---|---|
| Panel locked | 2026-07-31 |
| MSigDB access | via **msigdbr 26.1.0**; gene set `HALLMARK_IL6_JAK_STAT3_SIGNALING`. Underlying release **`msigdb.2026.1`**, read from the package's own hardcoded data URL (`msigdbr:::check_cache`) — see §15 for the caveat that this is the release the package targets, not an independent check on the 87 genes |
| ChEA3 version | **Not versioned by ChEA3 — this is the correct citation, not a gap.** Recorded fingerprint instead: file byte sizes, queried **2026-07-31** from `https://maayanlab.cloud/chea3/assets/tflibs/` (the public API is enrich-only, so GMTs were pulled directly) |
| TRRUST version | **v2 human**, queried **2026-07-31**, `https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv`, 297,659 bytes, 9,396 edges, 185 STAT3 edges → **142 unique STAT3 targets** (Activation 83, Repression 27, Unknown 75) |

**ChEA3 library fingerprints** (bytes / STAT3 terms / genes):

| Library | Bytes | STAT3 terms | Genes | Used for criterion B |
|---|---|---|---|---|
| ENCODE_ChIP-seq | 5,593,822 | 3 (HELAS3, GM12878, MCF10A; hg19) | 5,480 | yes |
| ReMap_ChIP-seq | 2,960,359 | 1 | 1,453 | yes |
| Literature_ChIP-seq (human) | 2,453,171 | 7 raw → 6 after duplicate collapse | 2,923 | yes (1 experiment: PMID 23295773, U87) |
| Literature_ChIP-seq (mouse) | — | 5 experiments | 3,995 | **no** — reported, not counted |
| ARCHS4_Coexpression | 3,237,629 | — | — | **excluded**: coexpression, not binding |
| GTEx_Coexpression | 3,186,651 | — | — | **excluded**: coexpression, not binding |
| Enrichr_Queries | 2,529,276 | — | — | **excluded**: user query lists, not primary |

One collapse applied: `STAT3_1855785_CHIPSEQ_MESC_MOUSE` removed as a
truncated-PMID duplicate of `STAT3_18555785_CHIPSEQ_MESC_MOUSE`.

**Counts.**

| Quantity | Value |
|---|---|
| Criterion A (HALLMARK_IL6_JAK_STAT3_SIGNALING) | **87** genes |
| Criterion B as originally written (≥2 independent sources, human-only) | 1,387 genes; 1,448 in union with A |
| Criterion B including mouse ChIP-seq | 2,974 genes |
| **Criterion B as amended** (Amendment 1: ≥1 human ChIP-seq source AND TRRUST v2 curation) | **75** genes |
| **Final panel** | **152** genes |

**Route breakdown** (`panel_provenance.txt`): both **10**, A_only **77**, B_only
**65**. 10 + 77 + 65 = 152.

**Original six.** Qualified: **SOCS3** (ENCODE, ReMap, TRRUST), **MYC** (ENCODE,
ReMap, TRRUST), **IL6** (ENCODE, Literature-human, TRRUST) — 3 human sources each.
Failed: **BCL2, MMP9, HGF**, TRRUST only. Recorded caveat: MMP9 and HGF each also
appear in mouse ChIP-seq and would pass under an include-mouse rule; BCL2 has no
ChIP-seq evidence in any group and fails under either counting scheme.

Pairwise Jaccard between the four retained evidence groups, as non-redundancy
evidence: ENCODE–ReMap 0.085, ENCODE–Literature 0.191, ENCODE–TRRUST 0.010,
ReMap–Literature 0.046, ReMap–TRRUST 0.016, Literature–TRRUST 0.012.

---

## 2. Exclusions — all 12, with the value that triggered each

Source: `data/panel/final_gene_list_exclusions.csv`. Rules are prespecification
§3: **rule 2** = below 1% of cells in every annotated compartment in every atlas;
**rule 3.3** = sex-chromosome gene (broad reading, Amendment 15).

| Gene | Rule | Trigger value | Chr | PAR | Atlases evaluable | Route |
|---|---|---|---|---|---|---|
| DNTT | 2 | max detection **0.027%** (GSE178341/lymphoid) | chr10 | no | 1 | A_only |
| GFAP | 2 | **0.5017%** (Peng/epithelial) | chr17 | no | 3 | B_only |
| LEP | 2 | **0.3032%** (GSE178341/myeloid) | chr7 | no | 1 | B_only |
| OPRM1 | 2 | **0.872%** (GSE178341/lymphoid) | chr6 | no | 2 | B_only |
| PAX3 | 2 | **0.1356%** (Peng/epithelial) | chr2 | no | 2 | B_only |
| CRLF2 | **2 + 3.3** | **0.1385%** (GSE178341/myeloid) *and* chrX/chrY | chrX/chrY | **yes** | 1 | A_only |
| IL13RA1 | 3.3 | chrX (detection 44.92%, not the trigger) | chrX | no | 3 | A_only |
| IL2RG | 3.3 | chrX (56.72%) | chrX | no | 3 | A_only |
| TIMP1 | 3.3 | chrX (98.38%) | chrX | no | 3 | B_only |
| CSF2RA | 3.3 | chrX/chrY pseudoautosomal (44.10%) | chrX/chrY | **yes** | 3 | A_only |
| IL3RA | 3.3 | chrX/chrY pseudoautosomal (61.10%) | chrX/chrY | **yes** | 3 | A_only |
| IL9R | 3.3 | chrX/chrY pseudoautosomal (4.17%) | chrX/chrY | **yes** | 2 | A_only |

CRLF2 is caught by **both** rules, which is why 12 exclusions (not 13) remove
12 genes: **152 − 12 = 140**.

- **Final list: 140** (`data/panel/final_gene_list_140.csv`, `list_id = primary_140`)
- **Sensitivity list: 143** (`data/panel/final_gene_list_143.csv`, `list_id = sensitivity_143`) — the narrow reading of rule 3.3, retaining the 3 non-PAR-only genes
- No origin-six gene is excluded under either list.

Two exclusions are worth stating in the text because they are not "undetectable"
in the ordinary sense: **GFAP is dominant in all 3 atlases** and **PAX3 in 2**, yet
both fall under rule 2. And three rule-2 genes (DNTT, LEP, CRLF2) were judged on a
**single** atlas.

---

## 3. Atlases

Declared in `03_compartments.R`; realised counts from `output/partA_run.err`.

| | GSE125449 | GSE178341 | Peng |
|---|---|---|---|
| Tissue | liver_biliary | colorectal | pancreatic |
| Publication | Ma et al. 2019 (GEO) | Pelka et al. 2021 (GEO) | Peng et al. 2019, Besca reprocessing |
| Accession / source | GSE125449 | GSE178341 | Zenodo **10.5281/zenodo.3969339**, `StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad` |
| Samples after filtering | **19** (all tumour, no filter) | **62** tumour patients | **24** tumour patients (11 normal excluded, 35 total) |
| Cells, deposit | 9,946 | 370,115 | 57,423 |
| **Cells analysed** | **9,946** | **257,251** (tumour only) | **41,964** (tumour only) |
| Bootstrap unit | patient (19) | patient (62) | patient (24) |
| Annotation column | `Type` | `clTopLevel` | **`celltype1`** (Amendment 11) |

**Compartments** (Amendment 4, six): `epithelial`, `fibroblast_stromal`,
`myeloid`, `lymphoid`, `endothelial`, `other`.

**Full label → compartment maps, as applied.**

GSE125449 (`MAP_GSE125449`):

| Label | Compartment |
|---|---|
| Malignant cell | epithelial *(judgement: malignant hepatocyte/cholangiocyte)* |
| HPC-like | epithelial *(hepatic progenitor-like)* |
| CAF | fibroblast_stromal |
| TEC | endothelial *(NB: fibroblast_stromal in the pilot's 5-compartment scheme)* |
| TAM | myeloid |
| T cell | lymphoid |
| B cell | lymphoid |
| unclassified | other |

GSE178341 (`MAP_GSE178341`):

| Label | Compartment |
|---|---|
| Epi | epithelial |
| Strom | fibroblast_stromal *(includes fibroblast/pericyte/endothelial subsets at the cl295 level)* |
| Myeloid | myeloid |
| Mast | myeloid *(judgement: granulocyte lineage, per A.c)* |
| TNKILC | lymphoid |
| Plasma | lymphoid |
| B | lymphoid |

**GSE178341 has no endothelial compartment** — endothelium is inside `Strom` at
this annotation level (recorded as a reporting limitation, `NOTES_FOR_REVIEW.md`).

Peng (`MAP_Peng`, keyed on **celltype1** per Amendment 11):

| Label | Compartment |
|---|---|
| pancreatic ductal cell | epithelial |
| pancreatic acinar cell | epithelial *(A.c judgement)* |
| enteroendocrine cell | epithelial *(endocrine/islet, A.c judgement)* |
| fibroblast | fibroblast_stromal |
| pancreatic stellate cell | fibroblast_stromal *(stellate = pancreatic CAF lineage)* |
| blood vessel endothelial cell | endothelial |
| myeloid leukocyte | myeloid |
| T cell | lymphoid |
| lymphocyte of B lineage | lymphoid |
| neural cell | other *(neural crest lineage, not stromal, per A.c)* |

Amendment 11 basis, verified in the run: the epithelial set is **identical** under
`celltype0` and `celltype1` (24,156 cells both ways); `celltype1` is a strict
refinement that separates myeloid from lymphoid, which `celltype0` collapses into
`hematopoietic cell`. The loader uses celltype0 for the epithelial/stromal/
endothelial branches and celltype1 to split haematopoietic.

---

## 4. Compartment method

All from `03_compartments.R`.

**Pseudobulk (A.d).** Raw counts summed within compartment, per gene:
`S[g, c] = rowSums(X[g, cells in c])` — **raw counts, not normalised**. Genes not
present in the atlas are `NA`, never 0. Detection is computed alongside as
`rowSums(X[g, cells] > 0) / n_cells[c]`.

**Purity grid.** `PI_GRID <- seq(0.30, 0.70, by = 0.01)` — **range 0.30–0.70, step
0.01, 41 grid points.**

**Fraction f(π) (A.e), as implemented in `sweep_f()`.** Per-cell intensity
`I = S / n_cells` per compartment; non-epithelial compartments are weighted by
their realised cell share `share = n_cells[non-epi] / sum(n_cells[non-epi])`; then
for each π:

```
w        = (epithelial = π,  non-epithelial = (1 − π) · share)
f(π)     = I[, epithelial] · w[epithelial]  /  (I %*% w)
```

**Evidence floor.** `EVIDENCE_MIN <- 20L` — a gene needs ≥ 20 summed counts in an
atlas to be evaluable; below that no fraction is reported.

**Dominance rule (`dominance_from_F`).** A gene is epithelial-dominant in an atlas
iff **f(π) > 0.50 at every one of the 41 grid points** (`all(r > 0.50)`), not at a
single point. Non-evaluable genes are `NA`, never `FALSE`.

**Bootstrap.** `B_RESAMPLES <- 2000L`, `SEED_BASE <- 20260731L`, resampled over
**patients** (the bootstrap unit above).

**k variants (`compute_k`), over the 152-gene inferential set only:**

| Variant | Definition | Value (95% CI) |
|---|---|---|
| **k** | genes dominant in **≥ 2** atlases | **46** (38–57) |
| k_all3 | dominant in **3** atlases **and** evaluable in all 3 | 24 (15–31) |
| k_evalall | dominant in ≥ 2 **and** evaluable in all 3 | 43 (35–55) |
| k_50 | dominant at the **midpoint π = 0.50 only** | 96 (86–103) |

Ordering asserted at runtime: `k_all3 ≤ k_evalall ≤ k`. CIs are bootstrap
percentile intervals (`output/k_estimates.csv`).

**Branch (Amendment 3).** `BRANCH_ABS_FLOOR <- 8L`; realised branch at k = 46:
**"epithelial subscore BUILT; subscore survival models run as SECONDARY."**

---

## 5. TCGA

Source: `data/cache/provenance.txt`, `output/clinical_cohort_summary.csv`,
`output/clinical_filter_flow.csv`, `05_clinical.R`.

| Item | Value |
|---|---|
| Retrieved | 2026-07-31 |
| TCGAbiolinks | 2.38.0 |
| GEOquery | 2.78.0 |
| **GDC data release** | **Data Release 45.0 — December 04, 2025** (`major 45, minor 0, release_date 2025-12-04`). Recorded at `04_lock_gene_list.R:135` and `PROJECT_SUMMARY.md:184` from the SummarizedExperiment objects; **re-queried live at `api.gdc.cancer.gov/status` 2026-08-06 and byte-identical**. Now written into `data/cache/provenance.txt`, whose placeholder at `01_download.R:183` had never been filled. |
| **Workflow** | **`STAR - Counts`** — the literal `workflow.type` argument at `01_download.R:68`. Assay column used downstream: **`tpm_unstrand`**, from the 4-assay, 60,660-row object (`unstranded, stranded_first, stranded_second, tpm_unstrand`). Now written into `data/cache/provenance.txt`. |

**Per cohort** (`clinical_cohort_summary.csv`):

| Cohort | Endpoint (primary) | Sensitivity | n | Events primary | Events sens | Stage missing | Redacted | Meta-eligible |
|---|---|---|---|---|---|---|---|---|
| TCGA-COAD | OS | PFI | 455 | 102 | 117 | 11 | 0 | yes |
| TCGA-READ | **PFI** | OS | 165 | 36 | 25 | 9 | 0 | yes |
| TCGA-STAD | OS | PFI | 406 | 157 | 130 | 24 | 0 | yes |
| TCGA-ESCA | OS | PFI | 184 | 77 | 87 | 23 | 0 | yes |
| TCGA-PAAD | OS | PFI | 178 | 93 | 104 | 3 | 0 | yes |
| TCGA-LIHC | OS | PFI | 370 | 130 | 178 | 24 | **1** | yes |
| TCGA-CHOL | OS | PFI | 35 | 18 | 20 | 0 | 0 | **no** |

**Endpoint basis.** Amendment 7 assigns the primary endpoint from the **TCGA-CDR
Table 3** recommendation per cohort: **PFI for READ**, OS for the other six.
Amendment 8 makes **CHOL descriptive-only** — its Table 3 explanation declares the
sample size too small for all endpoints — implemented as the flag
`meta_eligible = FALSE`, never by dropping the cohort. CDR source: Liu et al.,
`TCGA-CDR.xlsx` (PMC6066282, doi 10.1016/j.cell.2018.02.052).

**Administrative censoring.** `CENSOR_DAYS <- c(OS = 3650L, PFI = 1825L)` — 10
years OS, 5 years PFI, applied identically across cohorts.

**Filter flow, summed across the seven cohorts** (`clinical_filter_flow.csv`):

| Step | n |
|---|---|
| samples in | 1,998 |
| dropped, not primary tumour | 171 |
| dropped, duplicate aliquot | 23 (13 patients had >1 aliquot) |
| after sample filter | 1,804 |
| dropped, no CDR match | 3 |
| after CDR join | 1,801 |
| truncated at the primary horizon | 17 (1 event censored) |
| truncated at the sensitivity horizon | 80 |
| dropped, missing endpoint | 8 |
| **final** | **1,793** |

Aliquot tie-break: earliest **vial**, then **portion**, then **plate** (not a
full-barcode sort). The six meta-analysed cohorts contribute **1,755 patients and
595 events** at M1 (`meta_analysis.csv` + `survival_per_cohort.csv`); the
1,793 − 1,755 = 38 difference is CHOL (35) plus 3 patients lost to the
complete-case rule in §8.

---

## 6. Purity

Source: `06_purity.R`, `output/purity_summary.csv`, `data/cache/provenance.txt`.

| Item | Value |
|---|---|
| ESTIMATE package | **1.0.13**, loaded from the system library `/opt/homebrew/lib/R/4.6/site-library` (not installable from the configured conda channels) |
| Signature file | the package's own `extdata/SI_geneset.gmt` — 141 stromal, 141 immune genes |
| Genes scored | **9,911** in every TCGA cohort (asserted identical across cohorts) |
| **`platform` argument** | **`"illumina"`** |

**Why that platform argument.** Verified before use: the `platform` argument gates
**only** whether `estimateScore()` emits a `TumorPurity` column — the three scores
(Stromal, Immune, ESTIMATE) are computed identically on every platform. The
package computes `TumorPurity` for `affymetrix` only, and TCGA here is Illumina
RNA-seq, so purity is derived from the ESTIMATE score by the published
`cos(0.6049872018 + 0.0001467884 · S)` conversion, applied explicitly and labelled
`estimate_affymetrix_extrapolated` — i.e. an Affymetrix-calibrated transform used
off-platform, which is disclosed rather than hidden.

**CPE source.** Aran et al. 2015, *Nat Commun* **6**:8971, Supplementary Data 1
(`data/manual/aran_purity.xlsx`, md5 `c459e6a9…`, retrieved 2026-08-02): 9,364
sample rows, 21 TCGA cancer types.

**Coverage and the switch rule as applied** (plan §3.4: CPE is primary where
coverage is adequate, ESTIMATE where it is not):

| Cohort | Purity source | Calibration |
|---|---|---|
| TCGA-COAD | **CPE** | aran_consensus_cpe |
| TCGA-READ | **CPE** | aran_consensus_cpe |
| TCGA-LIHC | **CPE** | aran_consensus_cpe |
| TCGA-STAD | ESTIMATE | estimate_affymetrix_extrapolated |
| TCGA-ESCA | ESTIMATE | estimate_affymetrix_extrapolated |
| TCGA-PAAD | ESTIMATE | estimate_affymetrix_extrapolated |
| TCGA-CHOL | ESTIMATE | estimate_affymetrix_extrapolated |

The four ESTIMATE cohorts have **genuinely 0% CPE coverage**: STAD, ESCA, PAAD and
CHOL are **not in the Aran 2015 freeze at all** — verified by listing that table's
own cancer types, not inferred from a failed join. The switch therefore fires on
real absence, not on an acquisition gap. CPE–ESTIMATE agreement where both exist:
r = 0.70 (COAD), 0.66 (READ), 0.74 (LIHC).

---

## 7. Score construction

`score_cohort()` in `07_score.R`, in the order coded:

1. **Halt** if any scoring gene is absent from the expression matrix (the list was built to be present in every cohort).
2. Subset the expression matrix to the gene list.
3. Per-gene mean `mu` and SD `sdv` across the cohort's samples.
4. **Halt** if any gene has `sdv ≤ ZERO_VAR_TOL` — these must have been excluded **globally, never per cohort**.
5. **z-score each gene within cohort**: `z = (x − mu) / sdv`.
6. **Mean across the gene list** per sample: `s = colMeans(z)`.
7. **Rescale to unit SD within cohort**: `s = s / sd(s)`, asserted `sd(s) == 1`.

Input scale: **log2 TPM+1** (`expression_log2tpm`, which sums duplicate symbols on
the linear TPM scale *before* the log transform).

**Zero-variance rule.** `ZERO_VAR_TOL <- 1e-8`; a gene is zero-variance iff
`!is.finite(sd)` or `sd ≤ 1e-8`, evaluated on log2 TPM across each cohort's
primary-tumour samples. The registered wording is `sd > 1e-8` (`analysis_plan.md`
lines 662–663) — a **tolerance, not exact equality**. Excluded **globally** across
all cohorts, so the scoring set is identical everywhere.

---

## 8. Models

`08_survival.R`. All are Cox proportional-hazards, `survival::coxph`, Efron ties.

| Model | Formula as fitted | Parameters |
|---|---|---|
| M1 | `Surv(time, event) ~ score` | 1 |
| M2 | `Surv(time, event) ~ score + age + sex + stage_group` | 5 |
| M3 | `Surv(time, event) ~ score + age + sex + stage_group + purity` | 6 |
| M4 | `Surv(time, event) ~ score + age + sex + stage_group + purity + stromal_score` | 7 |

`score` is per **1 SD within cohort**. `stage_group` is I/II vs III/IV with an
explicit **missing** level (3 levels, 2 parameters). Parameter counts are **fixed
by registration**, not derived from the realised fit (`analysis_plan.md` line
1186); if `sex` is dropped in a cohort the count is reduced by exactly 1, and only
where a `sex` term was actually present.

**Sex rule.** Dropped in any cohort with < 10 patients of either sex
(`SEX_MIN_PER_LEVEL <- 10L`). Realised: **`use_sex = TRUE` in all seven cohorts**,
so no cohort lost the term.

**EPV rule.** `EPV_FIT <- 10`, `EPV_FLOOR <- 5`; EPV = events / parameters.
≥ 10 → `fit_and_pool`; 5–10 → `fit_pool_flag_LOO`; < 5 → `do_not_fit`.

**Realised band per cohort** (`output/clinical_epv.csv`):

| Cohort | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| COAD | fit_and_pool | fit_and_pool | fit_and_pool | fit_and_pool |
| **READ** | fit_and_pool | **fit_pool_flag_LOO** (7.2) | **fit_pool_flag_LOO** (6.0) | **fit_pool_flag_LOO** (5.14) |
| STAD | fit_and_pool | fit_and_pool | fit_and_pool | fit_and_pool |
| ESCA | fit_and_pool | fit_and_pool | fit_and_pool | fit_and_pool |
| PAAD | fit_and_pool | fit_and_pool | fit_and_pool | fit_and_pool |
| LIHC | fit_and_pool | fit_and_pool | fit_and_pool | fit_and_pool |
| **CHOL** | fit_and_pool | **do_not_fit** | **do_not_fit** | **do_not_fit** |

CHOL's M2–M4 are therefore **not fitted at all** (EPV 3.6 / 3.0 / 2.6), which is
why it is absent from every attenuation figure.

**Complete-case rule.** One analysis set per cohort, shared by M1–M4, requiring
non-missing `time, event, score, age, stage_group, purity, stromal_score` (+ `sex`
where used). Because the same set is used for all four models, the attenuation
difference is taken on identical patients.

---

## 9. Estimand and bootstrap

| Quantity | Definition |
|---|---|
| **attenuation_total** | β(M2) − β(M4) — **the primary result** |
| attenuation_purity | β(M2) − β(M3) |
| attenuation_stroma | β(M3) − β(M4) |
| prop_attenuated (secondary) | (β2 − β4) / β2, **reported only if** β2 ≠ 0 **and** sign(β2) = sign(β4); otherwise `NA` |

The registered direction is **positive** = adjustment reduces the association.

**Bootstrap.** `B_RESAMPLES <- 2000L`, `SEED_BASE <- 20260731L`, **paired** —
all four models refitted on each resample through the same code path, so the
difference is taken within resample. Interval type: **percentile**
(`quantile(v, c(0.025, 0.975))`, `boot_ci`); returns `NA` if fewer than 100
finite replicates. Every cohort achieved 2,000/2,000 usable resamples
(`boot_failed = 0` throughout).

---

## 10. Meta-analysis

`metafor::rma`, version **5.0.1**.

| Item | Value |
|---|---|
| Estimator | **REML** random effects |
| Primary interval | **Hartung–Knapp** (`test = "knha"`) |
| Also computed | Wald (`rma` REML default) and fixed-effect (`method = "FE"`), for cross-check |
| Heterogeneity reported | **τ², I², Q with its df and p** — all four, per pooled row |
| Multiplicity (B.n) | Benjamini–Hochberg **within model**, family = the six meta-analysed cohorts; CHOL carries `NA` and the label "not in family (descriptive, Amendment 8)" |

A recorded caveat for the text: at k = 3 with τ² = 0 (the exploratory colorectal
pool) the **Hartung–Knapp interval is narrower than Wald** — the known small-k
behaviour — so the Wald/FE bound is quoted as the conservative one there.

---

## 11. Null benchmark

`09_null.R`, per registered section B.m.

| Item | Value |
|---|---|
| N | **10,000 null signatures per cohort** |
| Set size | **152** (the registered panel size; all 152 panel genes are present in every cohort, so the "realised panel size" clause never bites) |
| Matching variables | **decile of mean log2 expression × decile of expression variance**, computed **per cohort** |
| Matching tolerance | set-level `MATCH_TOL <- 0.25` SD; paired per-gene `PAIRED_TOL <- 0.75` SD |
| Exclusion | null sets exclude **every panel gene** |
| Seed | `NULL_SEED <- BASE_SEED + 1000L`, i.e. `withr::with_seed(20260731 + 1000 + cohort_index)` |
| Model fitted | **M2** (as B.m specifies) |

**p-value formula** (`p_emp`), the add-one empirical p:

```
p = (1 + #{|null| >= |observed|}) / (1 + N)
```

with `N` the **requested** count (10,000), not the number of usable draws — so
failed draws cannot inflate significance. The two-sided form is primary; the
one-sided `v >= observed` is reported alongside as the directional companion.

**Three configurations and their roles** (`output/null_pvalues.csv`):

| Config | Role | N | p (two-sided, M1) |
|---|---|---|---|
| `registered_152 \| registered` | **primary_registered** | 10,000 | **0.3953** |
| `tested_140 \| registered` | size_sensitivity | 10,000 | 0.2652 |
| `tested_140 \| EXPLORATORY_POSTHOC` | **exploratory_posthoc** — additionally matched on stromal-score correlation, requested *after* the primary result was known | 1,000 | (labelled exploratory throughout) |

Roles were declared **before** the run, so no configuration could be selected
after seeing its p-value.

---

## 12. Validation

### GSE39582 (Amendment 16 SECONDARY) — `output/validation_gse39582_provenance.txt`

| Item | Value |
|---|---|
| Platform | **GPL570** (Affymetrix HG-U133 Plus 2.0) |
| Series matrix | `GSE39582_series_matrix.txt.gz`, md5 `9192d6561b8724caab6a5554c076eab6` |
| Samples | 585 → **19 non-tumour dropped** → 566 → 43 truncated at 10 y → **n = 561** |
| Events | **188** (OS) |
| Endpoint | OS from `os.event` / `os.delay..months.`, converted at 365.25/12 days per month; censored at 3,650 days, identical to 05. **Amendment 7's OS/PFI rule is defined on TCGA-CDR Table 3, which has no row for this cohort**; Amendment 16 fixes the estimand without naming an endpoint. |
| **Probe collapse** | drop **multi-symbol probes**, then take the **median** across a gene's probes. 54,675 probes → 21,656 symbols. |
| Cohort-independence | verified: max \|half-cohort − full-cohort\| = **0.000e+00**. The rejected max-mean-probe rule differs by up to 3.756 across 51 genes. |
| Probe rescue | **KRT17** has no single-symbol probe on GPL570 (both its probes are shared with JUP); rescued, and r(with rescue, without) = **0.9997** |
| Pooled with discovery | **FALSE** |

Note recorded in the provenance: `07`'s reader sums TPM before log2 and this
deposit is already log2 RMA, so that reader is **inapplicable, not merely unused**
— the collapse was re-implemented.

### FU-iCCA (Amendment 16 PRIMARY, executed under Amendment 18)

| Item | Value |
|---|---|
| Source | **Dong et al. 2022, *Cancer Cell*, supplementary Tables S1 (mmc2) and S5 (mmc6)** — **not** the NODE deposit OEP001105 (access unresolved) |
| md5 | mmc2 `f41b36e31e61c5ecc049ad636af35b1d`; mmc6 `6b331e884f97fd5cdd2423bd981bf3c8` |
| Sheets used | `S1A. Clinical info` (hdr row 2), `S1C. mRNA expression` (hdr row 2, log2 TPM+1), `S1D. proteins expression` (hdr row 3, median-norm log2), `S1E. 18,347 phosphosites` (hdr row 3, median-norm log2) |
| **n at each step** | **208** patients with both assays → **120** with a numeric STAT3:Y705 value → **114 joined** into the analysis (94 excluded for a missing Y705, **none imputed**) |
| Phosphosite | `STAT3:Y705`, a distinct row from `STAT3:S727` |
| Nomenclature | Amendment 19: **IL8 → CXCL8** for this cohort only, one symbol renamed, confined to panel scoring (ESTIMATE was run on the pre-rename symbols, whose reference set contains IL8) |

### TCGA RPPA (exploratory, post-hoc) — `data/validation/RPPA/provenance.txt`

| Item | Value |
|---|---|
| Source | **GDC Proteome Profiling**, 1,282 open-access files, retrieved via `api.gdc.cancer.gov` |
| Antibody, phospho | **AGID00388**, catalog **9131**, `peptide_target = "STAT3_pY705"` |
| Antibody, total | **AGID00185**, catalog **4904**, `peptide_target = "Stat3"` — a **different** antibody in the same files |
| Sample filter | primary tumour (`01`) only: 1,278 of 1,282 rows |
| Total STAT3 present | 1,259 of 1,278 (19 files lack it) |

**n per cohort** (RPPA patients / overlap with the analysis set):

| Cohort | RPPA | Overlap | Analysed |
|---|---|---|---|
| TCGA-COAD | 360 | **356** | yes |
| TCGA-READ | 131 | **127** | yes |
| TCGA-STAD | 357 | **327** | yes |
| TCGA-ESCA | 126 | **125** | yes |
| TCGA-PAAD | 120 | **113** | yes |
| TCGA-LIHC | 184 | **181** | yes |
| TCGA-CHOL | 0 | 0 | **no** — not meta-eligible (Amendment 8) **and** prior-knowledge disclosure: the author's prior ESMO Asia work used TCGA RPPA STAT3_pY705 in this cohort |

---

## 13. Software

**R version 4.5.3 (2026-03-11)** — recorded in `data/cache/provenance.txt` and
`data/panel/panel_provenance.txt`, and confirmed as the running version.

| Package | Version | Recorded in a committed file? |
|---|---|---|
| msigdbr | 26.1.0 | **yes** — `panel_provenance.txt` |
| TCGAbiolinks | 2.38.0 | **yes** — `data/cache/provenance.txt` |
| GEOquery | 2.78.0 | **yes** — `data/cache/provenance.txt` |
| estimate | 1.0.13 | **partly** — hand-written in `06_purity.R` (lines 35, 100) as a comment and a halt message, **not** captured at runtime and **not** in any run log |
| CMScaller | 2.0.1 | **yes** — `output/cms_provenance.txt` |
| survival | 3.8.9 | **no — read from the live library** |
| metafor | 5.0.1 | **no — read from the live library** |
| Matrix | 1.7.5 | **no — read from the live library** |
| rhdf5 | 2.54.1 | **no — read from the live library** |
| SummarizedExperiment | 1.40.0 | **no — read from the live library** |
| readxl | 1.5.0 | **no — read from the live library** |
| car | 3.1.5 | **no — read from the live library** |
| ggplot2 | 4.0.3 | **no — read from the live library** |
| cowplot | 1.2.0 | **no — read from the live library** |
| withr | 3.0.3 | **no — read from the live library** |
| matrixStats | 1.5.0 | **no — read from the live library** |

⚠ **The eleven "no" rows are the current library state, not a provenance record.**
`estimate` is a middle case and is marked "partly": v1.0.13 **is** written into
committed source (`06_purity.R` lines 35 and 100), but as a hand-authored comment
and halt-message string, not a version captured from the loaded package. It is
better evidence than the eleven live-library reads and weaker than `CMScaller`
2.0.1, which `output/cms_provenance.txt` records as a run-time output. Treat it as
an author's assertion, not a machine record.
No `sessionInfo()` is captured in any committed output. They are correct as of
2026-08-05 but are **not** evidence of what was loaded when each script ran, and a
methods section that cites them is citing today's library, not the run's. **The
one-line fix, if you want them defensible: add `sessionInfo()` to the tail of each
script's provenance block and re-run.** Also note `estimate` resolves from the
system library at `/opt/homebrew/lib/R/4.6/site-library` — an **R 4.6** path while
the analysis runs on **R 4.5.3** (it is pure R, so it loads, but the mismatch
should be stated).

---

## 14. Registration

**Correction (2026-08-06).** An earlier version of this document reported these as
`[NOT FOUND]`. That was wrong: they are committed in `README.md` lines 7–11 and
`HANDOFF.md` §5. My search missed them because it grepped `analysis_plan.md` and
`panel_definition.md` only. All four verified below.

| Item | Value | Verified in |
|---|---|---|
| **OSF registration** | **`tcvgb`** — https://osf.io/tcvgb/ — `registration: true`, public, not withdrawn | `README.md:7`, `HANDOFF.md:178,180` |
| **Registration date** | **2026-08-01**, timestamped **2026-08-01T04:26:17** | `README.md:7-8`, `HANDOFF.md:181` |
| **DOI** | **10.17605/OSF.IO/RKA4F** — this is the DOI of the mutable **parent project** `rka4f`, *not* of the registration | `README.md:9`, `HANDOFF.md:182-183` |
| **Registered commit** | **`468c7b4`** — confirmed present in git history (2026-07-31, *"Untrack installed R library and download byproducts"*) | `README.md:10`, `HANDOFF.md:140,186` |

⚠ **Two distinctions the methods section must not blur**, both stated in
`HANDOFF.md` §5:

1. **`tcvgb` is the registration; `rka4f` is not.** `rka4f` is the mutable parent
   project (`registration: false`) and it is the object carrying the DOI. Citing
   the DOI alone therefore points at a *mutable* object, not the frozen
   registration. Cite **both**: registration `tcvgb`, parent DOI
   10.17605/OSF.IO/RKA4F.
2. **The frozen snapshot does not match HEAD.** It archives the **pre-debug**
   `01_download.R` (6,147 B) and `02_panel.R` (4,386 B) rather than the working
   versions (8,224 B and 10,045 B) — the freeze took the 04:23–04:25 file state and
   corrected scripts went to the parent at 05:11, after it. It also contains no
   `commit.txt` (added to the parent at 05:12). The eight non-script files in the
   snapshot match HEAD exactly by SHA-256. `468c7b4` as the registered commit was
   verified by SHA-256 match against `commit.txt` on the parent.
| Analysis plan version | **1.5** (`analysis_plan.md` header) |
| Panel locked | 2026-07-31 |
| Registration status recorded | `analysis_plan.md` line 1117: *"Registration-blocking items: **none remaining**"* (as of 2026-08-01); `panel_definition.md` §8 states the intent to *"preregister after the panel is locked and after the compartment sweep is [run]"* |
| **Amendments** | **19 amendments**, numbered 1–20 with **no Amendment 17** — the numbering skips from 16 to 18. Verified by enumerating the headings. |

~~**The registration section cannot be written from this repository.**~~
**This paragraph is withdrawn.** It asserted that no OSF identifiers were recorded
anywhere; they are recorded in `README.md` and `HANDOFF.md`, and my search did not
cover those files. The registration exists, is frozen and is public.

---

## 15. Gap register — revised 2026-08-06

Of the seven gaps first reported, **four are now closed**, two stand, and one was
never a gap. Where a gap closed because the value was already committed somewhere I
had not searched, that is stated as **my search error**, not a new discovery.

| # | Fact | Status |
|---|---|---|
| 1 | **MSigDB release** | **CLOSED.** **`msigdb.2026.1`** — read from msigdbr 26.1.0's own hardcoded data URL, `https://zenodo.org/records/18968178/files/msigdb.2026.1.zip` (in `msigdbr:::check_cache`). See the note below on how this was obtained. |
| 2 | **ChEA3 version** | **NOT A GAP.** ChEA3 does not version its GMTs. The correct citation is the byte-size fingerprint plus the query date (2026-07-31), which §1 already gives. |
| 3 | **GDC data release** | **CLOSED.** **"Data Release 45.0 - December 04, 2025"**. Recorded at `04_lock_gene_list.R:135` and `PROJECT_SUMMARY.md:184`, read from the SummarizedExperiment objects. **Re-queried live at `api.gdc.cancer.gov/status` on 2026-08-06 and byte-identical** (`major 45, minor 0, release_date 2025-12-04`). Now written into `data/cache/provenance.txt`. My search missed it because I read only the unfilled placeholder at `01_download.R:183`. |
| 4 | **GDC workflow** | **CLOSED.** **`STAR - Counts`** — literal argument at `01_download.R:68` (`workflow.type = "STAR - Counts"`). Now in `data/cache/provenance.txt`. Same search error as gap 3. |
| 5 | **11 package versions** | **STANDS.** **There is no `renv.lock`** — not in the working tree, not tracked, and never in git history (0 commits touching it). A `renv/` directory exists but contains only `library-local/` (CMScaller, TCGAbiolinksGUI.data), no lockfile. Versions therefore **cannot be pinned to run time**; §13's eleven "no" rows are the library as of the reporting date. The manuscript should say so in those words. |
| 6 | **OSF ID, DOI, date, commit** | **CLOSED — my error, not a gap.** All four are committed in `README.md:7-11` and `HANDOFF.md` §5: registration **`tcvgb`**, registered **2026-08-01** (04:26:17), parent project **`rka4f`** DOI **10.17605/OSF.IO/RKA4F**, registered commit **`468c7b4`** (present in git history). I had grepped only `analysis_plan.md` and `panel_definition.md`. See §14 for the two distinctions that matter. |
| 7 | **Amendment 17** | **STANDS.** 19 amendments, numbered 1–16 and 18–20; **Amendment 17 was withdrawn**. A manuscript saying "20 amendments" is wrong. |

### How the MSigDB release was obtained, and its one caveat

msigdbr 26.1.0 does not ship the gene sets: it downloads them from Zenodo on first
use and caches them. `msigdbr_collections()` and `msigdbr()` therefore both fail
here (Zenodo is outside the sandbox allowlist; HTTP 403), and no local cache
exists. The release string was read instead from the **download URL hardcoded in
the package's own `check_cache()` function**, which names
`msigdb.2026.1.zip` — so this is the package's declared release, obtained without
any network call and without inference from the version number.

**Caveat for the methods section:** this establishes the release that msigdbr
26.1.0 *targets*. It is not an independent verification that the 87 criterion-A
genes were drawn from that release, because the objects that produced them are not
re-readable here. The 87-gene count itself is committed
(`data/panel/criterionB_provenance.txt`) and is the reportable figure.

### Two gaps that remain, stated for the manuscript

- **Package versions are as of the reporting date, not pinned at run time** (gap 5). No lockfile exists.
- **Amendment 17 does not exist**; the count is 19 (gap 7).

*Every value above was transcribed from the file named beside it. Nothing was
computed for this document, and nothing was filled in from memory.*
