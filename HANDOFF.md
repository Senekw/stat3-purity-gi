# Project handoff — stat3-purity-gi

**Written 2026-08-01. Complete state transfer for continuing this project in a new session.**

Read this file first, then `panel_definition.md` (the locked prespecification, ten
amendments) and `analysis_plan.md` v1.5 (the preregistration document). Those two
are authoritative; this file is a map.

---

## 1. The research question

Does a STAT3 transcriptional activity score carry **tumour-intrinsic** prognostic
information in GI cancers, or does it primarily index **stromal and myeloid
content** in the bulk sample?

Testable hypothesis, committed in advance: after adjustment for tumour purity, the
epithelial-restricted component of the score retains an independent association
with overall survival. Report whichever way it falls.

### Why this is publishable, and the framing that makes it so

The gap is real but narrow. A Europe PMC sweep (2026-07-30, eight query families)
found:

- **STAT3 + tumour purity: 11 records, none on point.** Nobody has asked whether a
  STAT3 activity score is a purity artefact.
- **Prespecified epithelial/stromal decomposition of a named clinical score: 1
  record, unrelated.**

What is *not* novel — and must not be the headline — is "prognostic signal in bulk
tumour transcriptomes is substantially stromal." That was established in 2015:

| Prior work | Consequence for this paper |
|---|---|
| Isella 2015, *Nat Genet* (~849 cites); Calon 2015, *Nat Genet* (~495) | Stromal/CAF programs carry much of the prognostic signal in bulk CRC. Not your discovery. |
| Guinney 2015, *Nat Med* (~4,020 cites) | CMS4 mesenchymal = stroma-rich, poor prognosis. Reviewers will ask whether the stromal subscore is CMS4 relabelled → the B.o CMS orthogonality analysis exists to pre-empt this. |
| Venet 2011, *PLoS Comput Biol* (~470 cites) | Most random signatures associate with breast cancer outcome. The B.m null-signature procedure IS this method — cite it, don't claim it. |
| Aran 2015, *Nat Commun* (~962 cites) | Pan-cancer purity, source of CPE. Purity-as-confounder is generic prior art; the contribution must be STAT3-specific. |
| Stroma admixture in serous ovarian cancer, *CEBP* 2020 | The same design in a different disease. "Their paper, in GI" is the criticism to pre-empt. |
| STAT3 prognostic meta-analyses (*Oncotarget* 2016; *PLoS One* 2015, digestive pSTAT3) | pSTAT3 **IHC** is cell-resolved, so a bulk-stromal finding does **not** contradict them. Say so before a reviewer does. |

**The framing that works:** *a biomarker already used to select patients in a Phase
3 programme is, on inspection of its constituent genes, largely a readout of
non-tumour compartments — which is a concrete, testable explanation for why the
CO.23 pSTAT3 signal did not replicate.* This makes the paper a sequel to the
author's own napabucasin meta-analysis rather than a re-run of Isella 2015.

**Realistic venue:** specialty journal, not high-impact. The methodological core is
established; the contribution is a well-executed application to a clinically
anchored target. The napabucasin link is what lifts it above a routine signature
paper.

**Self-citation caution:** the author's ESMO Asia 2026 abstract reports this score
as a *validated prognostic biomarker*. This paper argues it is substantially a
stromal readout. That is legitimate self-correction and reviewers respect it, but
it must be explicit in the discussion, not left for someone else to notice.

### The author's related prior work

- **CCRC-D-26-00363** (submitted, *Clinical Colorectal Cancer*): napabucasin
  complete-enumeration meta-analysis across four registrational Phase 3 trials.
  Pooled OS-HR **1.028**. Tested non-replication of the pSTAT3 biomarker signal:
  CO.23's pSTAT3-positive OS-HR 0.41 (upper 95% bound 0.73) vs CanStem303C's 0.969;
  the two subgroup log-HRs differ significantly (**P = 0.006**). BRIGHTER's pSTAT3
  HR 1.32 and CanStem303C's 1.518 are **control-arm prognostic contrasts**, not
  treatment effects — do not read them as napabucasin effects. Trials: CO.23,
  BRIGHTER, CanStem303C, CanStem111P.
- **ESMO Asia 2026 abstract**: the six-gene pSTAT3 transcriptomic score in biliary
  tract cancer. Score = **mean log2(expr+1)** of **SOCS3, BCL2, MYC, MMP9, HGF,
  IL6**.
- **CCF conference** poster (photo, `IMG_3797.jpeg`).

---

## 2. Where the project stands

**Part A is written but NOT RUN.** `03_compartments.R` awaits review. Nothing
downstream has been computed.

| Stage | State |
|---|---|
| Panel construction (`02_panel.R`) | **DONE** — 152 genes locked 2026-07-31 |
| Data download (`01_download.R`) | **DONE** — all cohorts and atlases on disk |
| Part A compartment sweep (`03_compartments.R`) | **WRITTEN, NOT RUN — awaiting review** |
| Part B survival | **NOT STARTED.** Scripts `07_score.R`, `08_survival.R`, `09_meta.R` do not exist |

**Hard stop currently in force:** no compartment fractions, no purity sweep, no
final gene list locked, no score constructed, no survival model, no outcome data
merged with expression.

---

## 3. The governing working rule

Every session in this project has operated under a scope-discipline rule supplied
by the author. It is the single most important thing to carry forward:

> Do exactly what is asked and stop. Do not do adjacent work because it seems
> useful or because you are already close to it.

Specifics:
- Do not add analyses, plots, QC steps, or summary statistics not asked for.
- Do not get a head start on later scripts.
- Do not refactor beyond what makes the current step run.
- **Do not generate exploratory figures.** "In a preregistered study, looking at
  data I have not authorised looking at is a real cost, not a free bonus."
- If a step finishes early, stop. Do not fill the time.
- Out-of-scope observations go to `NOTES_FOR_REVIEW.md` rather than being acted
  on: "Surfacing an idea is helpful; executing it unasked is not."
- **Ask before proceeding whenever the next action is ambiguous.** "A question
  costs a minute; unrequested work costs the preregistration."
- **If you think a step is wrong, say so before executing** rather than working
  around it.

Amendments are supplied by the author as text to be added **verbatim**. Do not
edit, summarise, or append subsections to amendment text. Compute the consequences
and report them in chat instead.

---

## 4. Repository

`~/Documents/stat3-gi`, remote `https://github.com/Senekw/stat3-purity-gi.git`.

**HEAD = `2a2314234e33eceb89e06590a72514a36b7c22bc` (`2a23142`).**

**FOUR COMMITS ARE UNPUSHED** — `725c1f1`, `3eb1549`, `12e058e`, `2a23142`. The
sandbox has no GitHub credential; the author pushes from their own terminal.

### Commit history

| Hash | Date | Subject |
|---|---|---|
| `2a23142` | 08-01 | Write 03_compartments.R from analysis_plan.md v1.5; document removed file |
| `12e058e` | 08-01 | Remove unattributed 03_compartments.R; provenance not established |
| `3eb1549` | 08-01 | README: correct registration pointers and status to v1.5 / three atlases |
| `725c1f1` | 08-01 | Amendment 10: three-atlas compartment set; analysis_plan v1.5 |
| `db93583` | 07-31 | Create README.md |
| `468c7b4` | 07-31 | Untrack installed R library and download byproducts ← **registered commit** |
| `53ba089` | 07-31 | Untrack installed R library and download byproducts |
| `dfc9072` | 07-31 | Amendment 8: prespecified EPV model-degradation rule; plan v1.3 |
| `abfe6bc` | 07-31 | Amendment 7: endpoint designation as one deterministic rule; COAD PFI → OS |
| `4d7d40c` | 07-31 | Amendments 5-6; analysis_plan v1.1; endpoints filled from CDR Table 3 |
| `f089219` | 07-31 | Add analysis_plan.md; the document to be preregistered |
| `e53999e` | 07-31 | LOCK panel_definition.md: locked 2026-07-31 by Sean GP Lee at ac9c5e0 |
| `ac9c5e0` | 07-31 | Amendment 4: estimand changed to lineage-level epithelial fraction |
| `f4cc07b` | 07-31 | NOTES_FOR_REVIEW: record tracked-library and byproduct hygiene issue |
| `2fdc91b` | 07-31 | Amendments 2 and 3; atlas malignancy-annotation audit |
| `361d287` | 07-31 | Starting state: brief, prespecification, feasibility assessment, scripts, panel v1 |

### Tracked files

```
.gitignore  CLAUDE.md  README.md  HANDOFF.md
panel_definition.md      # LOCKED prespecification, 10 amendments
analysis_plan.md         # v1.5, the preregistration document
NOTES_FOR_REVIEW.md      # 10 sections of surfaced-but-not-acted-on items
feasibility_assessment.md
01_download.R  02_panel.R  03_compartments.R
data/panel/  panel_locked.csv, criterionB_evidence_table.csv,
             criterionB_provenance.txt, origin_six_evidence.csv,
             panel_provenance.txt, criterionA_genes.txt,
             6 ChEA3 .gmt files, trrust_rawdata.human.tsv
output/      01_download.log, atlas_malignant_annotation_audit.csv,
             atlas_tumour_designation_audit.csv, gastric_atlas_comparison.csv
```

**Known hygiene issue:** `renv/library-local/` (23 MB of installed R package
binaries) entered history in the author's own starting commit `361d287`. Untracked
going forward but still in history; removing it would require rewriting that
commit. Recorded as §8 of `NOTES_FOR_REVIEW.md`.

---

## 5. OSF registration — read this carefully

**The registration is `tcvgb`, NOT `rka4f`.**

- `https://osf.io/tcvgb/` — the actual registration. `registration: true`, public,
  not withdrawn, registered **2026-08-01T04:26:17**. A frozen snapshot.
- `https://osf.io/rka4f/` — the mutable parent **project**, DOI
  **10.17605/OSF.IO/RKA4F**. `registration: false`. This is *not* the registration,
  despite the DOI.

**`468c7b4` is the registered commit**, recorded in `commit.txt` on the parent
project. Verified by SHA-256 match (OSF's file download redirects to a denylisted
storage host, but the API publishes the checksum, and a 41-byte file is uniquely
identified as a 40-char hash + newline).

**Two defects in the registration snapshot:**

1. It archives the **pre-debug** `01_download.R` (6,147 B) and `02_panel.R`
   (4,386 B) — not the working HEAD versions (8,224 B and 10,045 B). The frozen
   snapshot took the 04:23–04:25 file state; corrected scripts were uploaded to the
   *parent project* at 05:11, after the freeze.
2. It contains no `commit.txt` (added to the parent at 05:12, after the freeze).

The eight non-script files in the snapshot match HEAD exactly by SHA-256.

**Also note:** `panel_definition.md` was *removed* from the parent project's file
list at 05:11. The registration snapshot still has it. Worth restoring on the
parent — it is the locked prespecification the analysis plan references throughout.

Network note: `api.osf.io`, `osf.io`, `files.osf.io` are allowlisted;
`storage.googleapis.com` (the final download redirect) is denylisted and cannot be
granted. Use the API's published checksums rather than downloading files.

---

## 6. The locked panel — 152 genes

Locked 2026-07-31 at commit `ac9c5e0`. File: `data/panel/panel_locked.csv`.

**Criterion A** (n=87): MSigDB `HALLMARK_IL6_JAK_STAT3_SIGNALING`.
**Criterion B** (n=75): human ChIP-seq (ChEA3) **AND** TRRUST v2 curation.
**Union = 152.** Routes: both 10, A_only 77, B_only 65.

### How criterion B was settled (Amendments 1 and 2)

As originally written ("at least two independent lines of evidence") criterion B
admitted **1,387 genes** — ChIP-seq occupancy is far broader than direct
transcriptional regulation, and an any-two rule counts multiple re-slices of the
same occupancy data as independent. **Amendment 1** replaced it with the
intersection rule. **Amendment 2** clarified B(i) as **human** ChIP-seq only.

Independence collapse performed before counting:
- Excluded as not binding evidence: `ARCHS4_Coexpression`, `GTEx_Coexpression`,
  `Enrichr_Queries`.
- Retained: `ENCODE_ChIP-seq` (5,480 genes), `ReMap_ChIP-seq` (1,453),
  `Literature_ChIP-seq` human (2,923).
- **Collapsed a truncated-PMID duplicate**: `STAT3_1855785_CHIPSEQ_MESC_MOUSE` is
  the same experiment as `STAT3_18555785_CHIPSEQ_MESC_MOUSE`.
- TRRUST v2 human: 9,396 edges, **142 STAT3 targets**.

**Amendment 2 was NOT blind** and says so: MMP9 and HGF — the two genes it excludes
— were the two most stroma/myeloid-dominant in the pilot. The restriction removes
the clearest illustrations of the paper's own thesis, so it is conservative with
respect to it. Human-only gives 152; including mouse would have given 171.

### The origin six

| Gene | In panel? | Route | Evidence |
|---|---|---|---|
| **SOCS3** | yes | both | ENCODE + ReMap + TRRUST |
| **MYC** | yes | B_only | ENCODE + ReMap + TRRUST |
| **IL6** | yes | both | ENCODE + Literature-human + TRRUST |
| **BCL2** | **no** | — | no ChIP-seq evidence in *any* group; fails under either counting rule |
| **MMP9** | **no** | — | mouse ChIP-seq only |
| **HGF** | **no** | — | mouse ChIP-seq only |

BCL2, MMP9, HGF are reported as a **labelled non-qualifying subset** per §4 of the
prespecification — never silently dropped.

---

## 7. The ten amendments

All in `panel_definition.md`. Each records whether it was made blind and its
direction of bias.

1. **Criterion B → intersection rule** (ChIP-seq AND TRRUST). Blind; panel size not
   computed at decision time.
2. **B(i) requires HUMAN ChIP-seq.** NOT blind — disclosed. Excludes MMP9, HGF.
3. **§5 thresholds restated as proportions** of the final panel: k ≥ 20% (=31) AND
   ≥ 8 → subscore built, survival models SECONDARY / 5–20% (8–30) → EXPLORATORY
   only / < 5% (≤ 7) → DESCRIPTIVE only. Written anticipating a 30–50 gene panel;
   at 152 the original `k ≥ 8` would have been a 5% bar, not the ~20% intended.
4. **Estimand changed** from malignant-epithelial to **lineage-level epithelial
   fraction on tumour samples only**. (Its claim that labels exist "in all five"
   atlases turned out to be wrong for two — see 9 and 10.)
5. **Replication requirement** loosened from "two tissue-matched atlases" to "two of
   the five, any GI tissue."
6. **Full-band strictness disclosed** + `k_50` added. At π=0.30 dominance requires
   A/B > 2.333; at π=0.50, > 1.0. The full-band rule is much stricter and leans
   toward the paper's thesis.
7. **Endpoint designation as one deterministic rule** — OS unless Table 3 marks OS
   with a caution, else PFI. COAD changed PFI → OS.
8. **EPV model-degradation rule.** ≥10 fit and pool; 5–<10 fit, pool, flag, with a
   leave-one-out sensitivity; <5 not fitted for that model. Cohort-level exclusion
   for any cohort whose Table 3 explanation says the sample is too small for all
   endpoints (currently CHOL only).
9. **GSE155698 (Steele, pancreatic) removed** — deposits no cell-type annotation.
   Five → four atlases. `k_all5`→`k_all4`, `k_eval4`→`k_eval3`. **NOT
   conservative** — disclosed.
10. **GSE183904 (Kumar, gastric) removed** — same defect. Four → **three** atlases.
    `k_all4`→`k_all3`; `k_eval3` and `k_evalall` collapse to `k_evalall`. Gastric
    and oesophageal coverage lost. **NOT conservative on two counts** — the
    replication bar rose 40% → 50% → 67% without ever being deliberately restated
    as a proportion. Also records GSE178341's 129/128 channel discrepancy.

**The annotation-availability finding is itself reportable:** Amendment 4's claim
that lineage labels were "available and comparably defined in all five" atlases was
taken from the source publications' *methods sections* rather than the deposits.
It was wrong for two of five. Annotation availability in public single-cell
deposits is materially worse than the corresponding papers imply.

---

## 8. Part A specification (analysis_plan.md v1.5)

### The three atlases

| Atlas | Tissue | Labels | Notes |
|---|---|---|---|
| **GSE125449** (Ma) | liver/biliary | `Type` in Set1/Set2 `samples.txt.gz` | 19 samples, all tumour (10 iCCA, 9 HCC). No filter — explicit no-op. Only atlas with CNV-based malignancy calls (inferCNV, score >80th pct AND corr >0.4) |
| **GSE178341** (Pelka) | colorectal | `clTopLevel` in cluster file | 129 GEO T channels / 128 with cells / **62 patients** |
| **Peng** (Besca reprocessing) | pancreatic | `obs/celltype0`, `obs/celltype1` | 57,423 cells, 35 patients (24 T + 11 N) |

**Removed:** GSE155698 (Amdt 9), GSE183904 (Amdt 10) — both deposit only count
matrices. Their tars are still on disk; do not use them.

### Peng provenance (A.b)

Zenodo **10.5281/zenodo.3969339**, file
`StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad`, md5
**`41fb7b9f27b7bb613ff979baaac5272f`** (verified on disk). Reprocessed from BIGD
**PRJCA001063** (= GSA CRA001160). **Annotations are Besca-derived, NOT the
authors' own.** Raw GSA rejected because it carries no annotation.

**The "CRC" in the filename is a workflow naming artifact — the data are
pancreatic.** Verified: vocabulary is pancreatic ductal / acinar / stellate /
enteroendocrine, no colorectal type at any level. The script halts if a colorectal
vocabulary appears.

### The GSE178341 counting trap (§0.2, and the cause of a real bug)

- **`PID`** = patient identifier → **62** (34 MMRd + 28 MMRp, matching Pelka).
- **`PatientTypeID`** = patient × specimen → **64**. Patients **C130** and **C171**
  each contributed two spatially distinct tumour specimens (`_TA`/`_TB`); 60 + 4 = 64.
- **`batchID`** (cluster file) = channel → 128 carry cells.
- GEO `specimen_type == "T"` → 129 GSMs.
- The missing channel is **`C144_T_1_1_12_c1_v2`** (`GSM5388094`,
  CD45pCD3nCD19nMACS fraction, zero cells, QC failure). C144 retains two other
  channels.

**The bootstrap unit is `PID`.** Using `PatientTypeID` would resample 64
pseudo-patients as independent, splitting C130 and C171 across draws and
understating the interval.

### A.c compartment map

Six compartments: `epithelial`, `fibroblast_stromal`, `myeloid`, `lymphoid`,
`endothelial`, `other`.

GSE125449: `Malignant cell`→epithelial, `HPC-like`→epithelial, `CAF`→stromal,
`TEC`→endothelial, `TAM`→myeloid, `T cell`/`B cell`→lymphoid,
`unclassified`→other.
GSE178341: `Epi`→epithelial, `Strom`→stromal, `Myeloid`/`Mast`→myeloid,
`TNKILC`/`Plasma`/`B`→lymphoid.
Peng: `celltype0` for most; `hematopoietic cell` collapses myeloid+lymphoid so
`celltype1` splits that branch only.

Judgement calls recorded: acinar and endocrine/islet → epithelial (inflates the
epithelial fraction, conservative under Amdt 4); enteric glial and Schwann → other
(neural crest, not fibroblast). **`TEC` moved** from stromal (pilot's
five-compartment scheme) to endothelial — pilot stromal fractions are therefore not
directly comparable to Part A output.

### A.d pseudobulk — the estimand rule

Sum **RAW UMI counts** per compartment, normalise **only afterwards**.
**FORBIDDEN:** CP10K, log, scaling, per-compartment means, or any library-size
correction *before* the `rowSums`. Normalising first discards the
transcripts-per-cell term and reweights every compartment to equal RNA content —
the exact error that made the HPA pilot measure the wrong quantity. Genes absent
from an atlas are `NA`, never 0.

### A.e purity sweep

`I[g,c] = S[g,c] / n_c`; `w_epi = π`, `w_c = (1−π)·n_c/Σn_c` for non-epithelial;
`f[g](π) = w_epi·I[g,epi] / Σ_c w_c·I[g,c]`.
**Grid: π 0.30→0.70 step 0.01 (41 points).** Dominance iff `f(π) > 0.50` at
*every* grid point. f is monotone increasing in π, so this binds at π=0.30 — but
monotonicity is **asserted numerically**, not assumed.
Evidence threshold: genes with <20 summed counts in an atlas are
insufficient-evidence, not given a fraction.

### A.f bootstrap

Unit **patient**. B = **2000**, `set.seed(20260731)`, per-atlas streams via
`withr::with_seed(20260731 + atlas_index)`. The whole A.d→A.e chain is recomputed
inside each resample, not just the final ratio. 95% percentile CIs.

### A.g — k and variants

- **k** = panel genes epithelial-dominant in **≥2 of 3** atlases (any GI tissue)
- **k_all3** = dominant in all three
- **k_evalall** = k restricted to genes evaluable in all three
- **k_50** = as k, evaluated at π=0.50 only (Amdt 6)

Ordering `k_all3 ≤ k_evalall ≤ k` **verified exhaustively over all 3³ = 27
patterns: zero violations, 4 distinct outcome triples** (so the sensitivities can
actually differ from the primary).

**Amendment 3 branch, on primary k against 152 genes** — the amendment's own
wording, which matters because the middle band is weaker than "with caveats":

| k | Branch |
|---|---|
| **≥ 31** (20% of 152, and ≥ 8) | epithelial subscore **BUILT**; subscore survival models run as **SECONDARY** |
| **8–30** (5% ≤ k < 20%) | subscore reported as **EXPLORATORY only** |
| **≤ 7** (< 5%) | decomposition **DESCRIPTIVE only** |

**Expectation:** the pilot found only MYC majority-epithelial of the original six,
so a low k and the descriptive-only branch is the likely outcome. Amendments 9 and
10 both push k downward and both disclose this as non-conservative.

---

## 9. Part B specification — NOT STARTED

Scripts `07_score.R`, `08_survival.R`, `09_meta.R` do not exist. Full spec in
`analysis_plan.md` §B.h–B.o.

**B.h score:** primary tumour only, one aliquot per patient; `tpm_unstrand`,
log2(x+1); subset to the final gene list (halt if any missing); **z-score each gene
within cohort**; mean across the list; scale to unit SD within cohort. The **final
gene list is identical for every cohort** — not re-derived per cohort.

**B.i endpoints** (Amendment 7, from TCGA-CDR Table 3):

| Cohort | Endpoint | Why |
|---|---|---|
| COAD, STAD, ESCA, PAAD, LIHC | **OS** | OS mark uncautioned |
| **READ** | **PFI** | OS carries `✓*` |
| CHOL | OS, **descriptive only** | Table 3: sample too small for all endpoints |

Verified against **raw table markup**, not parsed cells — LIHC's and CHOL's
asterisks are on **DSS**, not OS, and a parser that collapses columns would flip
both cohorts to PFI on a misread.

**B.j models & EPV** (Amendment 8). Parameters, not covariates:
M1=1, M2=**5**, M3=**6**, M4=**7** (`stage_group` has three levels — I/II, III/IV,
`missing` — so contributes **two** parameters; Cox estimates no intercept).

| Cohort | Endpt | Events | M2 | M3 | M4 | Disposition |
|---|---|---|---|---|---|---|
| COAD | OS | 102 | 20.4 | 17.0 | 14.6 | fitted |
| **READ** | PFI | 39 | **7.8** | **6.5** | **5.6** | **low-EPV, flagged** |
| STAD | OS | 172 | 34.4 | 28.7 | 24.6 | fitted |
| ESCA | OS | 77 | 15.4 | 12.8 | 11.0 | fitted |
| PAAD | OS | 100 | 20.0 | 16.7 | 14.3 | fitted |
| LIHC | OS | 132 | 26.4 | 22.0 | 18.9 | fitted |
| CHOL | OS | 22 | 4.4 | 3.7 | 3.1 | **excluded at cohort level** |

Provisional — realised complete-case n will be smaller.

**Amendments 7 and 8 interact:** READ survives only because Amdt 7 designated PFI
(39 events) over OS (26). On OS its M3/M4 EPV would be 4.33 and 3.71 — below 5 —
and READ would drop out of both models.

**B.l pooling — the load-bearing requirement:** because EPV is applied per model,
contributing cohort sets can differ between M2 and M4. If they do, a **matched-cohort
comparison re-pooling M2 and M4 over the intersection is REQUIRED, not optional** —
otherwise apparent attenuation between the two pools could be an artefact of
different cohorts rather than of adjustment, corrupting the paper's primary claim.

**B.k estimand:** `attenuation_total = β₂ − β₄`, pooled only over cohorts fitting
both.

**B.m null signatures:** Venet 2011 method. **B.n:** multiplicity across seven
cohorts. **B.o:** CMS orthogonality in CRC — pre-empts the "is this just CMS4"
objection.

---

## 10. Data on disk

Not in git (`.gitignore`d). ~4 GB.

```
data/raw/
  GSE178341_crc10x_full_c295v4_submit.h5                    1.1G
  GSE178341_..._cluster.csv.gz                              2.8M   ← clTopLevel labels
  GSE178341_..._metatables.csv.gz                           2.4M   ← PID, SPECIMEN_TYPE
  StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad              1.6G   ← Peng
  GSE155698_RAW.tar                                         881M   ← REMOVED (Amdt 9), no labels
  GSE183904_RAW.tar                                         329M   ← REMOVED (Amdt 10), no labels
  TCGA-{COAD,READ,STAD,ESCA,PAAD,LIHC,CHOL}/                       ← GDC raw
data/tcga/   TCGA-{...}_se.rds        7 cohorts, 60,660 genes, 6 assays, 1,998 samples total
data/geo/    GSE125449/{Set1,Set2}_{barcodes,genes,matrix,samples}
             GSE39582_eset.rds  (CRC validation, ~560 w/ survival)
             GSE66229_eset.rds  (ACRG gastric, ~300)
             GPL570.soft.gz
data/manual/ TCGA-CDR.xlsx        ← Liu et al. Cell 2018, curated survival endpoints
```

**Missing, needed before Part A runs:** `data/raw/GSE178341_geo_tumour_channels.txt`
— the committed list of 129 GEO tumour channel IDs. The script reads a committed
artefact rather than re-fetching from GEO at analysis time, and halts without it.

**Also relevant, not downloaded:** Dong et al. iCCA proteogenomics (PMID
**34971568**, *Cancer Cell* 2022, n=262) — BioSino NODE **OEP001105** + Table S1.
Paywalled, no PMC copy. Carries a **measured STAT3 pY705 phosphosite**, the
out-of-sample validation the ESMO score never got. Caveat: measured in **bulk
tissue**, so it carries the same compartment ambiguity as the mRNA score.

---

## 11. Environment and infrastructure gotchas

**Conda env `stat3-gi`**, R 4.5.3. Installed: tidyverse, data.table, survival,
meta/metafor, jsonlite, BiocManager, TCGAbiolinks, GEOquery, msigdbr, R.utils,
**rhdf5**, **withr**, Matrix.

Hard-won fixes, all of which cost real time:

1. **R's default URL method does not honour the HTTP proxy.** TCGAbiolinks
   misreports the resulting failures as GDC being down while `curl` reaches it
   fine. Fix: `options(download.file.method="libcurl", url.method="libcurl")` and
   `timeout=3600`. This is at the top of `01_download.R`.
2. **The GDC API fails intermittently through the proxy** (~1 call in 5).
   `01_download.R` has retry-with-backoff, hardened specifically for `GDCprepare`'s
   clinical-enrichment step which the package's own retry cannot recover. TCGA-COAD
   succeeded only on retry 7.
3. **bioconda's `TCGAbiolinksGUI.data` ships a 0-byte tarball** — its post-link
   download never runs. The conda R library is read-only, so install to a
   project-local library from bioconductor.org directly.
4. **`git init` is blocked in the sandbox** — `.git/config` is a protected path. The
   author ran it manually.
5. **No GitHub credential in the sandbox** — the author pushes.
6. **The python kernel cannot read host-granted paths if it started before the
   grant.** Use `bash` (spawns fresh) or the env's interpreter directly.
7. Allowlisted domains added during this project: `bioconductor.org`,
   `maayanlab.cloud`, `www.grnpedia.org`, `api.osf.io`, `osf.io`, `files.osf.io`.

---

## 12. `NOTES_FOR_REVIEW.md` — ten open observations

Surfaced, not acted on. §10 is the most important.

1. `git init` blocked in sandbox (resolved by the author).
2. Used the corrected `feasibility_assessment.md`, not the pasted copy.
3. Criterion B yielded ~1,400 genes → became Amendment 1.
4. Three of the original six fail criterion B → BCL2, MMP9, HGF.
5. **`01_download.R` §2 — the ACRG gastric accession is still ambiguous.**
   GSE66229 is fetched; GSE62254 is the commonly cited alternative. Unresolved.
6. Sci Data 2026 integrated GI atlas has no malignant-cell annotation. Rejected as
   a gastric substitute in Amdt 10 (post-registration dataset selection is a
   researcher degree of freedom).
7. **Malignancy is inferred by three different methods across atlases** — Ma:
   inferCNV (>80th pct AND corr >0.4); Pelka: a *trained classifier* cross-checked
   against inferCNV; Steele: none at all.
8. `renv/library-local/` in git history (23 MB).
9. Three hazards in tumour-vs-normal designation.
10. **The unattributed `03_compartments.R`** — see §13.

---

## 13. The unattributed script incident — record-integrity item

A file named `03_compartments.R` (6,721 B, dated 2026-08-01 00:41) appeared in the
working tree, untracked. **It was not written by the assistant** — the two
preceding turns both ended with an explicit statement that it had not been written.
Its provenance could not be established. It was swept into `725c1f1` by `git add
-A` and removed in `12e058e`.

Three defects, all confirmed by direct test:

1. **Line 79, two errors.** `sum(keep178) != 129L` counted *cells* (257,251), not
   channels — the comparison could never pass. And `PatientTypeID` gave 64, not 62.
   `PID` was never referenced.
2. **Transposed tissue labels** — GSE125449 tagged `gastric`, GSE183904 tagged
   `liver_biliary`. Hardcoded, so nothing could catch it.
3. **`map_compartment` mapped `Malignant cell` → `lymphoid`.** `"T cell"` matches
   case-insensitively inside "Malignan**T CELL**" and the lymphoid pass overwrote
   the epithelial one. **Every malignant cell in the biliary atlas would have been
   counted as a lymphocyte, inverting the study's own estimand** — and no assertion
   could catch it, because `lymphoid` is a valid compartment. The same map returned
   `NA` for 3 of 7 real colorectal labels and 2 of 5 real Peng labels.

Removed rather than repaired: a preregistered analysis cannot rest on code of
unknown origin.

**Design lesson carried into the replacement:** defects 2 and 3 share a root cause —
compartment and tissue assignments were *asserted by the author* rather than
*derived and checked against the data*.

---

## 14. The current `03_compartments.R` — written, awaiting review

39 KB, 101 top-level expressions, **44 halts, 15 counted assertions**. Parses
clean. Written from v1.5 from scratch, no code reused.

- **`atlases`** declares tissue, but `verify_tissue()` checks each claim against the
  atlas's own gene universe and label vocabulary via `TISSUE_WITNESS`, and requires
  the claimed tissue to score **better than every other declared tissue**. Tested: a
  transposition halts.
- **`map_labels()`** is **exact-match against enumerated per-atlas vocabularies**,
  not regex. No fallback, no default. Unmapped → halt. A regression guard asserts
  `Malignant cell` → `epithelial` before any file opens.
- **GSE178341** asserts all three counts plus the **identity** of the absent
  channel, and asserts `PatientTypeID == 64` separately so the distinction is
  logged rather than rediscovered.
- **Peng** asserts md5 and has an explicit anti-colorectal vocabulary guard.
- **`verify_k_ordering(3)`** runs before any data is opened, with a collapse guard.

**Two open questions for the author:**

1. **Peng's `hematopoietic cell`** collapses myeloid+lymphoid at `celltype0`, so
   `celltype1` is used for that branch only. The epithelial/non-epithelial split
   that f(π) depends on is identical either way, so the estimand doesn't move — but
   it refines A.c's "level-1 label" wording and may want an amendment.
2. **`GSE178341_geo_tumour_channels.txt` does not exist yet** (see §10).

---

## 15. Immediate next steps

1. **Push the four unpushed commits** (`725c1f1`, `3eb1549`, `12e058e`, `2a23142`).
2. **Review `03_compartments.R`** — specifically `atlases`, `map_labels`, and the
   assertion set.
3. **Decide the Peng `celltype1` question** — amendment or not.
4. **Create `data/raw/GSE178341_geo_tumour_channels.txt`** (129 channel IDs, parsed
   from GEO `specimen_type == "T"` titles).
5. **Consider restoring `panel_definition.md`** to the OSF parent project.
6. **Then run Part A** — and if any assertion fires, stop and report rather than
   adjusting the assertion, the filter, or the expected count.
7. Part B remains unwritten and unstarted.
