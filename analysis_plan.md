# Analysis plan — compartment decomposition of STAT3 activity scores in GI cancers

**Version:** 1.3 (2026-08-01 — applies Amendment 8: EPV model-degradation rule in
B.j, per-model pooling in B.l, §3.3 rule resolved)
**Prior versions:** 1.2 (Amendment 7; B.i re-derived, COAD PFI→OS; `k_eval4`);
1.1 (Amendments 5–6; final gene list; §0.1/§0.2/§3.2 resolved); 1.0 (initial)
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

## 0. Two specification gaps — BOTH RESOLVED 2026-07-31

Both gaps below were open when this plan was first written. Both are now closed;
the original statements are retained for the record, each with its resolution.

### 0.1 — RESOLVED by Amendment 5

Amendment 3's "at least two tissue-matched atlases" is amended to "at least two of
the five compartment atlases, of any GI tissue". This is option (a) of the three
set out below. Option (b) — adding the Sci Data 2026 integrated atlas — was
considered and **rejected**: it is an integration of constituent GEO series, and if
GSE183904 or GSE178341 are among them, the two atlases would share cells, making
the replication requirement appear satisfied while supplying no independent
evidence. That is a sharper objection than the one raised in the original text
below, which noted only that Amendment 4 had made the atlas eligible.

Compensating reporting is required by Amendment 5 and specified in A.g: the
per-atlas dominance matrix, the per-tissue breakdown, and `k_all5`. The per-tissue
breakdown is what preserves the information the tissue-matched reading was trying
to protect.

Original statement, retained:

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
*(Now settled — see the resolution above. A.g no longer carries a pending marker.)*

### 0.2 — RESOLVED

The plan now asserts **both** counts for GSE178341: 129 GSM channels and 62 unique
patients, with a halt on either mismatch (A.a). The bootstrap unit remains
**patient** (A.f). The unit mismatch that prompted this note is therefore explicit
in the script rather than latent: the two assertions state plainly that they count
different things, and A.a records why channels outnumber patients roughly two to
one. A.a additionally halts if no patient identifier is present in the per-cell
metatable, since without it the specified bootstrap cannot be run at all.

Original statement, retained:

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

Expected: **129 GSM channels** with `specimen_type == "T"` (52 are `N`), resolving
to **62 unique patients** (28 MMRp + 34 MMRd, per Pelka et al. 2021). Both counts
are asserted; a mismatch in either halts the script.

```r
keep_gsm <- specimen_type == "T"

# Per-cell propagation must be verified on open, not assumed:
if (!any(c("specimen_type","channel","sample_type") %in% colnames(cell_meta))) {
  stop("HALT [GSE178341]: no tumour/normal column in per-cell metatable. ",
       "Join via GSM before filtering; do not proceed.")
}
if (!any(c("patient","patient_id","PatientTypeID","donor") %in% colnames(cell_meta))) {
  stop("HALT [GSE178341]: no patient identifier in per-cell metatable. ",
       "The A.f bootstrap unit is patient; do not proceed without it.")
}

assert_n(n_distinct(gsm_id[keep_gsm]),     129, "GSE178341", "specimen_type == T (GSM channels)")
assert_n(n_distinct(patient_id[keep_cell]), 62, "GSE178341", "specimen_type == T (unique patients)")
```

Asserting both counts is the point: 129 and 62 are counts of *different things*,
and only the second is the independent unit. The same specimen appears under
multiple `processing_type` values (`unosrted` [sic], `LiveMACS`), so channels
outnumber patients roughly two to one. If a future release of the series changes
either number, the script must stop rather than silently switch the effective
sample size or the bootstrap's independence assumption.

The published per-cell file is `GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz`.
If it carries no specimen column, the GSM-level designation must be joined onto
cells through the sample identifier before filtering — and the script halts rather
than guessing.

**Bootstrap unit remains patient** (A.f), not GSM and not cell.

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

### Strictness of the full-band rule, and `k_50` (Amendment 6)

Writing `A` for epithelial intensity per cell and `B` for the abundance-weighted
non-epithelial intensity, `f(π) = πA / (πA + (1−π)B)`, so dominance at π is
exactly `A/B > (1−π)/π`. The full-band requirement therefore binds entirely at the
lower boundary:

| π | Required `A/B` |
|---|---|
| 0.30 (band floor, binding) | **> 2.333** |
| 0.50 (typical purity) | > 1.000 |
| 0.70 (band ceiling) | > 0.429 |

The full-band bar is thus substantially stricter than dominance at typical purity:
it demands an epithelial cell carry more than twice the transcript mass of the
average non-epithelial cell. Per Amendment 6 this is the one criterion in the
design that leans toward the paper's thesis — a stricter bar yields a smaller *k*
and makes the descriptive-only branch more likely — and it is disclosed as such
rather than presented as neutral.

**`k_50` is therefore computed and reported alongside primary *k***, defined
identically (same panel, same two-atlas replication rule per Amendment 5, same
bootstrap) but evaluated at π = 0.50 only:

```r
dominance_50 <- f_at(0.50) > 0.50      # per gene per atlas
```

Both *k* and `k_50` are reported with 95% bootstrap intervals (A.f). **The branch
decision uses primary *k*.** If *k* and `k_50` fall in different branches of the
A.g table, that fact is stated explicitly in the paper — it is the cleanest single
indicator of how much the branch depends on the strictness of the band rule rather
than on the data.

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

*k* = number of panel genes that are epithelial-dominant (A.e) in **at least two
of the five compartment atlases, of any GI tissue** (Amendment 5).

### Required outputs (Amendment 5)

Three tables accompany *k* and are reported whatever branch obtains:

1. **Per-atlas dominance matrix.** Panel genes × five atlases, each cell holding
   the dominance indicator and `f[g](0.30)` (the binding value, per Amendment 6),
   with `NA` where the gene is absent from that atlas's annotation or falls below
   the 20-count evidence threshold. This is the primary evidence table for *k* and
   makes every gene's contribution auditable.
2. **Per-tissue breakdown.** Dominance counts by tissue — pancreatic (2 atlases),
   liver/biliary (1), colorectal (1), gastric (1) — so a reader can see whether *k*
   is carried disproportionately by any one tissue, which the pre-Amendment-5
   tissue-matched reading would have forced.
3. **Strict all-five count.** Number of panel genes epithelial-dominant in **all
   five** atlases, reported as `k_all5`. This is the conservative counterpart to
   *k*: it cannot be inflated by two agreeing atlases and is the number to quote
   when the claim is that a gene is epithelial regardless of tissue context.

4. **`k_eval4`: *k* restricted to well-evaluated genes.** *k* recomputed over only
   those panel genes **evaluable in at least four of the five atlases** — that is,
   with a non-`NA` dominance call in ≥ 4 atlases. Reported as an additional
   sensitivity alongside *k*, `k_50` and `k_all5`.

   This addresses a specific weakness of *k*: because absence encodes as `NA`
   rather than `FALSE`, a gene evaluable in only two atlases can reach *k* on those
   two alone, satisfying the replication requirement on the thinnest possible
   evidence. `k_eval4` asks what *k* would be if only genes with broad evaluability
   could contribute. Note it is not simply a subset relation on the count — a gene
   can be in *k* but not `k_eval4` (dominant in 2 of 2 evaluable) and the reverse
   cannot happen, so `k_eval4 ≤ k` always, and the gap between them measures how
   much of *k* rests on sparsely evaluated genes.

```r
# dominance: genes x atlases logical matrix, NA-aware
n_dom  <- rowSums(dominance, na.rm = TRUE)      # atlases calling dominance
n_eval <- rowSums(!is.na(dominance))            # atlases where the gene is evaluable

k          <- sum(n_dom >= 2)                              # primary (Amendment 5)
k_all5     <- sum(n_dom == 5 & n_eval == 5)                # strict
k_eval4    <- sum(n_dom >= 2 & n_eval >= 4)                # evaluability sensitivity
k_50       <- sum(rowSums(dominance_50, na.rm = TRUE) >= 2) # Amendment 6
per_tissue <- tapply(seq_len(ncol(dominance)), atlas_tissue,
                     function(j) sum(rowSums(dominance[, j, drop = FALSE],
                                             na.rm = TRUE) >= 1))

stopifnot(k_eval4 <= k, k_all5 <= k_eval4)     # ordering must hold by construction
```

The evaluability distribution itself is reported: a table of how many panel genes
are evaluable in 5, 4, 3, 2, 1 and 0 atlases. If most genes are evaluable in all
five, `k_eval4` ≈ *k* and the sensitivity is uninformative — which is itself worth
stating.

**The branch decision uses primary *k*** (A.g table). `k_50`, `k_all5` and
`k_eval4` are reported alongside with bootstrap intervals and do not override it.

A gene absent from an atlas contributes `NA`, not `FALSE`: absence is missing
evidence, not evidence against dominance. Consequently a gene present in only two
atlases can reach *k* on those two, and the dominance matrix's `NA` pattern is
reported so this is visible rather than hidden inside the count — with `k_eval4`
quantifying the consequence.

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

### The FINAL GENE LIST, and when it is fixed

The score is computed on a single **final gene list** used **identically in every
cohort**. It is not "the panel genes present in that cohort": a gene set that
varies by cohort makes the per-cohort hazard ratios estimates of different
quantities, and pooling them in B.l would be pooling non-comparable estimands.

**Sequencing, stated explicitly:**

1. The locked panel is 152 genes (`panel_definition.md`, locked at `ac9c5e0`).
2. The final list is the locked 152 **minus** the prespecification's exclusion
   criteria §3.1–§3.3:
   - **§3.1** — not present in the gene annotation of *every* TCGA cohort in the
     study. Evaluated from the seven `SummarizedExperiment` objects already on
     disk; requires no outcome data.
   - **§3.2** — not detected in the single-cell atlases at ≥ 1% of cells in at
     least one annotated compartment. **Evaluated from Part A output**, since Part
     A is what opens the atlases. Genes excluded here are *reported*, not silently
     dropped (prespecification §3).
   - **§3.3** — sex-chromosome genes.
3. The list is **locked and written to disk at the end of Part A**, to
   `data/panel/final_gene_list.csv`, with one row per locked-panel gene recording
   its exclusion status and the criterion that excluded it, plus a provenance
   header giving the Part A commit hash.
4. **No survival model is fitted before that file exists.** `08_survival.R` reads
   it and halts if it is absent:

```r
final <- readr::read_csv("data/panel/final_gene_list.csv", show_col_types = FALSE)
if (!nrow(final)) stop("HALT: final_gene_list.csv missing or empty. ",
                       "Part A must complete and lock the list before survival models.")
FINAL_GENES <- final$gene[is.na(final$excluded_reason)]
```

Because §3.2 depends on Part A, the final list is **not knowable now** — only its
derivation is. That is the intended order: the list is fixed by prespecified
criteria evaluated on compartment and annotation data, never on outcome data.

### Zero-variance genes: handled GLOBALLY, not per cohort

A gene with zero variance within one cohort but not others is **dropped from the
final list entirely, for all cohorts**, and recorded in
`final_gene_list.csv` with `excluded_reason = "zero_variance_in_<cohort>"`.

The alternative — dropping it only where it is flat — would reintroduce exactly the
cohort-varying gene set this section exists to prevent. Global removal keeps one
estimand across cohorts at the cost of a slightly smaller panel, which is the right
trade: a gene that is invariant in a cohort contributes nothing to that cohort's
score anyway (its z-score is undefined), so retaining it elsewhere buys
comparability-breaking heterogeneity for no signal. The rule is applied uniformly
and the count of genes lost this way is reported.

Zero variance is evaluated on log2 TPM across each cohort's primary-tumour samples,
using the tolerance `sd > 1e-8` rather than exact equality.

### Computation

```r
# se : SummarizedExperiment for one cohort, STAR-Counts
# FINAL_GENES : the single final gene list, identical for every cohort (above)
# 1. primary tumour samples only, one aliquot per patient
#    (sample_type == "Primary Tumor"; TCGA-CHOL and others carry normals)
# 2. expression matrix: tpm_unstrand, log2(x + 1)
# 3. subset to FINAL_GENES; halt if any is missing (the list was built to be present)
# 4. z-score EACH GENE WITHIN COHORT (mean 0, sd 1 across that cohort's patients)
# 5. score = mean across the final gene list
# 6. scale the resulting score to unit SD within cohort

score_cohort <- function(se, final_genes) {
  x <- log2(SummarizedExperiment::assay(se, "tpm_unstrand") + 1)
  missing <- setdiff(final_genes, rownames(x))
  if (length(missing)) {
    stop("HALT: final-list genes absent from this cohort: ",
         paste(missing, collapse = ", "),
         ". The final list is built to be present in every cohort; ",
         "do not proceed with a cohort-specific subset.")
  }
  x <- x[final_genes, , drop = FALSE]          # same genes, same order, every cohort
  z <- t(scale(t(x)))                          # per-gene z within cohort
  s <- colMeans(z)                             # mean across the final list
  as.numeric(scale(s))                         # per-SD scaling
}
```

Note the halt replaces the earlier `rownames(x) %in% panel_genes` intersection.
After the final list is locked, a missing gene is a defect in list construction, not
a condition to be silently absorbed — and `na.rm = TRUE` has been removed from the
mean for the same reason: with a properly built list there is nothing to drop.

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
- **Gene set identical across cohorts.** Superseding an earlier draft of this
  bullet, genes are *not* dropped per cohort. The final gene list is fixed once
  (above) and used unchanged everywhere; annotation absence and zero variance are
  handled globally at list construction, and a gene missing at score time is a halt,
  not a silent drop. The size of the final list and every exclusion is reported in
  `final_gene_list.csv`.

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

**Source, read 2026-07-31:** Liu et al. 2018, **Table 3** ("Recommended Use of
Survival Endpoints for Each Cancer Type"), obtained from the article
(PMC6066282; doi 10.1016/j.cell.2018.02.052). The `TCGA-CDR.xlsx` workbook's
`TCGA-CDR_Notes` sheet gives only the global recommendation — "we recommend the
use of PFI … and OS … Given the relatively short follow-up time, PFI is preferred
over OS" — and defers per-type detail to Table 3 of the paper, which is not
reproduced in the workbook. Table 3 was therefore read from the article. Only the
endpoint-usability recommendation table was consulted; no patient-level data was
opened, merged, or scored.

Table 3 legend: ✓ = recommended for use; × = not recommended; \* = caution, see
the explanation column; app. = approximate; acc. = accurate.

### The designation rule (Amendment 7)

> Primary endpoint = **OS wherever Table 3 marks OS usable without caution**.
> Where OS carries a caution mark, the primary endpoint is the CDR-preferred
> alternative for that cancer type. The non-primary endpoint is reported as a
> prespecified sensitivity analysis in every cohort.

This is the *only* input to the designation. Event counts do not enter it. No
assignment from v1.0 or v1.1 is carried forward; the table below is derived from
Table 3's marks alone.

### Table 3 marks, verbatim, and the endpoint the rule yields

OS marks read directly from Table 3's cells for these seven rows. `✓` = usable;
`✓ *` = usable with caution, see the explanation column.

| Cohort | Table 3 OS mark, verbatim | Table 3 PFI mark | OS cautioned? | **Rule yields** |
|---|---|---|---|---|
| TCGA-COAD | `✓` | `✓` | no | **OS** |
| TCGA-READ | `✓ *` | `✓` | **yes** — "need a longer follow-up for OS, DSS, and DFI; number of events for DFI is too small" | **PFI** (CDR-preferred alternative; `✓` unqualified) |
| TCGA-STAD | `✓` | `✓` | no | **OS** |
| TCGA-ESCA | `✓` | `✓` | no | **OS** |
| TCGA-PAAD | `✓` | `✓` | no | **OS** |
| TCGA-LIHC | `✓` | `✓` | no | **OS** |
| TCGA-CHOL | `✓` | `✓` | no | **OS**, but **descriptive only** — see below |

READ is the only cohort of the seven whose OS mark carries an asterisk. Verified
against the raw table markup, not only the parsed cells: the `✓ *` appears on
READ's OS cell, while COAD, LIHC and CHOL carry a bare `✓` on OS. (LIHC's and
CHOL's asterisks are on **DSS**, not OS, and DSS is not a candidate endpoint here.)

**CHOL.** The rule yields OS, and CHOL's OS mark is indeed uncautioned. But Table
3's explanation column for CHOL reads "sample size is too small for OS, DSS, DFI,
and PFI" — a caution against the *cohort*, not against a particular endpoint, and
one that no choice of endpoint can remedy. CHOL therefore keeps OS as its endpoint
under the rule and remains **descriptive only**, outside the meta-analysis as a
weighted stratum. This is a decision about cohort inclusion, not endpoint
designation, and the rule is applied to CHOL unchanged.

Cohort-level exclusion is now itself a rule rather than a judgement about this one
cohort: **Amendment 8** states that any cohort whose Table 3 explanation declares
the sample size too small for all endpoints is descriptive-only regardless of EPV.
CHOL is currently the only cohort meeting it.

### Change from v1.1

**COAD changes from PFI to OS.** In v1.1 COAD was assigned PFI, first (in v1.0) on
a nonexistent OS caution and then on the CDR's global preference plus event count.
Under Amendment 7 that reasoning is discarded entirely: COAD's OS is marked usable
without caution, so the rule yields OS. The row is changed, not silently replaced.

No other cohort's designation changes. READ remains PFI — now because its OS mark
is cautioned, which is the rule's own criterion, rather than because PFI had more
events. STAD, ESCA, PAAD, LIHC and CHOL remain OS.

The v1.1 justifications for STAD ("largest event count"), PAAD ("high event rate")
and LIHC ("PFI has more events, reported as sensitivity") are withdrawn as
justifications even where the endpoint is unchanged: they were event-count
reasoning, which Amendment 7 excludes. The endpoints now rest solely on Table 3's
marks.

### Event counts — descriptive context only

The counts below are the CDR's 2018 snapshot. They are reported for the reader's
orientation and played **no part** in the designation above.

| Cohort | N | OS events | PFI events |
|---|---|---|---|
| COAD | 459 | 102 | 123 |
| READ | 170 | 26 | 39 |
| STAD | 443 | 172 | 143 |
| ESCA | 185 | 77 | 87 |
| PAAD | 185 | 100 | 110 |
| LIHC | 377 | 132 | 185 |
| CHOL | 45 | 22 | 23 |

Note that in six of the seven cohorts (COAD, READ, ESCA, PAAD, LIHC, CHOL) PFI has
more events than OS; STAD is the sole exception. Under the v1.1 approach that fact
was doing work, and an event-count tiebreak applied consistently would have pushed
almost every cohort to PFI. Under Amendment 7 it does none, which is the point — an
event-count tiebreak is a choice the analyst makes after looking, and looking is
what the rule removes.

Realised event counts in the merged analysis set will differ from the above (the
expression cohorts are not identical to the CDR's clinical cohorts) and are
reported per cohort per §3.3.

One endpoint per cohort is designated primary in advance. The alternative endpoint
is reported as a prespecified sensitivity analysis in every cohort, never
substituted for the primary if the primary is null.

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

### Events per variable and the model-degradation rule (Amendment 8)

**Parameter counts.** Amendment 8's EPV denominator is the number of *estimated
parameters*, which exceeds the number of covariates because two covariates expand
into multiple dummies:

| Covariate | Levels | Estimated parameters |
|---|---|---|
| `score` | continuous | 1 |
| `age_at_diagnosis` | continuous | 1 |
| `sex` | 2 (reference + 1) | 1 |
| `stage_group` | **3** — I/II, III/IV, **`missing`** | **2** |
| `purity` | continuous | 1 |
| `stromal_score` | continuous | 1 |

`stage_group` carries an explicit `missing` level (B.j above retains
missing-stage patients rather than dropping them), so it contributes **two**
parameters, not one. Cox models estimate no intercept, so there is no additional
parameter for it.

| Model | Covariates | **Estimated parameters (EPV denominator)** |
|---|---|---|
| M1 | 1 | **1** |
| M2 | 4 | **5** |
| M3 | 5 | **6** |
| M4 | 6 | **7** |

Note this differs from the parenthetical in Amendment 8's reason clause, which
refers to "four covariates in M2 and six in M4" — correct as a covariate count, but
the EPV denominator specified in the same amendment is parameters, so M2 uses 5 and
M4 uses 7. Using parameters is the stricter reading and the one applied. It changes
no band assignment except CHOL's M2 (EPV 4.40 on parameters vs 5.50 on covariates),
and CHOL is excluded by the cohort-level rule regardless, so the distinction is
documented rather than consequential.

**Provisional EPV, from published CDR event counts.** Endpoint per Amendment 7;
events are the CDR's 2018 snapshot (B.i). **These are provisional pending the
realised complete-case n**, which will be smaller: the expression cohorts are not
identical to the CDR's clinical cohorts, and M1–M4 are fitted on a shared
complete-case set across all model-4 covariates (B.j), so missing purity, stage or
covariate data will reduce both n and events. Realised EPV is recomputed and
reported per cohort per model, and the rule is applied to the realised values, not
to these.

| Cohort | Endpoint | Events | M1 | M2 | M3 | M4 |
|---|---|---|---|---|---|---|
| COAD | OS | 102 | 102.0 | 20.4 | 17.0 | 14.6 |
| READ | **PFI** | 39 | 39.0 | **7.8** | **6.5** | **5.6** |
| STAD | OS | 172 | 172.0 | 34.4 | 28.7 | 24.6 |
| ESCA | OS | 77 | 77.0 | 15.4 | 12.8 | 11.0 |
| PAAD | OS | 100 | 100.0 | 20.0 | 16.7 | 14.3 |
| LIHC | OS | 132 | 132.0 | 26.4 | 22.0 | 18.9 |
| CHOL | OS | 22 | 22.0 | 4.4 | 3.7 | 3.1 |

Band assignment under Amendment 8:

| Cohort | M2 | M3 | M4 | Disposition |
|---|---|---|---|---|
| COAD, STAD, ESCA, PAAD, LIHC | fitted | fitted | fitted | all EPV ≥ 10; enter the meta-analysis as specified |
| **READ** | **low-EPV** | **low-EPV** | **low-EPV** | 5 ≤ EPV < 10 throughout: fitted, entered, **flagged**, and omitted in the leave-one-out sensitivity |
| **CHOL** | — | — | — | **Excluded at cohort level**, not on EPV: Table 3 states the sample size is too small for all endpoints. Descriptive only. (Its EPV would also fall below 5 for M2–M4.) |

Two observations worth recording:

- **READ sits in the low-EPV band on all three adjusted models but is never
  excluded**, because Amendment 7 designated PFI (39 events) rather than OS (26).
  Had OS been primary, READ's M3 and M4 EPV would be 4.33 and 3.71 — below 5 — and
  READ would have been dropped from the meta-analysis for both models. The two
  amendments interact, and the interaction favours retaining data; recorded so the
  dependency is visible rather than discovered later.
- **No cohort is expected to split across models** on these provisional numbers:
  every cohort falls in the same band for M2, M3 and M4. The per-model machinery in
  Amendment 8 and B.l is nonetheless implemented, because realised EPV can differ
  once complete-case attrition is applied — particularly for M3 and M4, which
  require purity (§3.4).

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

### EPV handling in the pooled estimates (Amendment 8)

**Flagging.** Every per-cohort row in the results table carries its realised EPV
and a band label (`fitted`, `low-EPV`, `not fitted`). Low-EPV cohorts are included
in the primary pooled estimate — Amendment 8 flags them, it does not exclude them —
so the primary meta-analysis is not conditioned on precision, which would bias it.

**Prespecified leave-one-out sensitivity.** A second pooled estimate is reported
with **every** low-EPV cohort omitted simultaneously (not one at a time — that is
the separate influence analysis above). On the provisional numbers this means the
pooled estimate excluding READ. Both estimates are reported side by side for each
model, and any material divergence is discussed rather than resolved in favour of
whichever is more congenial.

**Cohorts entering some models but not others.** Because the rule is applied per
model, the set of contributing cohorts can differ between M1–M4. Handling:

1. **Each model is pooled over exactly the cohorts that fitted it.** No imputation,
   no carrying a cohort's M2 estimate into the M4 pool.
2. **The contributing cohort set and its *n* are stated for every pooled estimate**,
   in the table itself, so no pooled figure appears without its denominator.
3. **The attenuation estimand (B.k) is pooled only over cohorts that fitted both
   M2 and M4**, since `attenuation_total = β₂ − β₄` is undefined for a cohort
   missing either. This set is stated explicitly and may be smaller than the M2
   pool.
4. **Where the contributing sets differ between M2 and M4, a matched-cohort
   comparison is additionally reported**: M2 and M4 both re-pooled over the
   intersection only. Without this, an apparent attenuation between the M2 and M4
   pooled estimates could be an artefact of different cohorts contributing to each,
   rather than of adjustment — which would corrupt the paper's primary claim. This
   is the single most important consequence of the per-model rule and is required,
   not optional.
5. τ² is estimated separately per model rather than shared across the sequence,
   since the contributing sets may differ.

On the provisional numbers no cohort splits across M2–M4 (B.j), so items 3–4 would
reduce to the same cohort set. They are specified because realised EPV after
complete-case attrition may differ, particularly for the purity-requiring models.

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

**Status as of 2026-08-01.** Registration-blocking items: **none remaining.** §3.2
is resolved (endpoints filled from CDR Table 3), §3.9 is resolved (Amendments 5
and 6, plus the dual-count assertion), and §3.3's *rule* is resolved (Amendment 8),
leaving only its realised numbers. Every other item is resolvable during Part A or
the Part B merge and is guarded either by a halt or by a prespecified rule, so that
an unresolved one stops the pipeline or is handled without discretion rather than
biasing the result.

| Item | Status | Blocks registration? | Resolved when |
|---|---|---|---|
| 3.1 atlas label strings | open | no | Part A, on atlas open; halt-on-unmapped |
| 3.2 endpoint recommendations | **RESOLVED** | — | done, CDR Table 3 |
| 3.3 event counts / EPV | **RULE RESOLVED** (Amendment 8); realised numbers open | no | Part B merge; rule leaves no discretion |
| 3.4 purity coverage | open | no | Part B merge; ≥80% switch rule prespecified |
| 3.5 GSE178341 per-cell propagation | open | no | Part A, on file open; two halts |
| 3.6 GSE183904 malignancy cutoff | open | no | Part A; may require an amendment |
| 3.7 Peng field names | open | no | Part A, on GSA deposit open |
| 3.8 whether *k*'s CI spans a branch | open by construction | no | after the sweep; disclosure prespecified |
| 3.9 §0.1 / §0.2 | **RESOLVED** | — | done, Amendments 5–6 and A.a |

## 3.1 Exact level-1 label strings for four of five atlases (A.c)

GSE125449's labels are exact — the per-cell table is on disk. GSE183904,
GSE178341, GSE155698 and Peng have their **mapping rule** fixed in A.c, but the
literal source strings must be read from each file. The published texts give
lineage vocabulary, not the annotation column's exact values, and Kumar's 34
lineage states and Pelka's 88 subsets are not enumerated verbatim in either paper's
main text. Mitigation: the halt-on-unmapped-label check in A.c makes an incomplete
map a loud failure rather than a silent miscount. Requires downloading the atlases,
which this session was instructed not to do.

## 3.2 Per-cohort endpoint recommendations (B.i) — RESOLVED 2026-07-31

`TCGA-CDR.xlsx` downloaded (GDC PanCanAtlas, 2,945,129 bytes). Its
`TCGA-CDR_Notes` sheet carries only the global recommendation and defers per-type
detail to **Table 3** of Liu et al. 2018, which is not reproduced in the workbook;
Table 3 was read from the article (PMC6066282). B.i is now filled with the actual
per-cohort recommendation, its event counts, and the CDR's caution flags, with the
designated primary endpoint and its basis for all seven cohorts.

Only the endpoint-usability recommendation table was consulted. No patient-level
data was merged with expression, no score was computed, nothing was fitted.

Two rounds of correction followed, both recorded in B.i rather than absorbed:

1. The v1.0 provisional table justified COAD→PFI by "the CDR notes OS underpowered
   in colon". Table 3 does not support that — COAD's OS is `✓` unqualified.
2. v1.1 corrected the fact but kept PFI on a substituted justification (global
   preference plus event count), while assigning LIHC to OS despite PFI having more
   events there — two contradictory tiebreaks. **Amendment 7** discards that
   approach and restates designation as one deterministic rule keyed only on Table
   3's caution marks. **COAD is now OS.** READ is the only cohort of the seven whose
   OS mark is cautioned, and so the only one taking the alternative endpoint.

## 3.3 Event counts, and whether any cohort is too small to model — RESOLVED by Amendment 8

The **rule** is now prespecified rather than left to a commitment to "report and
flag": EPV = events on the cohort's primary endpoint ÷ estimated parameters in the
model, applied **per model**, with EPV ≥ 10 fitted and pooled, 5 ≤ EPV < 10 fitted,
pooled and flagged with a leave-one-out sensitivity, and EPV < 5 not fitted for that
cohort and that model. Cohort-level exclusion is likewise a rule now: a cohort whose
Table 3 explanation states the sample size is too small for all endpoints is
descriptive-only regardless of EPV, which is what removes the earlier ad hoc
treatment of CHOL.

Parameter counts are fixed in B.j (M1 = 1, M2 = 5, M3 = 6, M4 = 7, with
`stage_group` contributing two parameters for its three levels). Provisional EPV
from the published CDR snapshot is tabulated there, and B.l specifies how differing
per-model cohort sets are pooled — including the required matched-cohort comparison
so that attenuation between the M2 and M4 pools cannot be an artefact of different
cohorts contributing to each.

**What remains** is only the realised numbers: complete-case *n* and events after
the CDR merge and purity join, which will be smaller than the published snapshot.
The rule is applied to the realised values. This is a Part B step, not a
registration blocker, since no discretion is left in how the numbers are used once
they exist.

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

## 3.9 §0.1 and §0.2 — RESOLVED 2026-07-31

Neither was a computation gap. §0.1 (Amendment 3's tissue-matching requirement) is
resolved by **Amendment 5**, which redefines replication as two of the five atlases
of any GI tissue and requires the per-atlas dominance matrix, per-tissue breakdown
and `k_all5` as compensating reporting. §0.2 (GSM channels vs patients) is resolved
by asserting both counts — 129 channels and 62 patients — with the bootstrap unit
fixed at patient and a halt if no patient identifier is present.

**Amendment 6** additionally closed a gap this section had not raised: the full-band
dominance rule is algebraically equivalent to `A/B > 2.33`, materially stricter than
dominance at typical purity (`A/B > 1`), and it leans toward the paper's thesis.
That is now disclosed, and `k_50` is reported alongside primary *k*.

