# Notes for review

Observations surfaced but NOT acted on, per the scope discipline rule.

## 1. `git init` is blocked in this sandbox (step 1, partially incomplete)

`git init` inside the granted host path fails with `Operation not permitted` on
`.git`. Ordinary directories create fine; the sandbox specifically refuses to
create git config paths. Everything else in step 1 is done (directory structure,
docs, scripts, `.gitignore` written).

To finish, run once in a terminal:

```
cd ~/Documents/stat3-gi
git init
git add -A
git commit -m "Starting state: brief, prespecification, feasibility assessment, untested scripts"
```

Nothing downstream depends on this.

## 2. `feasibility_assessment.md` — used the corrected version

The copy pasted into this session still says GSE125449 comprises "11 iCCA and 8
HCC" samples. The true counts from the GEO sample characteristics are **10 iCCA
and 9 HCC** (total 19 and 9,946 cells are both correct). The repo copy is the
corrected v3. If you have another copy elsewhere, it carries the stale numbers.

## 3. Criterion B yields ~1,400 genes, not 30-50 (decision needed — see chat)

`panel_definition.md` §5 anticipates a panel where *k* ≥ 8 is the top branch,
implying a panel of tens of genes. Criterion B as literally written (≥2
independent sources, ChEA3 + TRRUST at default thresholds) admits **1,387
genes**; the union with criterion A is **1,448**. ChIP-seq occupancy is simply
much broader than direct functional regulation. Not resolved unilaterally
because tightening the rule is an amendment to a document that is meant to be
binding.

## 4. Three of the original six fail criterion B under a human-only rule

BCL2, MMP9 and HGF are supported by TRRUST only among human sources — one source,
not two. MMP9 and HGF each additionally appear in mouse ChIP-seq, so an
include-mouse counting rule would pass those two with 2 sources. **BCL2 has no
ChIP-seq evidence in any group, human or mouse** — TRRUST curation alone — so it
fails under either counting scheme. §4 of the prespecification says a failing origin gene "is
dropped, and that is reported", so this is reported rather than worked around.
Note the interaction with the feasibility assessment: MMP9 and HGF are the two
most stroma/myeloid-dominant genes in the pilot, so dropping them would remove
the clearest examples of the paper's own thesis. That is an argument for
reporting them as a labelled subset, not for bending the criteria.

## 5. `01_download.R` section 2 — ACRG accession still ambiguous

The script's own comment flags this. GSE66229 is fetched; the script notes
GSE62254 as the commonly cited alternative. Not investigated — out of scope for
this session, and the gastric-atlas question in step 3 concerns single-cell
data, not this bulk cohort.

## 10. An unattributed `03_compartments.R` appeared in the tree and was removed

**This is a record-integrity item, not a code-quality one.** It is documented
because a preregistered study's audit trail has to account for every file that
influenced, or could have influenced, a registered analysis.

**When it appeared.** A file named `03_compartments.R` (6,721 bytes, dated
2026-08-01 00:41) was found in the working tree, untracked by git. It was **not
written by the assistant session**: the two turns immediately preceding its
discovery both ended with an explicit statement that the script had not been
written, pending (a) resolution of the OSF registration precondition and (b) the
GSE155698 missing-annotation halt. Its author and origin could not be
established. It was swept into commit `725c1f1` as an untracked add — an error on
the assistant's part, since `git add -A` picked it up — and removed in `12e058e`.

**What it contained.** A fail-closed Part A implementation written against
analysis-plan **v1.4** (the superseded four-atlas state): a `halt()` helper, a
four-row `atlases` frame, a locked-panel loader asserting 152 genes, an
exhaustive `3^4` k-ordering check, a regex `map_compartment`, per-atlas
assertions for GSE125449 / GSE178341 / Peng, a deliberate terminal halt at
GSE183904, and — below that halt, unreachable — a `pseudobulk_raw` stub and a
41-point purity grid.

**The three defects**, all confirmed by direct test before removal:

1. **Line 79 — two errors on one line.** `sum(keep178) != 129L` counts *cells*,
   not channels: `keep178` is a logical over the 370,115 rows of the per-cell
   metatable, so `sum()` returns **257,251**, and the comparison could never have
   passed. On the same line, `length(unique(m178$PatientTypeID[keep178])) != 62L`
   uses `PatientTypeID`, which is patient × specimen. Over tumour cells it takes
   **64** distinct values — 60 patients with a single `_T` specimen, plus C130 and
   C171 which each contributed two spatially distinct tumour specimens
   (`_TA`/`_TB`); 60 + 4 = 64. The correct field is `PID` (62), which the file
   never references. This produced a spurious "64 patients" halt. Had it
   propagated to A.f, the bootstrap would have resampled 64 pseudo-patients as
   independent units, splitting C130 and C171 across draws and understating the
   interval. It did not propagate — the file contains no bootstrap code.

2. **Transposed tissue labels.** In `atlases`, GSE125449 (liver/biliary) is
   tagged `"gastric"` and GSE183904 (gastric) is tagged `"liver_biliary"`. Every
   per-tissue breakdown in A.g would have carried the wrong tissue. Because the
   labels were hardcoded rather than derived from the atlas source, no assertion
   in the script could have detected it.

3. **`map_compartment` misassigns malignant cells to the lymphoid
   compartment.** The function applies six `grepl` passes in sequence, each
   overwriting the last, against unanchored case-insensitive substrings. The
   lymphoid pattern includes `"T cell"`, which matches inside **"MalignanT
   CELL"** — and runs *after* the epithelial pass. Applied to GSE125449's real
   `Type` vocabulary, the label **`Malignant cell` maps to `lymphoid`**. Every
   malignant cell in the liver/biliary atlas — the atlas the entire biliary arm
   of this study rests on — would have been counted as a lymphocyte, inverting
   the quantity the study exists to measure. This is the most serious of the
   three and would not have tripped any halt, because the label *did* map to a
   valid compartment.

   The same map is also a poor fit for its targets: it returns `NA` for three of
   seven real GSE178341 `clTopLevel` labels (`B`, `Epi`, `Strom`) and for two of
   five real Peng `celltype0` labels (`hematopoietic cell`, `neural cell`).

   *Correction to an earlier report:* the assistant initially speculated that the
   unanchored `"pit"` token would also match `"pericyte"`. Tested — it does not.
   The `Malignant cell` collision is the real defect, and a worse one.

**Removed, not repaired.** Repairing code of unknown provenance would leave a
preregistered analysis resting on a file nobody can account for; the three
defects are also evidence that it was not written against v1.5. The replacement
was written from `analysis_plan.md` v1.5 from scratch, reusing none of it. The
deleted file remains recoverable from git history at `725c1f1` for audit.

**Design consequence carried into the replacement.** Defects 2 and 3 share a
root cause: compartment and tissue assignments were *asserted by the author*
rather than *derived and checked against the data*. The replacement therefore
derives tissue from each atlas's own content, maps labels by exact match against
an enumerated per-atlas vocabulary rather than by substring, and halts on any
label not in that vocabulary — so a silent misassignment of the kind in defect 3
cannot recur.

## 9. Tumour-vs-normal designation: three hazards for Amendment 4

Full audit in `output/atlas_tumour_designation_audit.csv`. All five atlases can
support the tumour-only restriction, but three carry traps that would corrupt the
epithelial fraction silently rather than loudly:

1. **GSE155698 contains 21 PBMC samples out of 41.** The series is 17
   `PDAC_TISSUE`, 3 `AdjNorm_TISSUE`, 17 `PDAC_PBMC`, 4 `Healthy_PBMC`. Blood has
   essentially no epithelium. If the tumour filter keys on "PDAC" rather than on
   tissue-vs-blood, the PBMC samples pass and drive the epithelial fraction toward
   zero — manufacturing the paper's own conclusion. Filter must keep
   `PDAC_TISSUE*` only.
2. **GSE183904 label strings are not clean.** One GSM reads
   `Peritonium tissue  (Tumor)` with a double space. Exact-match filtering
   mis-sorts it. Normalise whitespace before comparing. Amendment 4 excludes
   peritoneal and normal, leaving the 26 `Primary Gastric Tissue (Tumor)` samples.
3. **Peng is in GSA, not GEO or Zenodo.** The brief calls it "the Peng Zenodo
   release", but the paper's data-availability statement gives
   **GSA: CRA001160**, project **PRJCA001063** (Genome Sequence Archive). 24
   primary PDAC tumours and 11 control pancreases. If a Zenodo mirror is the
   intended source, its provenance and processing should be recorded, because it
   is not the primary deposit.

Also worth noting: GSE125449 has no tumour/normal field in its per-cell table at
all, but this is moot — all 19 samples are tumours (10 iCCA, 9 HCC), so the
restriction is a no-op there.

For GSE178341, `specimen_type` is T (129 GSMs) / N (52 GSMs) at GSM level. Whether
the published per-cell metatable carries that column per cell still needs
confirming when the file is opened.

## 8. `renv/library-local/` is tracked in git (23 MB of package binaries)

`.gitignore` excludes `renv/library/` but the local library I created to work
around the broken bioconda `TCGAbiolinksGUI.data` is at `renv/library-local/`,
which that pattern does not match. 37 files, 23 MB, entered in the starting
commit `361d287`. Also tracked from that commit: `MANIFEST.txt` (a GDC download
byproduct), `df.rds`, `results.rds`.

Not acted on — removing them means rewriting your commit. To fix:

```
cd ~/Documents/stat3-gi
printf 'renv/library-local/\nMANIFEST.txt\n*.rds\n!data/**/*.rds\n' >> .gitignore
git rm -r --cached renv/library-local MANIFEST.txt df.rds results.rds
git commit -m "Untrack installed R library and download byproducts"
```

The blobs stay in history unless you rewrite it, which is probably not worth it
at this size.

## 7. Malignancy is inferred by three different methods across the five atlases

Full audit in `output/atlas_malignant_annotation_audit.csv`. The compartment
estimand is the *malignant* epithelial fraction, so how each atlas decides which
epithelial cells are malignant is not a detail — it is the definition of the
numerator, and it is not consistent:

| Atlas | Labelled? | Method |
|---|---|---|
| GSE125449 (liver) | yes | inferCNV, per-cell threshold (score >80th pct AND corr >0.40) |
| GSE183904 (gastric) | yes | inferCNV (per-cell cutoff not stated in main text) |
| GSE178341 (CRC) | yes | **trained classifier** (>0.75 tumour / <0.25 normal); inferCNV only as cross-check |
| Peng (PDAC) | partial | inferCNV at cluster level, then carried by a **ductal type-2 proxy** |
| GSE155698 (PDAC) | **no** | none — no CNV analysis of any kind in the paper |

Two consequences worth deciding on before the sweep:

1. **GSE178341 is method-discordant.** Its primary malignancy call is not
   CNV-based, and the paper reports that ~11% of its likely-malignant cells show
   no substantial copy-number difference from normal (8% MMRp, 15% MMRd). Those
   cells would probably fail GSE125449's inferCNV cutoffs. A malignant-epithelial
   fraction computed in CRC is therefore not measuring quite the same quantity as
   one computed in liver.
2. **GSE155698 cannot serve this estimand as published.** Deriving malignancy
   labels ourselves would make PDAC the only cohort where the numerator is our
   own construction — a reviewer-visible asymmetry. The Peng release is the better
   PDAC option despite its cluster-level proxy.

Amendment 3's "at least two tissue-matched atlases" requirement interacts with
this: for PDAC there may not be two atlases with comparable malignancy calls.

## 6. Sci Data 2026 atlas has no malignant-cell annotation

Relevant to step 3 and detailed in the chat recommendation. Its level-1 labels
are lineage-level (Epithelial, Fibroblast, Endothelial, T/NK, B/Plasma, Myeloid,
Mast, Enteric Glial) with no malignant-vs-normal epithelial split and no
inferCNV/CopyKAT step in the methods. For a paper whose estimand is the
*malignant* epithelial fraction, that distinction has to come from somewhere.
