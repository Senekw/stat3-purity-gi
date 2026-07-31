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
