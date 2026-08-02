# Notes for review

Observations surfaced but NOT acted on, per the scope discipline rule.

## 16. DISCREPANCY IN REGISTERED TEXT — Amendment 13's stated deviation is 4x low

**The amendment's procedure and its registered tolerance both stand. One prose
figure inside it is contradicted by the full run, and I have not edited it.**

Amendment 13, assertion (iii), as registered:

> every recovered value lies within 0.1 of an integer BEFORE rounding.
> **Maximum observed deviation is 2.9e-3**, consistent with float32 storage.

The full Part A run, over all 139,415,620 non-zero entries:

```
(iii) max deviation from integer before rounding = 1.170e-02  (tol 0.1)
```

**1.17e-2 is roughly 4x the 2.9e-3 the amendment records.** The 2.9e-3 figure
appears to derive from a 400-cell sample taken while the recovery was being
designed, before it was run over the whole matrix.

**Nothing about the procedure changes.** The registered *tolerance* is 0.1, and
1.17e-2 sits an order of magnitude inside it — and far below the 0.5 at which
rounding could recover a wrong integer. Assertion (iii) passed on its own terms.
Recovery remains exact: assertion (ii) is 0.000e+00, i.e. recovered per-cell sums
equal the deposited `n_counts` bit for bit.

**Why this is recorded rather than fixed.** Amendment text is registered and goes
in verbatim; I do not edit it, including to correct an observation, because the
value of a preregistration is that its text is what it was when written. The
correction belongs here and in the paper's methods, not in the amendment.

**It already caused one concrete error.** I took 2.9e-3 from this amendment as
the basis for a tightened pre-rounding guard at `dev_iii < 1e-2` (audit finding
F03). That guard halted a run which had already passed Amendment 13's own
criterion — the observed 1.17e-2 exceeds 1e-2. The guard has been removed and
only the registered 0.1 gate remains. A figure in registered prose is not a
specification, and should not be promoted to a threshold; the tolerance the
amendment actually sets is the specification.

**If a stricter bound is wanted**, it should be registered as a new dated
amendment based on the full-corpus measurement (1.17e-2), not back-fitted to a
sample figure after the data have been seen.

## 15. REPORTING LIMITATION — GSE178341 has no endothelial compartment

**Not an amendment. The code matches the registration; this is a limitation in
what the per-compartment breakdown can be used for.**

Surfaced by the Implementation Auditor (F08). GSE178341's `clTopLevel` vocabulary
is, in full:

```
B | Epi | Mast | Myeloid | Plasma | Strom | TNKILC
```

There is **no endothelial label**. That atlas's endothelial cells sit inside
`Strom`, which the registered A.c map sends to `fibroblast_stromal`. So
GSE178341 reports `n_cells['endothelial'] == 0` while GSE125449 (`TEC`) and Peng
(`blood vessel endothelial cell`) both report a real endothelial compartment.

**Why no amendment is required.** Amendment 4 specifies six compartments mapped
from *each atlas's own level-1 labels*. GSE178341 has no level-1 endothelial
label, so mapping `Strom` → `fibroblast_stromal` is what the registration
prescribes, not a deviation from it. Adding an endothelial target would require
descending to `cl295v11SubFull`, which is a different annotation level than the
one Amendment 4 names.

**Why the estimand is unaffected.** f(π) depends only on the epithelial /
non-epithelial split. Moving cells between two *non-epithelial* compartments
changes neither the numerator nor the denominator — the identical algebra to
Amendment 11's proven-neutral celltype0/celltype1 refinement. `k` and every
variant are unchanged.

**What it does limit.** The per-compartment breakdown is **not comparable across
atlases for the `endothelial` and `fibroblast_stromal` rows**. Any statement of
the form "endothelial contribution is lower in colorectal" would be an artefact
of annotation granularity, not biology. Recorded in
`output/compartment_dominance_matrix.readme.txt` so the constraint travels with
the data rather than living only here.

## 14. Automated Reviewer findings, and a tooling gap in how they reach the assistant

**The tooling gap first, because it is the reportable part.** During this session
the Reviewer panel accumulated **15 checks / 10 unresolved findings**, but only
**2** were ever surfaced in the working conversation — both `fail`-level, both
fixed immediately. The other **10 are all `warn`-level and were never injected**;
the assistant only saw them after querying the findings surface directly at the
author's prompting. Background reviews do not interrupt on warns by design, so a
warn-level finding can sit unaddressed indefinitely while work continues on top
of it.

Consequence for this project: **assume nothing about review status without
querying it explicitly.** Two of the ten (§9/§10 below) concerned a numerical
claim underpinning Amendment 13 and would have gone into the record uncorrected.
Recorded here as a process finding: a review signal that does not reach the
worker is a silent failure mode, and this project's fail-closed discipline should
extend to the review channel itself.

### Reconciliation of the panel's counters

| Panel shows | Reconciles to |
|---|---|
| 15 checks | 15 rows total |
| 10 findings | 10 `warn` / `unaddressed` |
| 3 fixed | 5 `resolved` = 2 `fail` + 2 `pass` + 1 `warn`. The "3 fixed" counter appears to count the 2 fails plus 1 resolved warn; the 2 `pass` rows are verifications that never required a fix. |

### The 10 unresolved findings, with disposition

All are `warn` / low severity. None invalidates a Part A result; all concern
prose, citation or scope precision in saved artifacts.

| # | Claim | Disposition |
|---|---|---|
| 1 | REVIEW artifact said `keep` in `load_GSE125449` is "never used to check anything" | **Correct finding.** `keep` *is* used (`if (!any(keep)) halt(...)` and to filter `meta`). The real gap was narrower: no assertion on `sum(keep) == nrow(meta)`. Moot in effect — B2 was later disproved (zero barcode collisions), so no code changed. Artifact wording stands as-written history. |
| 2 | R8 said hardcoded `152` appears in "three vapply calls"; there were four | **Correct finding.** Off-by-one in the count only; the recommendation was right and **all four** were replaced with `length(PANEL_GENES)`. |
| 3 | Cited HANDOFF §11.6 for an *R*-kernel grant failure; §11.6 names only the python kernel | **Correct finding.** Sourcing overreach in chat prose. The remedy applied (use `bash`) is what §11.6 prescribes. No artifact affected. |
| 4 | FIXES.md said the band identity was "checked exhaustively over n = 1…400" | **Correct finding.** The *shipped* `verify_branch_bands()` loops 8 panel sizes, not 1..400, so the artifact conflates the shipped guard with an exploratory check. On the separate question of whether that sweep happened at all: it did, and this row originally asserted so without checking — a second reviewer finding (66688cca) flagged the unverified assertion, and the transcript was then searched. The sweep is `/tmp/bands.R`, `for (n in 1:400) { lo <- ceiling(0.05*n); hi <- ceiling(0.20*n); for (k in 0:n) ... }`, output `Q1 band-definition violations over n=1..400, all k: 0`. Claim true; the process that produced it was not sound, and the correction stands as the record of that. |
| 5 | "24 T + 11 N patients ✓" for Peng presented as freshly reproduced | **Correct finding.** In that audit those figures were quoted from `assert_n` source lines, not recomputed. They *were* independently confirmed later, in the Peng dry run and the full Part A run. |
| 6 | Commit `dceeb1c`'s message describes the channel file, which gitignore excluded from that commit | **Correct finding.** The file landed in `c58b889`. The message overstates that commit's contents. |
| 7 | Said "explicit exception rather than a bare `-f`", then used `git add -f` | **Correct finding.** Both were done: the `.gitignore` exception makes tracking durable, `-f` was belt-and-braces for the initial add. The prose implied an either/or. |
| 8 | R10 comment pointed at `output/partA_run.log` for the benchmark, which that file does not contain | **Correct finding. FIXED** — the comment now points at this file (§14) and the introducing commit. |
| 9 | NOTES said `sum(expm1(raw/X))` = "exactly 10000 for the cells checked" | **Correct finding. FIXED** — see §12. The exploratory gate was *absolute* 1e-3 (= relative 1e-7, below float32 precision), so 5/400 cells tripped it. Amendment 13 asserts *relative* 1e-6: 0/400 fail, full run max 3.42e-07. |
| 10 | Same "exactly 10000 per cell" phrasing, unhedged, repeated in a commit message | **Correct finding.** Same root cause as §9. The NOTES text is fixed; the commit message is immutable history and is corrected here instead. |

Findings 1–7 are recorded rather than edited: they describe chat prose or the
as-written state of superseded artifacts, and rewriting history to look cleaner
than it was would be worse than leaving an accurate record with its corrections
attached. Findings 8 and 9 named live text in files still in use, so both were
fixed.

## 13. CANDIDATE FINDING for the paper's discussion — deposit defects in three of five atlases

**Not an incidental problem with this study. A general, reportable result about
single-cell deposit quality.**

Three of the five atlases originally selected for Part A carried defects that
were invisible from their publications and only surfaced on direct inspection of
the deposits:

| Atlas | Publication implies | Deposit actually contains | Consequence |
|---|---|---|---|
| **GSE155698** (Steele, pancreatic) | cell-type annotation | count matrices only, no annotation | Removed, Amendment 9 |
| **GSE183904** (Kumar, gastric) | cell-type annotation | count matrices only, no annotation | Removed, Amendment 10 |
| **Peng** (Besca reprocessing) | `raw/X` raw counts | log1p-CP10K; no count layer anywhere | Recovered by inversion, Amendment 13 |

Two further deposit-level surprises in the two surviving GEO atlases, both
requiring amendments:

- **GSE125449**: the two deposited sets have different gene universes (20,124 vs
  19,572 rows) with no note in the publication (Amendment 12).
- **GSE178341**: 129 GEO tumour channels but 128 in the deposited matrices; one
  channel carries zero cells (recorded under Amendment 10).

The pattern is consistent: **the methods sections describe the analysis the
authors performed, not the artefact they deposited.** Amendment 4's claim that
lineage labels were "available and comparably defined in all five" atlases was
taken from the source publications and was wrong for two of five. Peng's case is
sharper still — the field name `raw/X` positively implies raw counts and does not
contain them.

For a paper whose thesis is that a widely used biomarker behaves differently once
its compartment composition is examined, this is a thematically aligned secondary
finding: **claims about public single-cell data are frequently not checkable from
the publication, and are sometimes contradicted by the deposit.** Worth a short
discussion paragraph with the table above. Every instance here is documented with
the amendment that resolved it, so the evidence is already assembled.

## 12. BLOCKER — Peng's `raw/X` is log1p-CP10K, not raw counts

**Part A halts at A.a/Peng. This is a data finding, not a code defect, and the
decision is not mine to make.**

The B4 guard added this session fired on the real file:

```
HALT [A.a/Peng]: matrix 'raw/X' is not integer-valued and is therefore NOT raw
counts (139415620 of 139415620 non-zero entries fractional; e.g. 0.6134,
1.2641, 1.2641)
```

**Every one** of the 139,415,620 non-zero entries is fractional. Diagnosis:

- `sum(expm1(raw/X))` = **10000 to within 3.4e-07 relative** across all 57,423
  cells → the layer is **log1p of CP10K-normalised** expression.

  Precision note, because an earlier revision of this section said "exactly
  10000" and that was too strong. The exploratory check used an *absolute* gate
  of 1e-3 on a sum of 1e4 — a relative tolerance of 1e-7, which is below
  float32 storage precision (~1.19e-07) — and 5 of 400 sampled cells tripped it,
  printing `FALSE`. That was a badly chosen gate, not a failed premise.
  Amendment 13's assertion (i) uses a *relative* 1e-6 tolerance, which 0 of the
  same 400 cells fail, and the full run passed it at a maximum of 3.42e-07 over
  every cell. The CP10K conclusion is unaffected.
- The `.h5ad` contains **no other count layer**: `X` is 2,033 × 57,423 (HVG
  subset, also normalised), there is no `/layers` group, and the only other
  sparse matrices are the neighbour graph (`obsp/connectivities`, `distances`).
- `obs/n_counts` and `obs/n_genes` are integer-valued and survive, so per-cell
  library size is known even though the count matrix is not deposited.

This directly contradicts the A.d premise recorded in analysis_plan.md v1.5 and
HANDOFF §8: pseudobulk requires summing **raw UMI counts**, and A.d explicitly
forbids any library-size correction before the `rowSums`. Normalising first
reweights every compartment to equal RNA content — the error the feasibility
addendum attributes to the HPA pilot. The Besca release cannot satisfy A.d as
written.

### Counts are *approximately* recoverable, but this is a reconstruction

`expm1(raw/X) * n_counts / 1e4` returns near-integers and the per-cell sum
reproduces `n_counts` exactly. Over 400 randomly sampled cells:

| Check | Result |
|---|---|
| Max deviation from integer | **2.9e-3** |
| Cells exceeding 1e-4 tolerance | **63 / 400** |
| Recovered per-cell sum vs `n_counts` | matches |
| Minimum recovered count | 1 |

The residual is consistent with float32 storage of the normalised values. So the
counts are recoverable to rounding, but the result is **inferred, not
deposited**, and rounding to integer is an analyst decision that changes the
values A.d sums.

### The options, none of which I should pick unilaterally

1. **Amend A.d to accept a reconstructed count matrix for Peng**, specifying
   `round(expm1(X) * n_counts / 1e4)`, disclosing that Peng's counts are
   reconstructed rather than deposited, and stating the direction of bias.
2. **Obtain genuine counts** from GSA CRA001160 (PRJCA001063) and re-derive the
   annotation, or find another Besca release carrying a count layer. Note
   HANDOFF §8 records that raw GSA was rejected *because it carries no
   annotation* — so this means re-annotating, which is a much larger change and
   a post-registration researcher degree of freedom.
3. **Drop Peng**, leaving two atlases. This breaks Amendment 5's "two of three"
   replication requirement and would collapse `k`, `k_all3` and `k_evalall`
   into near-identical quantities. Amendments 9 and 10 already removed two
   atlases; a third removal is a substantial change to the registered design.

**Recommendation: option 1**, because it is the only one that neither
re-annotates data nor removes an atlas, and the reconstruction error is far
below the granularity that could change a dominance call. But it is an amendment
to the estimand's input and must be recorded as such, with the 2.9e-3 figure and
the 63/400 disclosed.

**Note this guard is doing exactly its job.** Without B4, `read_h5ad_counts`
would have silently pseudobulked log1p-CP10K values for one of three atlases,
and no downstream assertion — monotonicity, dominance, the k ordering — could
have detected it.

## 11. BLOCKER — GSE125449 Set1 and Set2 have different gene universes

**Part A cannot run until this is decided. It is a prespecification gap, not a
bug, and the fix changes a reported quantity — so it is not mine to choose.**

Discovered 2026-08-01 while proving the R10 speedup on real data. `load_GSE125449`
halts at its own guard:

```
Set1 gene rows: 20124 | Set2 gene rows: 19572
identical(rownames(Set1), rownames(Set2)) -> FALSE
HALT [A.a/GSE125449]: Set1 and Set2 gene rows differ; cannot combine without a join
```

The guard is correct and fired exactly as intended — `cbind` on mismatched row
spaces would silently misalign genes against cells. Per the standing rule the run
was stopped rather than the assertion adjusted.

`analysis_plan.md` v1.5 specifies combining Set1 and Set2 but does not say how to
reconcile their gene universes. The three options are not equivalent:

| Option | Panel coverage | Note |
|---|---|---|
| Intersection (18,367 genes) | **143 / 152** panel, 146 / 155 reporting | Loses 6 panel genes measured in one set |
| Union with 0-fill | 149 / 152 panel, 152 / 155 reporting | **Violates A.d**: absent genes must be `NA`, never 0 |
| Union with NA-fill | 149 / 152 panel, 152 / 155 reporting | Per-gene denominators differ between sets; changes what f(π) means |

(Union counts are genes present in Set1 **or** Set2: 146 + 146 − 143 = 149 of the
152-gene panel, and 149 + 149 − 146 = 152 of the 155-gene reporting set. An
earlier revision of this table quoted 149/155 for both union rows, which mixed
the panel-union count with the reporting denominator.)

The 6 panel genes present in one set but lost to intersection: **CCL7, CRLF2,
CSF2, IL9R, ITGB3, LEP**. None is an origin-six gene, so SOCS3, MYC and IL6 are
unaffected either way.

Intersection is the most defensible — it is the only option that keeps a single
coherent denominator per cell and does not invent zeros — but it drops 6 of 152
panel genes from one of three atlases, which mechanically lowers each gene's
evaluability and can only move `k` downward. Given Amendments 9 and 10 already
pushed `k` down and both disclosed it as non-conservative, this is a third
same-direction effect and should be disclosed on the same footing.

**Recommendation: an amendment specifying intersection, recording the 6 lost
genes by name and the direction of bias.** Not implemented — awaiting a decision.

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
