# REFERENCE VERIFICATION — all 25 entries checked against real records

**Report only. The manuscript was not edited.**

Method: each entry was resolved to a **DOI record via the Crossref REST API** and
**cross-checked against PubMed** (E-utilities). For every entry I compared the
resolved record's **title, first author, author count, journal, year, volume and
page range** against what the manuscript states. A DOI alone is not verification —
I recall DOIs imperfectly — so **every resolved record was title-matched against
the manuscript's own title text**, and an entry counts as verified only when the
titles agree. Token-overlap scores are given where they are below 1.00.

**Note on the napabucasin precedent.** The instruction refers to verifying "the
napabucasin paper's 12 identifiers" as I did before. **I have no record of that
work** — it is not in this session, not in my archived transcript, and not in any
committed file or artifact I can reach. I did not assume its results. Entries 7–11
were verified fresh here, on the same basis as the rest.

## Verdict

| | |
|---|---|
| Entries resolved to a real record and title-matched | **24 of 25** |
| Entry 25 | correctly *not* resolvable — under review, and **already formatted correctly** |
| Entries with a **wrong** stated value | **0** |
| Entries with a **year ambiguity** | **1** (R20) |
| Entries whose **title is rewritten, not merely shortened** | **1** (R9) |
| Entries with a **spelling change to a published title** | **1** (R14) |
| Entries **missing volume and pages** | **12 of 25** |
| Entries **missing a DOI** | **24 of 25** |

**No entry cites a paper that does not exist, and no entry attributes a paper to the
wrong first author.** The reconstruction is substantively sound; what it lacks is
completeness of fields, and it has three specific defects listed below.

---

## The three defects to fix

### R9 — the title is a rewrite, not the published title

| | |
|---|---|
| Manuscript | *"Paclitaxel with or without napabucasin in pretreated advanced gastric or gastro-oesophageal junction adenocarcinoma (BRIGHTER)"* |
| **Published** | *"Randomized, Double-Blind, Placebo-Controlled Phase III Study of Paclitaxel ± Napabucasin in Pretreated Advanced Gastric or Gastroesophageal Junction Adenocarcinoma"* |

The manuscript renders `±` as "with or without" and moves the design descriptors.
**"without" appears nowhere in the published title.** Everything else about the
entry is right — Clin Cancer Res 2022;28:3686-94 matches the record exactly
(`10.1158/1078-0432.CCR-21-4021`, PMID 35833783). This is the one entry where the
title itself must be replaced rather than extended.

### R14 — British spelling applied to a US journal's title

Manuscript: *"**Tumour** cell biodiversity drives…"* — published as *"**Tumor** Cell
Biodiversity Drives Microenvironmental Reprogramming in Liver Cancer"* (Cancer Cell,
US spelling). Also **missing volume and pages: 2019;36:418-430.e6**
(`10.1016/j.ccell.2019.08.007`).

Note this is a real distinction, not pedantry about house style: **R5 and R18
legitimately use "tumour"** because Nature Communications published them that way.
The spelling should follow each published title, not be normalised across the list.

### R20 — year is ambiguous, and both answers are defensible

TRRUST v2, Nucleic Acids Res 46:D380-D386. **PubMed says 2018** (cover date
2018 Jan 4); **Crossref says 2017** (online publication). This is the usual NAR
database-issue split. The manuscript says 2018, which matches PubMed and the
printed issue — **I would leave it at 2018** and note that a reference checker
comparing against Crossref may query it.

---

## Full verification table

`ov` = title token overlap with the published record (blank = 1.00, exact).
**Bold** marks a field the manuscript does not currently state.

| # | First author (n) | Published title | Journal | Year | Vol | Pages | DOI | PMID | ov |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Venet (3) | Most Random Gene Expression Signatures Are Significantly Associated with Bre | PLoS Comput Biol | 2011 | **7** | **e1002240** | `**10.1371/journal.pcbi.1002240**` | 22028643 |  |
| 2 | Isella (18) | Stromal contribution to the colorectal cancer transcriptome | Nat Genet | 2015 | **47** | **312-319** | `**10.1038/ng.3224**` | — |  |
| 3 | Calon (19) | Stromal gene expression defines poor-prognosis subtypes in colorectal cancer | Nat Genet | 2015 | **47** | **320-329** | `**10.1038/ng.3225**` | 25706628 |  |
| 4 | Guinney (40) | The consensus molecular subtypes of colorectal cancer | Nat Med | 2015 | **21** | **1350-1356** | `**10.1038/nm.3967**` | — |  |
| 5 | Aran (3) | Systematic pan-cancer analysis of tumour purity | Nat Commun | 2015 | 6 | — | `**10.1038/ncomms9971**` | 26848121 |  |
| 6 | Liu (740) | An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality S | Cell | 2018 | **173** | **400-416.e11** | `**10.1016/j.cell.2018.02.052**` | — |  |
| 7 | Jonker (23) | Napabucasin versus placebo in refractory advanced colorectal cancer: a rando | The Lancet Gastroenter | 2018 | 3 | 263-270 | `**10.1016/s2468-1253(18)30009-8**` | 29397354 |  |
| 8 | Shah (20) | Napabucasin Plus FOLFIRI in Patients With Previously Treated Metastatic Colo | Clinical Colorectal Ca | 2023 | 22 | 100-110 | `10.1016/j.clcc.2022.11.002` | 36503738 |  |
| 9 | Shah (22) | Randomized, Double-Blind, Placebo-Controlled Phase III Study of Paclitaxel ± | Clinical Cancer Resear | 2022 | 28 | 3686-3694 | `**10.1158/1078-0432.ccr-21-4021**` | — | 0.83 |
| 10 | Bekaii-Saab (19) | Napabucasin plus nab-paclitaxel with gemcitabine versus nab-paclitaxel with  | eClinicalMedicine | 2023 | 58 | 101897 | `**10.1016/j.eclinm.2023.101897**` | 36969338 |  |
| 11 | Froeling (16) | Bioactivation of Napabucasin Triggers Reactive Oxygen Species–Mediated Cance | Clinical Cancer Resear | 2019 | 25 | 7162-7174 | `**10.1158/1078-0432.ccr-19-0302**` | 31527169 |  |
| 12 | Johnson (3) | Targeting the IL-6/JAK/STAT3 signalling axis in cancer | Nat Rev Clin Oncol | 2018 | 15 | 234-248 | `**10.1038/nrclinonc.2018.8**` | — |  |
| 13 | Dong (32) | Proteogenomic characterization identifies clinically relevant subgroups of i | Cancer Cell | 2022 | 40 | 70-87.e15 | `**10.1016/j.ccell.2021.12.006**` | 34971568 |  |
| 14 | Ma (16) | Tumor Cell Biodiversity Drives Microenvironmental Reprogramming in Liver Can | Cancer Cell | 2019 | **36** | **418-430.e6** | `**10.1016/j.ccell.2019.08.007**` | 31588021 | 0.89 |
| 15 | Pelka (72) | Spatially organized multicellular immune hubs in human colorectal cancer | Cell | 2021 | 184 | 4734-4752.e20 | `**10.1016/j.cell.2021.08.003**` | — |  |
| 16 | Peng (23) | Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant pro | Cell Res | 2019 | 29 | 725-738 | `**10.1038/s41422-019-0195-y**` | 31409908 |  |
| 17 | Marisa (24) | Gene Expression Classification of Colon Cancer into Molecular Subtypes: Char | PLoS Med | 2013 | **10** | **e1001453** | `**10.1371/journal.pmed.1001453**` | — |  |
| 18 | Yoshihara (15) | Inferring tumour purity and stromal and immune cell admixture from expressio | Nat Commun | 2013 | 4 | — | `**10.1038/ncomms3612**` | — |  |
| 19 | Keenan (10) | ChEA3: transcription factor enrichment analysis by orthogonal omics integrat | Nucleic Acids Research | 2019 | **47** | **W212-W224** | `**10.1093/nar/gkz446**` | 31114921 |  |
| 20 | Han (20) | TRRUST v2: an expanded reference database of human and mouse transcriptional | Nucleic Acids Research | 2017 | **46** | **D380-D386** | `**10.1093/nar/gkx1013**` | 29087512 |  |
| 21 | Liberzon (6) | The Molecular Signatures Database Hallmark Gene Set Collection | Cell Systems | 2015 | **1** | **417-425** | `**10.1016/j.cels.2015.12.004**` | 26771021 |  |
| 22 | Viechtbauer (1) | Conducting Meta-Analyses in R with the metafor Package | J. Stat. Soft. | 2010 | 36 | — | `**10.18637/jss.v036.i03**` | — |  |
| 23 | Eide (4) | CMScaller: an R package for consensus molecular subtyping of colorectal canc | Sci Rep | 2017 | **7** | **—** | `**10.1038/s41598-017-16747-x**` | — |  |
| 24 | Hartung (2) | On tests of the overall treatment effect in meta‐analysis with normally dist | Statistics in Medicine | 2001 | 20 | 1771-1782 | `**10.1002/sim.791**` | 11406840 |  |
| 25 | Lee (1) | Napabucasin across four registrational Phase 3 trials in advanced gastrointestinal cancers | Clin Colorectal Cancer | — | — | — | — | — | n/a |

### Entry 25 — already correct, no change needed

The manuscript reads:

> Lee SGP. Napabucasin across four registrational Phase 3 trials in advanced
> gastrointestinal cancers. Clin Colorectal Cancer [under review].

This is exactly the required form: **journal named, no volume, no pages, no DOI, and
the status marked**. I confirmed it is **not** in PubMed or Crossref, which is what
"under review" should look like — an entry that resolved would mean it had actually
been published. Nothing to change. Two notes: if the journal accepts it before
submission, "[under review]" becomes "in press" and a DOI may exist; and this is a
**self-citation of unpublished work**, which some journals require be declared in
the cover letter.

---

## Missing-field inventory

**12 entries state no volume or pages** — 1, 2, 3, 4, 6, 14, 17, 19, 20, 21, 23,
and 25 (25 correctly so). The verified values are in the table above in **bold**.

**24 of 25 state no DOI.** Only R8 carries one. Every DOI is in the table above; all
24 were confirmed to resolve.

**Four titles are shortened by dropping a subtitle** — R7, R8, R10, R17. This is
acceptable in many house styles (e.g. R7 omits ": a randomised phase 3 trial"), so I
have not called it an error, but if the target journal requires full titles, all four
need extending. **R8 and R10 lose the trial names** CanStem303C and the
adaptive-superiority design descriptor from the title proper, though R8 keeps
"CanStem303C" in the manuscript's parenthetical.

**One record-side artifact, not a manuscript error:** Crossref returns R22's title
with italics markup collapsed (`inRwith themetaforPackage`). The manuscript's
rendering is correct; I have shown the corrected form in the table.

---
---

# PART 2 — In-figure claim-titles and journal submission

**Report only. Nothing was changed.** `11b_figures_refined.R` is untouched and the
11 figures in `figures/refined/` are as committed at `4187f8f`.

## The answer in one line

**Yes, the titles would need to come out, and it is a contained change — but the
`assert_plot` machinery needs no modification at all**, because it never touches
title text. What *does* need care is the two-line layout geometry and the material
that currently lives only in a title.

## What is actually there

**Ten claim-titles**, across nine figure functions plus the composite:

| Line | Figure | Title |
|---|---|---|
| 196 | `panel_null` | "The panel is indistinguishable from a matched random signature" |
| 270 | `panel_forest` | (conditional on `models == "M4"`) |
| 310 | `panel_fuicca` | "The score tracks measured STAT3 pY705 in FU-iCCA" |
| 364 | `fig_crc` | "Pooled across three colorectal cohorts, the data are compatible…" |
| 442 | `fig_rppa` | "The RPPA phosphosite shows no established agreement with total STAT3…" |
| 484 → 509, 540 | `fig3_bars`, `fig3_heat` | `TITLE3`, **shared by both variants** |
| 586 | `fig4_atten` | "Adjustment for purity and stroma is compatible with no reduction…" |
| 634 | `fig5_cms` | "The score is highest in CMS4, the mesenchymal subtype" |
| 679 | `figure_main` | "A preregistered STAT3 score tracks its target yet performs near the median…" |

Counts: **18 `title=`, 9 `subtitle=`, 9 `caption=`, 3 `draw_label`**.

## `assert_plot` needs no modification — verified, not assumed

Its signature is `assert_plot(fig, plotted, file, keys, cols, tol)`. I checked all
20 call sites: **not one passes a title, subtitle or caption string.** The
assertions compare *plotted data frames* against a fresh read of the source CSV.
Removing a title changes no argument to any of them.

Per-function counts confirm the separation — the assertions live in the same
functions as the titles but are independent of them:

| Function | assert_plot calls | Sets a title |
|---|---|---|
| `panel_null` | 2 | yes |
| `panel_forest` | 3 | yes |
| `panel_fuicca` | 1 | yes |
| `fig_crc` | 2 | yes |
| `fig_rppa` | 2 | yes |
| `fig4_atten` | 2 | yes |
| `fig5_cms` | 2 | yes |
| `compartment_data` | 1 | no |
| `verify_assert_can_fail` | 4 | no |

**All 29 `ok` lines (26 assertions + 3 other checks) would still pass unchanged.**

## What the change would actually involve

**1. The mechanical part is small.** Nine `labs(title=…)` removals plus one
`draw_label` in the composite. `TITLE3` is a single constant shared by both fig-3
variants, so that is one edit covering two figures. The cleanest form is a
parameter — the file already has the precedent in `panel_forest(compact = FALSE)`
and `standalone = FALSE` — so `titles = FALSE` could suppress them without
deleting the text, keeping the claim-titles available for talks and preserving them
in version control.

**2. The layout is coupled to title height, and this is the real work.** Four
titles are two-line `paste0(...\n...)` strings, sized to 180 mm at the current font.
`figure_main` composes with **`rel_heights = c(0.13, 0.45, 0.42)`** where the
`0.13` is the title strip — removing the title without re-tuning that vector leaves
a band of white space above panel A. Every standalone figure would need its
`save_fig(..., h)` height revisited; these were tuned by eye against the rendered
raster, so they need re-checking the same way, not arithmetic.

**3. Some titles carry content that exists nowhere else — this is the risk.** I
checked title/subtitle/caption number overlap per figure. Mostly the titles are
*qualitative* claims while the numbers live in the subtitle and caption, so deleting
a title loses no value. But two cases need the text moved rather than dropped:

- **`fig3`/`TITLE3`** is the only place stating the **whole-band rule** — that MYC is
  dominant *across the entire 30–70% band*, and only in colorectal tissue. Delete it
  and the 50% reference line becomes readable as a simple threshold, which is exactly
  the misreading the line was annotated to prevent (nine of 54 f values clear 0.50
  somewhere; only MYC/GSE178341 clears it across the band). **This sentence must move
  into the formal legend, not be deleted.**
- **`fig4`** and **`fig_crc`** titles carry the compatibility framing ("compatible
  with no reduction", "compatible with at most a 9% increase") that the corrected
  wording was specifically chosen for. Their captions already state the interval, so
  the information survives — but the *framing* is what stops a reader treating a
  point estimate as an established absence, and it should be reproduced in the legend.

**4. Do not simply delete `figure_main`'s title.** It is the paper's headline claim
and the only place the two halves of the result are stated together ("tracks its
target yet performs near the median"). In a journal that strips in-figure titles it
belongs in the legend's first sentence.

## Recommendation

Wait until the target journal is known. Requirements differ — some strip in-figure
titles, others accept a short bold lead-in — and the layout re-tuning is wasted if
done twice. When it is known: add a `titles = FALSE` switch rather than deleting
text, migrate the `TITLE3` whole-band sentence and the two compatibility framings
into the legends first, then re-tune `rel_heights` and the `save_fig` heights and
re-render, checking each raster by eye. **The assertions come along unchanged**, and
`FIGURE_DATA.md` would need its claim-title table updated to say the titles now live
in the legends.
