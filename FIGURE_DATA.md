# Figure data provenance

Which committed file, and which columns, each figure draws from. Written so a
reader can check any plotted value against its source without reading the code.

**The rule the script enforces on itself:** for every figure, the frame handed
to ggplot is compared against a **fresh read of the source file performed inside
the assertion** — the caller cannot pass an in-memory object as the source.
`assert_plot()` matches rows on their keys, compares the numeric values at a
tolerance matched to each column's written precision, and **halts** on any
disagreement.

Before any figure is drawn, `verify_assert_can_fail()` corrupts one real value
and confirms the assertion rejects it, then confirms it passes the true frame.
An earlier version took the source as an argument and every call site passed the
object it had just subset, so the comparison was a value against itself and
could not fail — verified: a frame with a value corrupted to 999 passed. The
self-test exists so that defect cannot recur silently.

No figure recomputes a reported number. Two derivations are presentation
transforms, both checked: `1 - f` for the stacked bar, and `beta ± 1.96·se` for
the cohort intervals in Figure 2 — the latter asserted equal to the file's own
`HR_lo`/`HR_hi` on the log scale (max difference 5e-05).

**One declared exception.** Figure 1 needs the 10,000 *per-draw* pooled M1
values. Those exist only in `output/null_replicates.rds`;
`null_distributions.csv` holds their five-number summary, not the draws.
Section 0 of `11_figures.R` exports them to `figures/null_pooled_draws.csv` —
deliberately **not** into `output/`, which holds the audited artefact set — and
asserts the export reproduces the committed summary exactly — min, q25, median,
mean, q75, max, the observed percentile, the draw **count**, and a 101-point
percentile grid — before any figure reads it. The seven moments alone would not
constrain the shape of a histogram, which is the one property the figure
displays; the grid does. The export copies stored values; it does not recompute
them.

---

## Figure 1 — `fig1_null_distribution.png` (headline)

Pooled M1 log-HR of the real panel against 10,000 matched null signatures.

| Source | Columns used | Role |
|---|---|---|
| `figures/null_pooled_draws.csv` | `pooled_M1_logHR` (10,000 rows) | the histogram |
| `output/null_distributions.csv` | `observed_reported_08`, `observed_percentile`, `median`, `min`, `max` | the panel's line, its percentile, asserted |
| `output/null_pvalues.csv` | `p_two_sided_PRIMARY`, `N_requested` | the annotated p-value |
| `output/meta_analysis.csv` | `n_total` (row `M1`) | patient count in the caption |
| `output/survival_per_cohort.csv` | `events`, `meta_eligible` | event count in the caption |

The export also carries `pooled_M2_logHR` and `pooled_attenuation_total`; neither
is plotted here.

Rows filtered to `role == "primary_registered"` and
`statistic == "pooled_M1_logHR"` / `"p_crude_M1"` — the registered configuration
only. The size sensitivity and the exploratory post-hoc draw are **not** plotted.

Shaded band is the central 95% of the drawn null, computed from the plotted
values themselves; it is a description of what is drawn, not a reported number.

## Figure 2 — `fig2_forest_logHR.png`

Per-cohort score log-HR for M1 and M4, with the pooled estimate and GSE39582.

| Source | Columns used | Role |
|---|---|---|
| `output/survival_per_cohort.csv` | `cohort`, `model`, `beta`, `se`, `HR`, `HR_lo`, `HR_hi`, `fitted`, `meta_eligible` | seven discovery rows for M1, six for M4 — CHOL is `fitted = FALSE` for M2–M4 |
| `output/meta_analysis.csv` | `analysis`, `est`, `ci_lo`, `ci_hi` | the pooled estimate (rows `M1`, `M4`) |
| `output/validation_gse39582_models.csv` | `cohort`, `model`, `beta`, `se`, `fitted` | the external cohort |

Intervals for cohort rows are `beta ± 1.96·se`, computed from the file's own
`beta` and `se`. The pooled interval is the file's Hartung-Knapp `ci_lo`/`ci_hi`,
**not** recomputed. `HR`/`HR_lo`/`HR_hi` are not plotted; they are the cross-check
— `exp()` of every drawn bound is asserted equal to the committed HR bound, so a
divergence between the drawn interval and the file's would halt.

**CHOL** carries `meta_eligible = FALSE`; it is drawn with an open marker,
labelled "descriptive, not pooled", and never enters the pooled row.
**GSE39582** sits in its own facet band labelled "External" and is explicitly not
pooled with discovery (Amendment 16).

## Figure 3 — `fig3_compartment_attribution.png`

Epithelial share for the origin six across the three atlases.

| Source | Columns used | Role |
|---|---|---|
| `output/origin_six_compartment.csv` | `gene`, `atlas`, `f_at_0.30`, `f_at_0.50`, `f_at_0.70`, `qualifying`, `dom_rate` | the stacked bars |

Each bar stacks `f_at_<pi>` (epithelial) against `1 - f_at_<pi>` (everything
else), at all three points of the registered 30–70% purity band rather than one
point inside it. The dashed line is the 50% dominance threshold.

`qualifying` drives the axis marking: BCL2, MMP9 and HGF carry ` *`, because they
failed criterion B under Amendment 2 and enter no k variant.

`dominant` drives a diamond above the atlas-gene cells the compartment analysis
calls epithelial-dominant, and the count of diamonds is asserted equal to the
file's `dominant` column. This is load-bearing: nine of the 54 bars clear the 50%
line at some grid point — MYC in all three atlases and BCL2 in all three at
π = 0.70 — while only MYC × GSE178341 is dominant across the whole band. Without
the diamond, the line and the dominance claim would contradict each other.

## Figure 4 — `fig4_attenuation_forest.png`

attenuation_total per cohort with the pooled estimate.

| Source | Columns used | Role |
|---|---|---|
| `output/attenuation_per_cohort.csv` | `cohort`, `attenuation_total`, `att_total_lo`, `att_total_hi`, `att_total_se` | the six cohort rows |
| `output/meta_analysis.csv` | `est`, `ci_lo`, `ci_hi`, `tau2`, `I2` (row `attenuation_total`) | the pooled row and the caption |

Cohort intervals are the file's own **percentile** intervals from the 2,000
paired bootstrap resamples — not recomputed, and not the Wald form. CHOL is
absent from this file by construction (M2–M4 are not fitted there under
Amendment 8) and so does not appear; the caption says why.

## Figure 5 — `fig5_score_by_cms.png` (the deferred B.o.1 item)

Score distribution by consensus molecular subtype, COAD and READ.

| Source | Columns used | Role |
|---|---|---|
| `output/cms_calls.csv` | `cohort`, `barcode`, `cms` | subtype per patient |
| `output/scores_per_patient.csv` | `barcode`, `score` | the score |
| `output/cms_distribution.csv` | `cohort`, `cms`, `n` | the per-group n, asserted |

Joined on `barcode`; the script halts unless the join is 1:1 across all 620 CMS
calls. The per-group counts drawn are asserted equal to `cms_distribution.csv`,
so the violins cannot be drawn from a different set of patients than the
committed distribution table. `unclassified` is retained as a level (B.o).

The χ² and R² figures in the caption are **read at render time** from
`output/cms_tertile_crosstab.csv` (`chisq`, `df`, `p_asymptotic`) and
`output/cms4_association.csv` (`r2_score_on_cms`) — not typed as literals, so
they cannot go stale if the analysis is re-run. An earlier draft named
`cms_adjusted_models.csv` for the R², which has no such column.

---

## Not plotted, deliberately

- **FU-iCCA, GSE62254, ICGC** — unresolved access blockers (`DATA_NEEDED.md`).
  Amendment 16's PRIMARY validation analysis (STAT3 pY705 concordance) has not
  been run, so no figure speaks to whether the score tracks STAT3
  phosphorylation.
- **The size sensitivity and exploratory post-hoc null configurations** — the
  registered primary is the one figure 1 shows. Both others are in
  `null_pvalues.csv` and `null_distributions.csv` under their `role`.
- **`cit.molecularsubtype`** — an additional, explicitly non-registered
  orthogonality check; kept out of the CMS figure so the two cannot be read as
  one analysis.

---

# Refined set — `figures/refined/`, produced by `11b_figures_refined.R`

`11_figures.R` and `figures/` are **untouched**; the refined set is additive. No
reported number differs between the two sets — the refinement is presentation
only, and the same `assert_plot()` contract binds every figure. **29 assertions
pass per run**, each comparing the ggplot frame against a *fresh read* of the
named file.

## The assertion, and why it can fail

`assert_plot(fig, plotted, file, keys, cols)` re-reads `file` **inside itself**.
The caller supplies only the plotted frame — it cannot pass an in-memory object
as the "source", which is what made the first version of this guard in
`11_figures.R` tautological. `verify_assert_can_fail()` runs before any figure is
drawn: it corrupts one real value, requires rejection, then requires the true
frame to pass. If either half misbehaves the script halts having drawn nothing.

Keys must be **unique in the source**. This caught a real error while the refined
set was being built: the chi-square statistic in `cms_tertile_crosstab.csv` is
repeated on all 12 tertile × CMS rows of each cohort, so `cohort` alone is not a
key. The assertion refused it; the figure now asserts on the full
`cohort × tertile × cms` key and confirms the statistic is constant within cohort
before collapsing.

## Per figure

| Figure | Source file(s) | Columns asserted |
|---|---|---|
| `figure_main.png` A | `null_distributions.csv`, `null_pvalues.csv`, `meta_analysis.csv`, `survival_per_cohort.csv` | `observed_reported_08`, `observed_percentile`, `median`, `min`, `max`; `events`, `n` per cohort |
| `figure_main.png` B | `survival_per_cohort.csv`, `meta_analysis.csv`, `validation_gse39582_models.csv` | `beta`, `se`; `est`, `ci_lo`, `ci_hi` |
| `figure_main.png` C | `fuicca_correlations.csv`, `fuicca_per_patient.csv` | `n`, `r`, `ci_lo`, `ci_hi` |
| `figure_crc_pool.png` | `exploratory_colorectal_pooled.csv`, `exploratory_colorectal_inputs.csv` | `est`, `HR`, `HR_lo_wald`, `HR_hi_wald`, `I2`, `tau2`; `beta`, `se`, `n`, `events` |
| `figure_rppa_concordance.png` | `exploratory_rppa_control_correlations.csv`, `exploratory_rppa_control_pooled.csv` | `n`, `r`, `ci_lo`, `ci_hi`; `r_pooled`, `ci_lo_wald`, `ci_hi_wald`, `I2` |
| `fig1_null_distribution.png` | as `figure_main` A | as above |
| `fig2_forest_logHR.png` | as `figure_main` B, both models | as above |
| `fig3_compartment_bars.png` | `origin_six_compartment.csv` | `f_at_0.30`, `f_at_0.50`, `f_at_0.70` |
| `fig3_compartment_heatmap.png` | `origin_six_compartment.csv` | same — both variants share one data-prep function |
| `fig4_attenuation_forest.png` | `attenuation_per_cohort.csv`, `meta_analysis.csv` | `attenuation_total`, `att_total_lo`, `att_total_hi`; `est`, `ci_lo`, `ci_hi`, `tau2`, `I2` |
| `fig5_score_by_cms.png` | `cms_calls.csv`, `scores_per_patient.csv`, `cms_distribution.csv`, `cms_tertile_crosstab.csv`, `cms4_association.csv` | group counts vs `n`; `chisq`, `df`, `p_asymptotic`, `n`; `r2_score_on_cms` |
| `fig6_fuicca_y705_scatter.png` | as `figure_main` C | as above |
| `figS_forest_M4.png` | as `fig2`, M4 only | as above |

## Numbers in titles, subtitles and captions are asserted too

A number interpolated into a caption is still a reported number. Everything
`sprintf`'d into text is read from a committed file and covered by an
`assert_plot` call on the same frame — including the CMS chi-square, degrees of
freedom, p-values and R², the colorectal pooled HR and its interval, the RPPA
pooled correlations, and the FU-iCCA correlation and n.

**One number had no committed column**: the event total in panel A's subtitle
(595). `meta_analysis.csv` carries `n_total` but no event total, so the
per-cohort `events` are asserted individually and the sum is taken from those,
with `sum(n)` cross-checked against the committed `n_total`.

## The one declared exception, unchanged from the original set

The 10,000 per-draw pooled M1 values exist only in `null_replicates.rds`;
`null_distributions.csv` holds their summary. `export_null_draws()` writes them to
`figures/refined/null_pooled_draws.csv` and asserts the export reproduces the
committed `min`, `q25`, `median`, `mean`, `q75`, `max`, `n_null` **and** the
observed percentile, and that the percentile grid is monotone. The histogram then
reads the CSV. This is the only figure input not itself a committed output.

## Two run-order discrepancies, resolved in favour of the files

- The run order says the null is "centred at 0.118". `null_distributions.csv`
  gives median **0.117486** and mean 0.117481, both rounding to **0.117**. Panel A
  reads the file and prints 0.117.
- The run order names `exploratory_crc_pooled.csv`; the committed file is
  `exploratory_colorectal_pooled.csv`. Its HR 0.9766 (0.8756–1.0892) matches the
  quoted "0.977 (0.876–1.089)" exactly — that is the **Wald** interval, which the
  figure quotes as the conservative one.

## Claim-titles: two were overstated and were corrected

Every title must be true of every plotted element. Testing them against the data
found two that were not:

- The RPPA figure said the phosphosite **"disagrees with"** total STAT3. The
  pooled r is 0.139 with a Wald interval of **−0.064 to 0.331**, which spans
  zero. An interval spanning zero cannot establish disagreement — it fails to
  establish agreement. Now: *"shows no established agreement with"*.
- The colorectal figure said pooling **"excludes"** a per-SD hazard increase above
  9% while its own caption said a bound is "not a proof of absence" — the figure
  contradicted itself. Now: *"the data are compatible with at most a 9% per-SD
  hazard increase"*.

The compartment title ("MYC is the only gene epithelial-dominant across the whole
30–70% band") is true of the **dominance calls**, which is what the mark shows:
`origin_six_compartment.csv` has `dominant = TRUE` for MYC in GSE178341 and
nothing else. Nine of the 54 fractions clear 50% at *some* grid point — MYC in all
three atlases, BCL2 in all three at π = 0.70 — and the bars show that, which is
why both variants mark the file's own call rather than relying on the 50% line,
and why the caption states the whole-band rule explicitly.

## Which figure-3 variant reads better

**The heatmap.** With 54 values it prints every number in-cell, so nothing has to
be estimated against an axis; the colour is centred on the 50% dominance
threshold, so "above or below the rule" is visible at a glance; and the dominance
call is an outlined cell rather than a diamond, which occludes nothing. The
stacked-bar version spends most of its ink on the "other compartments" remainder,
which is 1 − f and carries no independent information, and its dominance diamond
has to float above the bar to avoid overprinting. Both are produced; the heatmap
is the one I would ship.
