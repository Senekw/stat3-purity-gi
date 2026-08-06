# FINDINGS_REVIEW — all unresolved reviewer findings, verified against committed files

`host.findings()` returns **5 unresolved**: 3 warn/medium, 2 warn/low. **No
fail-level finding is open.** Each was verified against the committed files or by
execution, not conceded or disputed from memory.

## Bottom line first

> **⚠ One finding touches a number that appears in `PROJECT_SUMMARY.md` and
> `NOTES_FOR_REVIEW.md` — F2. It is CORRECT that the number was never computed at
> the time it was asserted. I have now computed it, and the assertion is TRUE:
> k = 46/43/43 over the 152/143/140 gene lists → 30.3%, 30.1%, 30.7%, all above
> Amendment 3's 20% floor. The branch is BUILT under every denominator. The claim
> stands; what was missing was the check, which is now on record below.**

**No finding changes any reported number, figure value, or the abstract.** Four of
the five are about *how a claim was evidenced or counted*, not about whether it is
right.

| # | ID | Sev | Status | Touches a reported number? |
|---|---|---|---|---|
| F1 | `a9460823` | warn/medium | **CORRECT BUT NOT ACTED ON** | no |
| F2 | `bdb89f85` | warn/medium | **AFFECTS A REPORTED NUMBER** — now verified true | **yes** — `PROJECT_SUMMARY.md`, `NOTES_FOR_REVIEW.md` |
| F3 | `876d6896` | warn/low | **CORRECT — now fixed** | no |
| F4 | `2c22af6f` | warn/low | **INCORRECT** — refuted by execution | no |
| F5 | `d93d3361` | warn/medium | **CORRECT BUT NOT ACTED ON** | no |

---

## F1 — `a9460823` · warn/medium · CORRECT BUT NOT ACTED ON

**Claims.** The PAR-1 audit finding's biological claim — that PAR1 genes are
male-biased by roughly 40%, and IL9R specifically shows strong male bias — rests on
web-search results whose text is redacted from the persisted transcript. Only
titles and URLs survive, so the specific figures are not verifiable from the
record. Called load-bearing because PAR-1 reversed the script's stated rationale
for the narrow reading.

**Correct?** **Yes, on the point it makes.** Verified in the committed files:

- `panel_definition.md:603` states the *mechanism* — *"PAR1 genes have a functional Y homolog and escape X-inactivation, so both sexes…"* — and `NOTES_FOR_REVIEW.md:1305-1309` records the correction to my original rationale. **The mechanism is committed and is sound.**
- The **quantitative** claim (~40% male bias) appears **nowhere** in any committed file. The only "40%" in the repository is the unrelated replication-bar figure (`analysis_plan.md:500`, `panel_definition.md:416`, two-of-five = 40%).

So the argument that survives in the repository is mechanistic, not quantitative —
which is the weaker and more defensible form. The redacted number was never
promoted into a committed claim.

**Why it is not acted on.** The finding asks for provenance on a figure that does
not appear in any deliverable. Nothing to correct; recorded here so that if anyone
later wants to *add* the quantitative claim to the manuscript, they know it needs a
citable source rather than a redacted search result.

**Does it touch a reported number?** No. But note the *downstream* connection: the
PAR-1 reasoning informs the broad-vs-narrow rule-3.3 reading, and that choice
selects `primary_140` over `sensitivity_143` as the scoring list
(`07_score.R:56-57`). The **choice** is consequential; the **~40% figure** is not
what justifies it in the committed record.

---

## F2 — `bdb89f85` · warn/medium · ⚠ AFFECTS A REPORTED NUMBER — verified true

**Claims.** The summary prose asserted *"43/140 = 30.7% — BUILT under every
denominator in play"* as established fact, while my own PAR-2 finding said k under
the broad reading *"is never computed… unverified"*. An internal contradiction:
presenting as verified exactly what a submitted finding flagged as unchecked.

**Correct?** **Yes — the contradiction was real.** Verified:

- `output/k_estimates.csv` commits **k = 46** over the **152-gene inferential set only**. It contains no k for the 140- or 143-gene lists.
- `04_lock_gene_list.R` asserts list *membership* (`stopifnot(length(final_narrow)==143L, all(final %in% final_narrow), …)` at line 270) but **never** asserts k under either reading.
- The claim is stated in `NOTES_FOR_REVIEW.md:1331-1333` — *"k is 43 under both readings, and the Amendment 3 branch is BUILT under 152, 143 and 140 alike (30.3%, 30.1%, 30.7% of panel, all above the 20% floor)"* — and `PROJECT_SUMMARY.md:121` states the 152-denominator form, *"k = 46 of 152 is 30.3%"*.

**So I computed it, from committed files only** — `output/compartment_dominance_matrix.csv`
(465 rows = 155 genes × 3 atlases) intersected with the two locked gene lists,
counting a gene as dominant where it is dominant in ≥ 2 atlases and treating `NA`
as not-evaluable rather than as `FALSE`:

| Denominator | k | % of panel | Above Amendment 3's 20% floor? |
|---|---|---|---|
| 152 (inferential set) | **46** | **30.3%** | yes |
| 143 (sensitivity list) | **43** | **30.1%** | yes |
| 140 (primary list) | **43** | **30.7%** | yes |

The k = 46 over 152 **reproduces `k_estimates.csv` exactly**, which validates the
derivation. **The asserted numbers are all correct and the branch is BUILT under
every denominator.** The finding was right that the check was missing; the claim it
questioned turns out to be true.

**Status.** The number is confirmed, not changed. The gap was evidential, and this
document closes it.

---

## F3 — `876d6896` · warn/low · CORRECT — now fixed

**Claims.** A COHORT-ID finding attributed *"seven named cohorts (COAD, READ, STAD,
ESCA, LIHC, CHOL, PAAD per Amendments 8 and 10)"* to the prespecification, but
`panel_definition.md` names only six and **contains no occurrence of PAAD at all**.
The finding's substantive point — that cohort identity is checked by count, not by
name — is otherwise sound.

**Correct?** **Yes, exactly right.** Verified directly:

- `grep -c "PAAD" panel_definition.md` → **0**.
- The TCGA codes present in `panel_definition.md` are exactly **CHOL, COAD, ESCA, LIHC, READ, STAD** — six.
- PAAD *is* prespecified, but in `analysis_plan.md` (line 767 endpoint table, lines 798 and 800), **not** in `panel_definition.md`, and not under Amendments 8 or 10.
- The seven-cohort vector is defined in code at `05_clinical.R:47`.

The attribution was wrong even though the seven-cohort set itself is right — the
seventh cohort is registered in a different document.

**Fixed.** `METHODS_FACTS.md` §5 cites `analysis_plan.md` for the endpoint table and
`05_clinical.R` for the cohort vector, and does not attribute PAAD to
`panel_definition.md`. The underlying substantive point — identity checked by count
rather than name — remains a live robustness item for anyone touching the cohort
loop.

**Does it touch a reported number?** No. All seven cohorts are correctly present in
every output; only the citation was misplaced.

---

## F4 — `2c22af6f` · warn/low · INCORRECT — refuted by execution

**Claims.** That my assertion "the guarded driver runs at Rscript top level
(`sys.nframe() == 0` → `DRIVER RAN`), verified by execution" is unsupported,
because the confirming stdout does not appear in the traceable transcript — only
the `source()` half is visible. The finding concedes the claim is ordinary R
semantics and calls it a presentation gap, not a wrong conclusion.

**Correct?** **No — refuted by running it.** Both halves, just now:

```
--- executed via Rscript:      DRIVER RAN (sys.nframe()==0 at Rscript top level)
--- the same file source()d:   driver skipped
```

The behaviour is exactly as claimed: the driver fires under `Rscript` and is
skipped under `sys.source()`. This is the pattern every script `05`–`14` relies on,
and it is independently confirmed by the fact that **all of them ran and produced
their committed outputs** — impossible if the guard suppressed the driver under
`Rscript`.

The finding is right that the confirming line was missing from the persisted
transcript, and that is a fair observation about evidence. But the claim itself is
correct, so the disposition is INCORRECT rather than CORRECT-BUT-UNACTED: the
conclusion it questions holds.

**Does it touch a reported number?** No.

---

## F5 — `d93d3361` · warn/medium · CORRECT BUT NOT ACTED ON

**Claims.** A final message described a saved audit memo as containing *"20 clean
checks"*, but the artifact's `## CHECKED AND CLEAN` section lists only **16**. Four
items present in the submitted JSON never reached the saved file. The finding notes
the companion count — 17 findings — **is** correct, and that the substantive verdict
and the two blocking findings are unaffected.

**Correct?** **Yes.** Verified by reading the artifact itself
(`audit9_09_null_findings.md`, version `f944b678`): its `CHECKED AND CLEAN` section
contains **exactly 16 bullets**, programmatically counted. So "20 clean checks"
overstated the saved file by four.

**Why it is not acted on.** The artifact is a **sub-agent audit memo, never
committed to this repository** — `git log --all -- "*audit9_09_null_findings*"`
returns nothing. It is not a deliverable, it is not cited by any manuscript
document, and the four missing items are all *clean* checks, so none of them
records a defect that went unaddressed. Editing a historical audit artifact to
match a prose claim about it would be rewriting the record rather than correcting
it.

**What matters for the manuscript**: the two blocking findings from that audit (B1,
B2) were both fixed and are consistent across the JSON, the artifact and the
execution log; the 17-findings count was right. **Nothing in `09_null.R`'s
committed outputs is affected.**

**Does it touch a reported number?** No.

---

## What a reader should take from this

1. **Nothing here changes a committed number.** F2 asked whether a number had been checked; it had not, it now has, and it is correct.
2. **Two findings (F1, F5) concern evidence quality on claims that never entered a deliverable** — a redacted search figure, and a count of clean checks in an uncommitted audit memo.
3. **F3 was a real citation error** and is fixed.
4. **F4 is wrong on its conclusion**, and I refuted it by execution rather than argument.
5. **No fail-level finding is open.** All five are warn-level.

The one item worth carrying into the manuscript as a genuine caveat is **F1's
downstream connection**: the broad reading of rule 3.3 selects the 140-gene primary
list, and the biological rationale for that reading is committed in mechanistic
form only. If the manuscript wants to state a magnitude of male bias, it needs a
citable source — not the redacted search result.
