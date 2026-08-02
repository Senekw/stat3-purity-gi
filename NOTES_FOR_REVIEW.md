# Notes for review

Observations surfaced but NOT acted on, per the scope discipline rule.

## 21. RESOLVED — Amendment 15 settles rule 3.3 on the broad reading; primary list is 140

**Supersedes §17, which recorded this as open.** Amendment 15 (2026-08-02) applies
rule 3.3 on its literal text: any panel gene annotated on chrX or chrY is
excluded, PAR genes included. **Primary scoring set = 140**; the 143-gene narrow
list is retained as a prespecified sensitivity set and both are written to disk.

### What I verified rather than accepted

Amendment 15's reasoning turns on the four PAR genes spanning *both* PAR regions,
so a single dosage argument cannot cover them. I recomputed the assignment from
the annotation's own chrX coordinates against the GRCh38 PAR intervals rather
than relying on the amendment or on memory:

| Gene | chrX coordinates | Region |
|---|---|---|
| CRLF2 | 1,187,549–1,212,723 | **PAR1** |
| CSF2RA | 1,268,800–1,310,381 | **PAR1** |
| IL3RA | 1,336,616–1,382,689 | **PAR1** |
| IL9R | 155,997,581–156,010,817 | **PAR2** |

(GRCh38 PAR1 chrX:10,001–2,781,479; PAR2 chrX:155,701,383–156,030,895. Each gene
also carries a `_PAR_Y` row at the corresponding chrY coordinate.) **Three PAR1,
one PAR2 — the mixture Amendment 15 describes.** The lock script now recomputes
this at run time and halts if the assignment or the assembly changes.

### Both of the earlier justifications were unsound, and neither is used

Amendment 15 records this explicitly, and it is worth restating because the
sequence is the point:

1. **Mine** (for retaining PAR genes): *"they escape X-inactivation and are
   present in two copies in both sexes."* Asserted from memory, unchecked,
   written into the lock script.
2. **The audit's replacement**: *"XCI escape is the mechanism that produces
   female-biased expression."* Correct for non-PAR X-linked escapees — but not
   for PAR1, which has a functional Y homolog.

So the objection that corrected me was itself not uniformly right, and the
resolution rests on **neither**: Amendment 15 decides on the registered text,
which is categorical and names no exception. Defining a scoring set on a
contested empirical claim about four genes would have been the worse failure.

### Direction of bias — now conservative

The three additional exclusions (CSF2RA, IL3RA, IL9R) are all non-epithelial.
Retaining them would make the score a *more stromal* readout, the direction
consistent with this study's thesis; excluding them makes the hypothesis **harder**
to support. The reading that was implemented before (narrow/143) was the
hypothesis-friendly one. This reverses that.

### One correction to the arithmetic as commonly stated

The three rules are **not disjoint** under the broad reading. Rule 2 catches 6
genes and rule 3.3 catches 7, but **CRLF2 is caught by both**, so the union is
**12**, not 13, and the primary list is 140. The retained and excluded sets do
partition the 152 exactly (140 + 12), and `reason()` records both rules for
CRLF2 rather than letting the first one silence the second. Under the earlier
narrow reading the rules *were* disjoint — that property was specific to it.

k is untouched (Amendment 14: computed over the locked 152), and the Amendment 3
branch is BUILT under 152, 143 and 140 alike.

## 29. THE 09 RUN ORDER MISSTATED THREE B.m PARAMETERS; THE REGISTERED TEXT WAS FOLLOWED

The Part B run order of 2026-08-02 specified three B.m parameters from memory that
do not match the registered text. The author confirmed the same day that the
registered text governs and the run order was NOT to be followed. A first draft of
09 had been written to the run order and was discarded; its audit was cancelled
mid-flight.

| Parameter | Run order said | B.m registers | Used |
|---|---|---|---|
| N | ~1,000 | 10,000 per cohort | **10,000** |
| Set size | 140 | 152, or the cohort's realised size | **152** |
| Matching | mean expression x stromal correlation | decile of mean log2 expression x decile of expression variance, within cohort | **expression x variance** |

**The realised-size clause is inert here, verified:** all 152 panel genes survive
the annotation and zero-variance filters in all seven cohorts (152/152 each), so
the realised size IS 152 everywhere. That clause concerns genes missing from a
cohort, which is a different thing from the 140-gene post-exclusion list; the two
are not conflated.

Three configurations were run, each with an explicit declared role so that two
full-N p-values could not arrive with equal standing:

- `registered_152|registered` — **primary_registered**, N = 10,000
- `tested_140|registered` — size_sensitivity, N = 10,000 (the signature 08 tested)
- `tested_140|EXPLORATORY_POSTHOC` — exploratory_posthoc, N = 1,000, additionally
  matched on stromal correlation. **Added after the primary result was seen**, at
  the author's request. Labelled EXPLORATORY_POSTHOC in every output row.

## 30. THE NULL COMPARISON WAS NOT LIKE-FOR-LIKE UNTIL CORRECTED

Found by the Implementation Auditor before 09 ran. 08's pooled attenuation_total
(-0.024766) was meta-analysed with a **paired bootstrap SE**; every null signature
is pooled with **sqrt(se2^2 + se4^2)**, since no per-signature bootstrap is run.
Those are different estimators: the naive form ignores the positive covariance
between b2 and b4 and is **1.72x larger** on this data (mean 0.1835 vs 0.1069).
Different weights give a different pooled POINT estimate, not merely a different
interval, so comparing 08's value against a naive-pooled null distribution would
not have been like-for-like.

**Corrected:** p_atten is computed against the observed value re-pooled with the
nulls' own estimator, **-0.022577**. 08's -0.024766 remains the reported estimate
and is carried alongside. The difference is 0.0022 in log-HR units — small, but it
is the decisive analysis and the two numbers are not interchangeable.

## 31. AN AUDIT FINDING I CHECKED AND REJECTED — the matching gate is not vacuous

The auditor reported (R9) that the tolerance gate could not fail, simulating a
completely unmatched draw at median 0.103 SD against a 0.25 SD threshold. My own
measurement on COAD contradicts this: **unmatched draws give median 1.62 SD**
(p95 1.74, max 1.84), far above the gate, which therefore fires correctly.

The underlying point was still worth acting on. A set-level difference of means
over 152 genes is shrunk by averaging, so a **paired per-gene** statistic is more
discriminating: matched draws give 0.32 SD against 1.84 SD unmatched, and a second
gate at PAIRED_TOL = 0.75 SD now tests it. Both thresholds are unregistered and
labelled as such in the source. The paired statistic immediately earned its place
— it shows the exploratory scheme cuts stromal mismatch from 1.52 SD to 0.25 SD,
a difference the set-level statistic could not see.

## 26. PARAMETER COUNTS KEPT FIXED, AGAINST AN AUDIT RECOMMENDATION

The Implementation Auditor proposed deriving each model's EPV denominator from
the realised fit (`sum(!is.na(coef(fit)))`) so a cohort whose `missing` stage
level is empty would count one fewer parameter. **Not adopted.**
analysis_plan.md:1186-1187 fixes the counts: "Parameter counts are fixed in B.j
(M1 = 1, M2 = 5, M3 = 6, M4 = 7, with `stage_group` contributing two parameters
for its three levels)." Deriving them would deviate from the registered
denominator and make the EPV band depend on realised level occupancy — a cohort
could cross a band boundary because one stage cell happened to be empty. The one
adjustment the plan itself specifies (dropping `sex` below 10 of either sex) is
applied, and only where a sex term exists.

Realised: all seven cohorts have use_sex = TRUE, so the adjustment is dormant.

## 27. BOOTSTRAP CONVERGENCE WARNINGS TRACED TO A NUISANCE TERM, NOT THE SCORE

The 08 run emitted 50+ "Loglik converged before variable 5; coefficient may be
infinite" warnings. Traced rather than suppressed: variable 5 is
`stage_groupmissing`, and the source is READ, whose `missing` stage level holds 9
patients with 1 event. In 300 test resamples **zero** had any |coefficient| > 20
and the score coefficient stayed within [-0.79, 0.90].

`boot_nuisance_unstable` is now written to `attenuation_per_cohort.csv`, counting
resamples where any coefficient diverged or aliased: READ 1/2000, PAAD 93/2000,
all others 0. The score coefficient is screened independently, and 0 of 12,000
resamples across the six cohorts failed. The warnings reflect coxph's own
threshold, which is more sensitive than the divergence screen; the estimand is
unaffected.

## 28. THE ALTERNATIVE-ENDPOINT SENSITIVITY SPLITS ITS COHORT SET

In the PRIMARY analysis the M2 and M4 pools rest on the same six cohorts, so
Amendment 8 item 4's matched re-pool is not needed. **It IS needed inside the
alternative-endpoint sensitivity:** switching READ from PFI to OS gives 25 events,
EPV 3.57 at M4, so the EPV < 5 rule drops READ from that model — M2 pools 6
cohorts and M4 pools 5.

Reported in `sensitivity_matched_repool.csv`, both models re-pooled over the
intersection (COAD+STAD+ESCA+PAAD+LIHC, n = 1,590): M2 = 0.0826 (95% CI -0.0066
to 0.1718), M4 = 0.1013 (-0.0337 to 0.2364). Without this the variant's apparent
M2-to-M4 change would partly reflect READ leaving, not adjustment.

## 24. RESOLVED — the instructed `final_140` label was erroneous; the artefact is correct

The Part B run order asked to assert `list_id == "final_140"`. **The instructing
label was written from memory and is wrong**; `04_lock_gene_list.R` wrote
`primary_140` (and `sensitivity_143`), and the locked artefact committed at
744d84c is correct and must not change. Confirmed by the study author 2026-08-02.

No gene was ever in question. The real identity check is the **md5 of the sorted
gene symbols**, asserted on every run:

| List | list_id (canonical) | n | md5 of sorted symbols |
|---|---|---|---|
| Primary | `primary_140` | 140 | `eb167e8c7a33b4202bd609a17defa629` |
| Sensitivity | `sensitivity_143` | 143 | `8a54835eedcfe3b23dcb56ce74805a87` |

plus a check that the 140 is a subset of the 143. A file with the right row count
and label but the wrong genes fails the digest; that is what closes stop condition
6, not the label string.

`07_score.R` retains a one-line note recording that the run order's label differed
from the artefact's, downgraded from a warning since there is nothing to decide.

## 23. STROMAL FALLBACK TRIGGERED IN ALL SEVEN COHORTS — and it moots an open question

B.h prespecifies that if the panel-derived stromal subscore is collinear with the
main score (|r| > 0.9), the ESTIMATE stromal score is substituted. Realised
correlations: COAD 0.982, READ 0.981, STAD 0.965, ESCA 0.962, PAAD 0.978, LIHC
0.986, CHOL 0.957. The threshold is exceeded in **every** cohort, so `stromal_score`
in model M4 is the ESTIMATE StromalScore throughout. The plan anticipated this
("model 4 is collinear with model 1 by construction... stated now so the choice is
not made after seeing k"). The panel-derived subscore is retained as
`stromal_score_subscore` for inspection.

**This resolves an audit finding by making it moot.** The Implementation Auditor
objected (F2, blocking) that the subscore was built over the locked 152-gene panel
and so readmitted genes that exclusion rules 2 and 3.3 had removed. Realised, after
the evaluability fix, **six** are readmitted: TIMP1, IL3RA, IL2RG, CSF2RA, IL9R and
OPRM1. **Five of the six are excluded only under rule 3.3** — TIMP1 (chrX, detected
in 98.4% of stromal cells), IL3RA 61.1%, IL2RG 56.7%, CSF2RA 44.1%, IL9R 4.2% —
with OPRM1 the one excluded as undetectable. (An earlier version of this note said
four sex-chromosome-only genes; IL9R was omitted from that count.) The registered text
(analysis_plan.md:718-719) says "panel genes", which is the 152; the auditor's
reading would use the 140. Both were computed rather than one being chosen:

| Domain | Stromal set size | r with main score |
|---|---|---|
| Panel 152 (registered wording) | 103 | 0.956-0.986 |
| Final 140 (audit reading) | 97 | 0.965-0.987 |

The two variants correlate at r >= 0.997 with each other and BOTH exceed 0.9 in
every cohort, so both trigger the same fallback and neither reaches the model. The
question does not need deciding for Part B. **It would need deciding if the
fallback were ever not triggered** — e.g. in an external validation cohort — so the
alternative remains computed as `stromal_score_140` and the readmitted genes are
named in `score_gene_sets.txt`.

Separately, 3 panel genes (CRLF2, DNTT, LEP) are evaluable in fewer than 2 atlases
and so can neither satisfy nor refute the two-of-three dominance rule. They are now
excluded from BOTH the dominant and the stromal sets rather than defaulting to
non-dominant, which is why the panel-domain subscore is 103 and not 106.

## 22. PURITY SOURCE SPLITS 3/4 ACROSS COHORTS — the switch rule is doing real work

`aran_purity.xlsx` was absent at the start of Part B. The plan (3.4) anticipated
this and prespecified a switch rule, but letting a MISSING DOWNLOAD flow through
that rule would have inverted B.j's registered CPE-primary ordering for all seven
cohorts on an acquisition gap rather than on the data. The file was therefore
obtained (Aran et al. 2015, Nat Commun 6:8971, Supplementary Data 1; md5
c459e6a965789b96860fc77bd346c681, 9,364 rows, 21 cancer types) and the rule
applied to genuine coverage. **06 now HALTS if the file is absent** rather than
silently taking the ESTIMATE branch.

Realised coverage over the analysis set:

| Cohort | CPE coverage | Primary source |
|---|---|---|
| COAD | 99.8% | CPE |
| READ | 100% | CPE |
| LIHC | 100% | CPE |
| STAD | 0% | ESTIMATE |
| ESCA | 0% | ESTIMATE |
| PAAD | 0% | ESTIMATE |
| CHOL | 0% | ESTIMATE |

The four zeroes are **not a join failure**: STAD, ESCA, PAAD and CHOL are absent
from the Aran 2015 freeze entirely, verified by listing the table's own 21 cancer
types. The switch rule applies exactly as registered.

**Consequence for the meta-analysis, for the reviewer to weigh:** the purity
covariate is not one quantity across the seven cohorts. Three use a consensus of
four orthogonal methods; four use an ESTIMATE-derived value whose conversion,
cos(0.6049872018 + 0.0001467884 * S), was calibrated on Affymetrix arrays and is
applied here to Illumina RNA-seq TPM. `purity_calibration` in
`purity_summary.csv` records which is which per cohort. Where both exist they
agree at r = 0.70 (COAD), 0.66 (READ), 0.74 (LIHC) — the registered cross-method
sensitivity, computable only in those three.

## 21. ESTIMATE PURITY CONVERSION FOLDS BACK BELOW ESTIMATEScore = -4121.5

Found by the Implementation Auditor on 06 before it ran; confirmed
mathematically. cos(a + bS) is monotone-decreasing in S only while the angle
a + bS is in [0, pi/2], i.e. S in [-4121.5, +6579.6]. Below -4121.5 the angle
goes negative and the cosine **folds**: purity decreases as the tumour gets purer,
and two distinct scores map to one purity — S = -6000 and S = -2243 both give
0.962223, verified numerically.

The original guard tested the COSINE for [0,1], which cannot catch this because
folded values are in range. The guard now tests the ANGLE.

**Latent, not active, on this data:** observed ESTIMATEScore spans about -2689 to
+4891 across the cohorts checked, so no sample sits in the fold region and no
reported number changes. It is recorded because it would have become active
silently in a validation cohort with more stroma-rich samples, and because the
affected quantity is the primary model's adjustment covariate.

## 20. AMENDMENT 14's STATED REASON IS PARTLY FALSE — conclusion holds on other grounds

**The decision Amendment 14 makes is sound. One clause of its justification is
factually wrong, and it is the clause carrying the argument.** Registered text is
not edited; the correction is recorded here.

Amendment 14's reason paragraph says the 143-gene final list *"is derived from
the compartment output and did not exist when k was defined."* The **second
clause is false**; the first is broadly right but overstated:

| Claim | Check | Verdict |
|---|---|---|
| "derived from the compartment output" | **6 of the 9** exclusions come from rule 2, which reads `detection_rate_by_compartment.csv` — Part A output. The other **3** (IL13RA1, IL2RG, TIMP1) come from rule 3.3, TCGA annotation only, independent of Part A. | **Overstated, not false.** A majority of exclusions *are* compartment-derived, so the clause is broadly right; it is only imprecise in implying the list is *wholly* so. This is the weaker of the two objections. |
| "did not exist when k was defined" | Section 3's exclusion rules are present in the **first commit**, `361d287` ("Starting state"), verbatim as applied. | **False.** The rules predate k's definition; only their *evaluation* is recent. |

**Why the conclusion survives anyway.** The literal reading holds on its own:
Amendment 3 defines k over *"panel genes"*, and the panel is the 152 locked at
`ac9c5e0`. k is a property of the panel by definition; the 143-gene list is the
**scoring set**, a different object serving a different purpose in Part B. That
argument needs neither the false chronology clause nor the overstated one.

**Why it matters even though the answer is unchanged.** The chronology clause is
the load-bearing one, and it is false: it makes the decision look forced by
sequence when it is actually a definitional choice
— and a definitional choice deserves the scrutiny Amendment 14's own disclosure
paragraph invites. The disclosure itself is unusually good and stands: it states
the choice is not blind, that both values were computed first, that the option
*not* taken (43) is the hypothesis-friendly one, and it names the three genes
driving the difference.

**If an Amendment 15 is written** (see §17 on the PAR question), it should restate
the reason on the definitional ground alone and drop the "derived from the
compartment output / did not exist" clauses.

## 19. UNEQUAL EVIDENCE BASE — half the rule-2 exclusions rest on one atlas

**Surfaced by the `n_atlases_evaluable` column added after the lock audit. It
does not change the locked list, but it changes how confidently three of the six
exclusions can be stated.**

Rule 2 excludes a gene whose maximum detection across all atlases is under 1%.
The number of atlases contributing to that maximum is **not equal across genes**:

| Gene | Max detection | Atlases contributing | Where the max was |
|---|---|---|---|
| CRLF2 | 0.1385% | **1** | GSE178341 / myeloid |
| DNTT | 0.0270% | **1** | GSE178341 / lymphoid |
| LEP | 0.3032% | **1** | GSE178341 / myeloid |
| OPRM1 | 0.8720% | 2 | GSE178341 / lymphoid |
| PAX3 | 0.1356% | 2 | Peng / epithelial |
| GFAP | 0.5017% | 3 | Peng / epithelial |

For contrast, **138 of the 143 retained genes were evaluated against all three
atlases** (the other 5 against two). A gene tested on one atlas has fewer chances
to clear 1% than one tested on three, so **the exclusion threshold is not applied
with equal power across the panel** — three genes were dropped on strictly less
evidence than the rule's design assumes.

**Cause.** Two mechanisms, both already registered: Amendment 12's intersection
removes CRLF2 and LEP from GSE125449 entirely, and the `evidence_ok` floor (≥20
summed counts) suppresses others per atlas. Neither is a defect; the interaction
with rule 2 was simply never stated.

**How much it could matter.** For CRLF2 and LEP the missing atlas is GSE125449 —
*not* where their maximum came from — and the gene is absent from that deposit
altogether, so no rescue value can be computed even in principle. OPRM1 at
0.872% is the one to watch: it is the closest of any excluded gene to the 1%
boundary and was judged on two atlases, so a third might plausibly have cleared
it.

**Not acted on.** Changing rule 2 to require equal atlas coverage, or to
condition on the number available, would be a post-hoc redefinition of a
registered rule after seeing which genes it catches — precisely what
preregistration forbids. The rule is applied as written. This is recorded so the
paper can state the exclusions with the right confidence, and so a reader can see
which of the six rest on one atlas. `n_atlases_evaluable` and `max_at` are now
columns in both `final_gene_list.csv` and `final_gene_list_exclusions.csv`.

## 18. CORRECTION — the rule-2 genes are not "MSigDB pathway membership" entrants

The instruction accompanying Amendment 14 asked me to log that GFAP, PAX3 and
OPRM1 "are neural-lineage genes entering via MSigDB pathway membership." The
first half is right and worth reporting; **the second half is contradicted by the
panel's own provenance**, so I logged the verified route instead.

All three entered by **criterion B (direct transcriptional target)**, not
criterion A (MSigDB `HALLMARK_IL6_JAK_STAT3_SIGNALING`):

| Gene | `route` | In HALLMARK set | STAT3 term — library it came from | TRRUST v2 (PMID) |
|---|---|---|---|---|
| GFAP | `B_only` | **no** | `STAT3_HELAS3_HG19` — **ENCODE_ChIP-seq.gmt**; `STAT3` — **ReMap_ChIP-seq.gmt** | Activation (21833841), Unknown (11740937) |
| PAX3 | `B_only` | **no** | `STAT3_23295773_CHIPSEQ_U87_HUMAN`, `STAT3_24763339_CHIPSEQ_IMNESCS_MOUSE` — both **Literature_ChIP-seq.gmt** | Repression (19074888) |
| OPRM1 | `B_only` | **no** | `STAT3_23295773_CHIPSEQ_U87_HUMAN` — **Literature_ChIP-seq.gmt** | Activation (15448191) |

Every identifier above is reproducible from committed files: the GMT libraries in
`data/panel/`, `data/panel/trrust_rawdata.human.tsv`, and the `ev_*` flag columns
of `panel_locked.csv`. An earlier revision of this table listed GFAP's terms
without naming their libraries and grouped all three genes' evidence together,
which implied they shared a source; they do not — GFAP's evidence is ENCODE +
ReMap, while PAX3's and OPRM1's is the Literature library.

The two rule-2 genes that *did* enter via criterion A are **DNTT** and **CRLF2**
(both `A_only`, both in the HALLMARK set) — neither is neural.

This sharpens rather than weakens the observation. Three of the six undetectable
genes are neural-lineage, and they entered on **direct STAT3 binding evidence**:
PAX3 and OPRM1 both trace to `STAT3_23295773_CHIPSEQ_U87_HUMAN` — a **U87
glioblastoma** experiment — and GFAP to a HeLa STAT3 ChIP. So the neural
signal enters through the cell lines behind the ChIP-seq libraries, not through
pathway annotation. A panel built from ChIP-seq in glioma and HeLa lines carries
their lineage context into a GI study, and these genes then prove undetectable in
GI tumour tissue. That is a cleaner statement of the deposit-provenance caution
than "pathway membership" would have been, and it belongs in the discussion.

## 17. UNREGISTERED INTERPRETIVE CHOICE — rule 3.3 and the pseudoautosomal genes

> **SUPERSEDED by §21.** Amendment 15 (2026-08-02) decided this on the broad
> reading: primary scoring set 140, with 143 retained as a prespecified
> sensitivity set. This section is kept as the record of how the question was
> found and what was wrong with my first answer to it — not as current guidance.

**Surfaced while locking the list. It does not change any reported Part A number,
but it does change the scoring set, so it is recorded rather than absorbed.**

Rule 3.3 reads, in full: *"Sex-chromosome genes, since cohorts differ in sex
composition."* Applying it needs a decision the registered text does not make.

Four panel genes lie in the **pseudoautosomal region** and are annotated on
**both** chrX and chrY — the annotation emits a second `_PAR_Y` row for each:
**CRLF2, CSF2RA, IL3RA, IL9R**. (These are the same four `_PAR_Y` duplicates that
appeared in Part A's GSE178341 reader, where all four Y-copies carried zero
counts.) Two readings are available:

| Reading | Basis | Excludes | Final list |
|---|---|---|---|
| **NARROW** (implemented, **provisional**) | — see below | genes annotated *only* on X or Y | **143** |
| BROAD | the rule's stated **text** | anything with any row on X or Y | 140 |

### My first justification for NARROW was wrong, and is withdrawn

An earlier revision of this section, and of the lock script, justified narrow on
the biology: that PAR genes *"escape X-inactivation and are present in two copies
in both sexes, so cohort sex composition does not bias them."* **That argument is
wrong at the operative step.** The Implementation Auditor caught it.

Escaping X-inactivation is precisely the mechanism that **produces**
sex-differential expression: an escapee is transcribed from *both* X copies in XX
individuals, so it tends to be expressed **higher in females**. Rule 3.3's stated
reason — *"since cohorts differ in sex composition"* — therefore argues **for**
excluding PAR genes, not for retaining them. I asserted the claim from memory
without checking it, and it went into a lock script that defines the Part B
scoring set. That is the error worth recording here, more than the conclusion.

### So the question is open, and the implemented reading is provisional

Narrow is implemented **only because Amendment 14's registered text describes a
143-gene final list**, and 143 is what narrow produces. But Amendment 14 was
written *before* the PAR question was identified: it records the number without
deciding the question. A number is not an argument, and I no longer have a
biological one. `PAR_UNDECIDED <- TRUE` in `04_lock_gene_list.R` marks this.

**DIRECTION OF BIAS** — absent from my first version, and it has one. Narrow
retains CSF2RA, IL3RA and IL9R in the scoring set. None is epithelial-dominant,
so each adds a **non-epithelial** gene to the score — the direction that favours
this study's "substantially stromal" thesis. **The reading I implemented is the
hypothesis-friendly one.** That is a reason for you to decide it explicitly.

**What does not depend on it**: k is 43 under both readings, and the Amendment 3
branch is BUILT under 152, 143 and 140 alike (30.3%, 30.1%, 30.7% of panel, all
above the 20% floor). No Part A quantity moves.

**What does**: the Part B scoring denominator, 143 vs 140.

**Recommended resolution** — decide it in a dated Amendment 15 on the merits,
with a direction-of-bias statement, and have the script cite that rather than
Amendment 14. If it lands on broad/140, Amendment 14's "143-gene final list"
sensitivity figure needs restating. A defensible route to narrow that does not
rest on the withdrawn claim: measure sex-differential expression of these four
genes in the seven cohorts directly. That requires opening expression data, which
the Part A/B boundary currently forbids, so it belongs after the boundary or in
an explicit annotation-only carve-out.

`04_lock_gene_list.R` asserts the PAR set is exactly those four genes **and**
that each is still annotated on both sex chromosomes, so an annotation change
halts the lock rather than silently relocking a different scoring set.

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
