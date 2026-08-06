# RESULTS_FACTS — values for the manuscript results section

Every value is **transcribed** from the committed file named in its section
heading. **No statistic was computed for this document.** Two places derive a
*list* by filtering a committed table (the 24 all-three-atlas genes, §1.5); both
are labelled, and the derivation reproduces the committed count.

Assembled 2026-08-06 at commit `6f08bb9`. Numbers are shown at the
precision stored in the file, so a value written `0.1212` is 4 dp in the source, not
rounded here.

---

## 1. Part A — compartment attribution

### 1.1 k and variants — `output/k_estimates.csv`

| Quantity | Estimate | 95% CI (bootstrap percentile) | Branch |
|---|---|---|---|
| **k** | **46** | 38–57 | epithelial subscore BUILT; subscore survival models run as SECONDARY |
| **k_all3** | **24** | 15–31 | — |
| **k_evalall** | **43** | 35–55 | — |
| **k_50** | **96** | 86–103 | epithelial subscore BUILT; subscore survival models run as SECONDARY |

**Branch as realised:** `epithelial subscore BUILT; subscore survival models run as SECONDARY`.

Definitions (from `03_compartments.R`): **k** = dominant in ≥ 2 atlases;
**k_all3** = dominant in 3 and evaluable in all 3; **k_evalall** = dominant in ≥ 2
and evaluable in all 3; **k_50** = dominant at the midpoint π = 0.50 only.
Runtime assertion `k_all3 ≤ k_evalall ≤ k` held (24 ≤ 43 ≤ 46).

### 1.2 Per-atlas dominance counts — `output/per_tissue_dominance.csv`

| Atlas | Tissue | Genes dominant |
|---|---|---|
| GSE125449 | liver_biliary | **43** |
| GSE178341 | colorectal | **67** |
| Peng | pancreatic | **38** |

Sum of per-atlas counts = **148** gene-atlas dominance calls. This is *not* k:
k counts genes dominant in ≥ 2 atlases (**46**), so the two figures answer
different questions and should not be presented as if one checks the other.

### 1.3 Evaluability distribution — `output/evaluability_distribution.csv`

| Atlases in which the gene is evaluable | Genes |
|---|---|
| 3 | **140** |
| 2 | **9** |
| 1 | **3** |
| 0 | **0** |
| **total** | **152** |

152 = the 152-gene inferential set. Genes evaluable in 0 atlases contribute no
dominance call and are `NA`, never `FALSE`.
### 1.4 The 24 genes dominant in all three atlases

**Derived** by filtering `output/compartment_dominance_matrix.csv` to panel genes
that are evaluable in all three atlases and dominant in all three. The count
reproduces the committed `k_all3` = **24** (MATCH), which is what
validates the filter. The gene *names* are not themselves in a committed file.

`ACVR1B`, `BCL2L1`, `CXCL1`, `DDIT3`, `FGG`, `FGL1`, `GFAP`, `IL17RB`, `KRT17`, `LGALS3BP`, `LTBR`, `MMP7`, `MUC1`, `PF4`, `PHB`, `PPARA`, `REG1A`, `S100A11`, `SHH`, `STAM2`, `TNFRSF10B`, `TNFRSF12A`, `TNFRSF21`, `VEGFA`

### 1.5 Origin six and the three non-qualifying — `output/origin_six_compartment.csv`

f(π) is the epithelial fraction at that purity; **dominant** requires f > 0.50 at
**every** one of the 41 grid points, not at the values shown.

| Gene | Atlas | Tissue | f(0.30) | f(0.50) | f(0.70) | Evaluable | Dominant | Dom at 0.50 | Qualifying |
|---|---|---|---|---|---|---|---|---|---|
| **SOCS3** | GSE125449 | liver_biliary | 0.0571 | 0.1237 | 0.2478 | TRUE | **FALSE** | FALSE | TRUE |
| **SOCS3** | GSE178341 | colorectal | 0.1089 | 0.2219 | 0.3995 | TRUE | **FALSE** | FALSE | TRUE |
| **SOCS3** | Peng | pancreatic | 0.1295 | 0.2577 | 0.4475 | TRUE | **FALSE** | FALSE | TRUE |
| **MYC** | GSE125449 | liver_biliary | 0.3274 | 0.5318 | 0.7261 | TRUE | **FALSE** | TRUE | TRUE |
| **MYC** | GSE178341 | colorectal | 0.8826 | 0.9461 | 0.9761 | TRUE | **TRUE** | TRUE | TRUE |
| **MYC** | Peng | pancreatic | 0.1847 | 0.3459 | 0.5523 | TRUE | **FALSE** | FALSE | TRUE |
| **IL6** | GSE125449 | liver_biliary | 0.0182 | 0.0415 | 0.0918 | TRUE | **FALSE** | FALSE | TRUE |
| **IL6** | GSE178341 | colorectal | 0.0119 | 0.0273 | 0.0615 | TRUE | **FALSE** | FALSE | TRUE |
| **IL6** | Peng | pancreatic | 0.0381 | 0.0846 | 0.1775 | TRUE | **FALSE** | FALSE | TRUE |
| **BCL2** | GSE125449 | liver_biliary | 0.2185 | 0.3948 | 0.6035 | TRUE | **FALSE** | FALSE | FALSE |
| **BCL2** | GSE178341 | colorectal | 0.1838 | 0.3445 | 0.5509 | TRUE | **FALSE** | FALSE | FALSE |
| **BCL2** | Peng | pancreatic | 0.2014 | 0.3704 | 0.5785 | TRUE | **FALSE** | FALSE | FALSE |
| **MMP9** | GSE125449 | liver_biliary | 0.0871 | 0.1820 | 0.3418 | TRUE | **FALSE** | FALSE | FALSE |
| **MMP9** | GSE178341 | colorectal | 0.0142 | 0.0326 | 0.0729 | TRUE | **FALSE** | FALSE | FALSE |
| **MMP9** | Peng | pancreatic | 0.1202 | 0.2418 | 0.4266 | TRUE | **FALSE** | FALSE | FALSE |
| **HGF** | GSE125449 | liver_biliary | 0.1044 | 0.2139 | 0.3884 | TRUE | **FALSE** | FALSE | FALSE |
| **HGF** | GSE178341 | colorectal | 0.0107 | 0.0246 | 0.0556 | TRUE | **FALSE** | FALSE | FALSE |
| **HGF** | Peng | pancreatic | 0.0405 | 0.0896 | 0.1867 | TRUE | **FALSE** | FALSE | FALSE |

**Non-qualifying reasons** (criterion B, `non_qualifying_reason`):

- **BCL2**: no ChIP-seq evidence in any group (fails Criterion B)
- **MMP9**: mouse ChIP-seq only (excluded by Amendment 2, human-only)
- **HGF**: mouse ChIP-seq only (excluded by Amendment 2, human-only)

⚠ **A presentational trap, recorded because it is visible in the figure.** Nine of
the 54 f values exceed 0.50 somewhere, but only **MYC in GSE178341** exceeds it at
all three shown grid points and is the only atlas-gene pair with `dominant = TRUE`.
**BCL2 clears 0.50 in all three atlases at π = 0.70** while being non-qualifying and
non-dominant. A reader looking only at a 50% reference line will mis-read this;
`fig3` marks the file's own `dominant` call for that reason.

---

## 2. Survival models, per cohort — `output/survival_per_cohort.csv`

`score` is per 1 SD within cohort. `p_adj_BH` is Benjamini–Hochberg **within
model**, family = the six meta-eligible cohorts. CHOL carries no adjusted p (it is
outside the family, Amendment 8) and is reported descriptively in §2.2.


### 2.1.1 M1

| Cohort | n | Events | Params | EPV | Band | β | SE | HR | 95% CI | p | p (BH) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD | 455 | 102 | 1 | 102.00 | fit_and_pool | 0.005159 | 0.094338 | **1.0052** | 0.8355 to 1.2093 | 0.9564 | 0.9564 |
| READ | 165 | 36 | 1 | 36.00 | fit_and_pool | -0.035099 | 0.166886 | **0.9655** | 0.6961 to 1.3391 | 0.8334 | 0.9564 |
| STAD | 403 | 157 | 1 | 157.00 | fit_and_pool | 0.113122 | 0.083450 | **1.1198** | 0.9508 to 1.3187 | 0.1752 | 0.3504 |
| ESCA | 184 | 77 | 1 | 77.00 | fit_and_pool | 0.054184 | 0.115658 | **1.0557** | 0.8416 to 1.3243 | 0.6394 | 0.9564 |
| PAAD | 178 | 93 | 1 | 93.00 | fit_and_pool | 0.276873 | 0.117908 | **1.3190** | 1.0468 to 1.6619 | 0.0189 | 0.0566 |
| LIHC | 370 | 130 | 1 | 130.00 | fit_and_pool | 0.233098 | 0.088456 | **1.2625** | 1.0615 to 1.5015 | 0.0084 | 0.0505 |
| CHOL | 35 | 18 | 1 | 18.00 | fit_and_pool | -0.160616 | 0.216782 | **0.8516** | 0.5568 to 1.3025 | 0.4587 | — |

### 2.1.2 M2

| Cohort | n | Events | Params | EPV | Band | β | SE | HR | 95% CI | p | p (BH) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD | 455 | 102 | 5 | 20.40 | fit_and_pool | 0.004596 | 0.098124 | **1.0046** | 0.8288 to 1.2176 | 0.9626 | 0.9626 |
| READ | 165 | 36 | 5 | 7.20 | fit_pool_flag_LOO | -0.045714 | 0.157689 | **0.9553** | 0.7013 to 1.3013 | 0.7719 | 0.9263 |
| STAD | 403 | 157 | 5 | 31.40 | fit_and_pool | 0.075526 | 0.085496 | **1.0785** | 0.9121 to 1.2752 | 0.3770 | 0.5655 |
| ESCA | 184 | 77 | 5 | 15.40 | fit_and_pool | 0.146456 | 0.118195 | **1.1577** | 0.9183 to 1.4595 | 0.2153 | 0.4306 |
| PAAD | 178 | 93 | 5 | 18.60 | fit_and_pool | 0.275399 | 0.121129 | **1.3171** | 1.0387 to 1.6700 | 0.0230 | 0.0690 |
| LIHC | 370 | 130 | 5 | 26.00 | fit_and_pool | 0.233928 | 0.090408 | **1.2636** | 1.0584 to 1.5085 | 0.0097 | 0.0580 |
| CHOL | 35 | 18 | 5 | 3.60 | **do_not_fit** | *not fitted* | | | | | |

### 2.1.3 M3

| Cohort | n | Events | Params | EPV | Band | β | SE | HR | 95% CI | p | p (BH) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD | 455 | 102 | 6 | 17.00 | fit_and_pool | -0.019686 | 0.118769 | **0.9805** | 0.7769 to 1.2375 | 0.8684 | 0.9301 |
| READ | 165 | 36 | 6 | 6.00 | fit_pool_flag_LOO | 0.016465 | 0.187710 | **1.0166** | 0.7037 to 1.4687 | 0.9301 | 0.9301 |
| STAD | 403 | 157 | 6 | 26.17 | fit_and_pool | 0.012236 | 0.092029 | **1.0123** | 0.8452 to 1.2124 | 0.8942 | 0.9301 |
| ESCA | 184 | 77 | 6 | 12.83 | fit_and_pool | 0.211618 | 0.136952 | **1.2357** | 0.9448 to 1.6161 | 0.1223 | 0.2446 |
| PAAD | 178 | 93 | 6 | 15.50 | fit_and_pool | 0.457845 | 0.169187 | **1.5807** | 1.1346 to 2.2022 | 0.0068 | 0.0204 |
| LIHC | 370 | 130 | 6 | 21.67 | fit_and_pool | 0.321546 | 0.104778 | **1.3793** | 1.1232 to 1.6937 | 0.0021 | 0.0129 |
| CHOL | 35 | 18 | 6 | 3.00 | **do_not_fit** | *not fitted* | | | | | |

### 2.1.4 M4

| Cohort | n | Events | Params | EPV | Band | β | SE | HR | 95% CI | p | p (BH) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD | 455 | 102 | 7 | 14.57 | fit_and_pool | -0.071333 | 0.137779 | **0.9312** | 0.7108 to 1.2198 | 0.6046 | 0.9069 |
| READ | 165 | 36 | 7 | 5.14 | fit_pool_flag_LOO | -0.022076 | 0.221829 | **0.9782** | 0.6333 to 1.5109 | 0.9207 | 0.9273 |
| STAD | 403 | 157 | 7 | 22.43 | fit_and_pool | 0.008397 | 0.092021 | **1.0084** | 0.8420 to 1.2077 | 0.9273 | 0.9273 |
| ESCA | 184 | 77 | 7 | 11.00 | fit_and_pool | 0.207848 | 0.139230 | **1.2310** | 0.9370 to 1.6173 | 0.1355 | 0.2710 |
| PAAD | 178 | 93 | 7 | 13.29 | fit_and_pool | 0.411072 | 0.173273 | **1.5084** | 1.0741 to 2.1185 | 0.0177 | 0.0530 |
| LIHC | 370 | 130 | 7 | 18.57 | fit_and_pool | 0.362574 | 0.106751 | **1.4370** | 1.1657 to 1.7715 | 0.0007 | 0.0041 |
| CHOL | 35 | 18 | 7 | 2.57 | **do_not_fit** | *not fitted* | | | | | |

### 2.2 CHOL, descriptive only — `output/chol_descriptive.csv`

| Model | n | Events | Params | EPV | Band | Fitted | β | SE | HR | 95% CI | p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 | 35 | 18 | 1 | 18.00 | fit_and_pool | TRUE | -0.160616 | 0.216782 | 0.8516 | 0.5568 to 1.3025 | 0.4587 |
| M2 | 35 | 18 | 5 | 3.60 | **do_not_fit** | FALSE | *not fitted* | | | | |
| M3 | 35 | 18 | 6 | 3.00 | **do_not_fit** | FALSE | *not fitted* | | | | |
| M4 | 35 | 18 | 7 | 2.57 | **do_not_fit** | FALSE | *not fitted* | | | | |

**CHOL fits M1 only.** M2–M4 are `do_not_fit` (EPV below the floor of 5), so CHOL
contributes no attenuation estimate and appears in no attenuation figure.
---

## 3. Attenuation per cohort — `output/attenuation_per_cohort.csv`

`attenuation_total` = β(M2) − β(M4); `purity` = β(M2) − β(M3); `stroma` = β(M3) −
β(M4). **Registered direction is positive** = adjustment reduces the association.
Intervals are percentile from 2,000 **paired** resamples (all four models refitted
per resample). `prop_attenuated` is reported only where β2 ≠ 0 and sign(β2) =
sign(β4).

| Cohort | n | Events | β M1 | β M2 | β M3 | β M4 | **atten. total** | 95% CI | purity | 95% CI | stroma | 95% CI | prop_att | boot ok/failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD | 455 | 102 | 0.005159 | 0.004596 | -0.019686 | -0.071333 | **0.075929** | -0.1519 to 0.2967 | 0.024282 | -0.1234 to 0.1718 | 0.051647 | -0.1143 to 0.2125 | — | 2000/0 |
| READ | 165 | 36 | -0.035099 | -0.045714 | 0.016465 | -0.022076 | **-0.023638** | -0.4363 to 0.3453 | -0.062179 | -0.3464 to 0.1256 | 0.038541 | -0.2542 to 0.3841 | 0.5171 (-16.6559 to 0.8462) | 2000/0 |
| STAD | 403 | 157 | 0.113122 | 0.075526 | 0.012236 | 0.008397 | **0.067129** | -0.0089 to 0.1653 | 0.063290 | -0.0073 to 0.1559 | 0.003839 | -0.0257 to 0.0408 | 0.8888 (-12.0112 to 0.9423) | 2000/0 |
| ESCA | 184 | 77 | 0.054184 | 0.146456 | 0.211618 | 0.207848 | **-0.061391** | -0.2327 to 0.0631 | -0.065162 | -0.2354 to 0.0660 | 0.003771 | -0.0718 to 0.0527 | -0.4192 (-5.3482 to 0.5844) | 2000/0 |
| PAAD | 178 | 93 | 0.276873 | 0.275399 | 0.457845 | 0.411072 | **-0.135673** | -0.5065 to 0.1065 | -0.182446 | -0.5335 to 0.0293 | 0.046773 | -0.0319 to 0.1398 | -0.4926 (-2.6777 to 0.4054) | 2000/0 |
| LIHC | 370 | 130 | 0.233098 | 0.233928 | 0.321546 | 0.362574 | **-0.128646** | -0.2775 to -0.0367 | -0.087618 | -0.2168 to 0.0024 | -0.041028 | -0.1172 to -0.0017 | -0.5499 (-2.8434 to -0.1257) | 2000/0 |

**Cohorts absent from this table:** CHOL — CHOL has no M2–M4 fit, so no attenuation is defined for it.

**`prop_attenuated` is undefined in 1 of 6 cohorts** (COAD) — the sign condition fails there.

---

## 4. Meta-analysis — `output/meta_analysis.csv`

`metafor::rma`, REML. **Primary interval is Hartung–Knapp** (`ci_lo`/`ci_hi`); the
Wald interval and the fixed-effect estimate are given for cross-check.

| Analysis | k | n | **Estimate** | **HK 95% CI** | Wald 95% CI | FE est | FE 95% CI | HR | HR HK CI | p | τ² | I² | Q | df | Q p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **M1** | 6 | 1755 | **0.121174** | **0.0002 to 0.2421** | 0.0296 to 0.2127 | 0.122397 | 0.0400 to 0.2048 | 1.1288 | 1.0002 to 1.2740 | 0.0497 | 0.002212 | 16.83 | 6.0779 | 5 | 0.2987 |
| **M2** | 6 | 1755 | **0.122869** | **0.0001 to 0.2457** | 0.0300 to 0.2158 | 0.123111 | 0.0390 to 0.2072 | 1.1307 | 1.0001 to 1.2785 | 0.0499 | 0.002178 | 16.09 | 6.0369 | 5 | 0.3026 |
| **M3** | 6 | 1755 | **0.157205** | **-0.0414 to 0.3558** | 0.0059 to 0.3085 | 0.143900 | 0.0455 to 0.2423 | 1.1702 | 0.9595 to 1.4273 | 0.0975 | 0.018588 | 53.99 | 10.9672 | 5 | 0.0520 |
| **M4** | 6 | 1755 | **0.152248** | **-0.0620 to 0.3665** | -0.0142 to 0.3187 | 0.145752 | 0.0430 to 0.2485 | 1.1644 | 0.9399 to 1.4426 | 0.1273 | 0.023679 | 57.45 | 11.9518 | 5 | 0.0355 |
| **attenuation_total** | 6 | 1755 | **-0.024766** | **-0.1253 to 0.0757** | -0.1163 to 0.0668 | -0.008059 | -0.0667 to 0.0506 | 0.9755 | 0.8823 to 1.0787 | 0.5543 | 0.005195 | 44.35 | 8.5395 | 5 | 0.1289 |
| **attenuation_purity** | 6 | 1755 | **-0.022711** | **-0.1080 to 0.0626** | -0.0963 to 0.0509 | -0.005177 | -0.0574 to 0.0470 | 0.9775 | 0.8976 to 1.0646 | 0.5241 | 0.003006 | 37.86 | 7.5615 | 5 | 0.1821 |
| **attenuation_stroma** | 6 | 1755 | **0.001124** | **-0.0252 to 0.0275** | -0.0229 to 0.0251 | 0.001124 | -0.0229 to 0.0251 | 1.0011 | 0.9751 to 1.0279 | 0.9170 | 0.000000 | 0.00 | 3.5036 | 5 | 0.6229 |
| **drop_lowEPV_M1** | 5 | 1590 | **0.133444** | **-0.0046 to 0.2715** | 0.0374 to 0.2295 | 0.133058 | 0.0480 to 0.2182 | 1.1428 | 0.9954 to 1.3119 | 0.0550 | 0.002437 | 20.25 | 5.1270 | 4 | 0.2745 |
| **drop_lowEPV_M2** | 5 | 1590 | **0.138181** | **0.0019 to 0.2745** | 0.0409 to 0.2354 | 0.136615 | 0.0492 to 0.2240 | 1.1482 | 1.0019 to 1.3159 | 0.0481 | 0.002232 | 18.05 | 4.7990 | 4 | 0.3086 |
| **drop_lowEPV_M3** | 5 | 1590 | **0.176803** | **-0.0676 to 0.4212** | 0.0064 to 0.3472 | 0.153720 | 0.0516 to 0.2558 | 1.1934 | 0.9346 to 1.5238 | 0.1150 | 0.022781 | 61.73 | 10.4707 | 4 | 0.0332 |
| **drop_lowEPV_M4** | 5 | 1590 | **0.172278** | **-0.0872 to 0.4317** | -0.0111 to 0.3556 | 0.155683 | 0.0499 to 0.2614 | 1.1880 | 0.9165 to 1.5399 | 0.1390 | 0.027296 | 64.07 | 11.3456 | 4 | 0.0229 |
| **attenuation_total_drop_lowEPV** | 5 | 1590 | **-0.025414** | **-0.1501 to 0.0992** | -0.1224 to 0.0716 | -0.007677 | -0.0670 to 0.0517 | 0.9749 | 0.8606 to 1.1043 | 0.6016 | 0.005822 | 51.60 | 8.5329 | 4 | 0.0739 |

### 4.1 Matched re-pool — `output/sensitivity_matched_repool.csv`

| Analysis | k | n | Estimate | HK CI | Wald CI | HR | τ² | I² | Q (df) | Q p |
|---|---|---|---|---|---|---|---|---|---|---|
| alt_endpoint_matched_M2 | 5 | 1590 | 0.082606 | -0.0066 to 0.1718 | -0.0024 to 0.1676 | 1.0861 | 0.000000 | 0.00 | 2.1932 (4) | 0.7003 |
| alt_endpoint_matched_M4 | 5 | 1590 | 0.101348 | -0.0337 to 0.2364 | -0.0002 to 0.2029 | 1.1067 | 0.000006 | 0.04 | 3.5290 (4) | 0.4735 |

---

## 5. Registered sensitivities — `output/sensitivity_meta.csv`

Every registered sensitivity, pooled. The primary M1/M2/M3/M4 rows in §4 are the
reference.

| Analysis (variant) | k | n | Estimate | HK 95% CI | Wald 95% CI | HR | HR CI | p | τ² | I² | Q (df) | Q p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **alt_endpoint_M1** | 6 | 1755 | 0.073979 | -0.0450 to 0.1930 | -0.0072 to 0.1551 | 1.0768 | 0.9560 to 1.2129 | 0.1710 | 0.000001 | 0.01 | 6.2546 (5) | 0.2822 |
| **alt_endpoint_M2** | 6 | 1755 | 0.081354 | 0.0091 to 0.1536 | -0.0015 to 0.1642 | 1.0848 | 1.0092 to 1.1660 | 0.0340 | 0.000000 | 0.00 | 2.2097 (5) | 0.8194 |
| **alt_endpoint_M3** | 5 | 1590 | 0.093237 | -0.0290 to 0.2155 | -0.0043 to 0.1907 | 1.0977 | 0.9714 to 1.2404 | 0.1016 | 0.000004 | 0.04 | 3.1336 (4) | 0.5357 |
| **alt_endpoint_M4** | 5 | 1590 | 0.101348 | -0.0337 to 0.2364 | -0.0002 to 0.2029 | 1.1067 | 0.9668 to 1.2667 | 0.1057 | 0.000006 | 0.04 | 3.5290 (4) | 0.4735 |
| **estimate_purity_M1** | 6 | 1755 | 0.121174 | 0.0002 to 0.2421 | 0.0296 to 0.2127 | 1.1288 | 1.0002 to 1.2740 | 0.0497 | 0.002212 | 16.83 | 6.0779 (5) | 0.2987 |
| **estimate_purity_M2** | 6 | 1755 | 0.122869 | 0.0001 to 0.2457 | 0.0300 to 0.2158 | 1.1307 | 1.0001 to 1.2785 | 0.0499 | 0.002178 | 16.09 | 6.0369 (5) | 0.3026 |
| **estimate_purity_M3** | 6 | 1755 | 0.169627 | -0.0602 to 0.3994 | -0.0091 to 0.3483 | 1.1849 | 0.9416 to 1.4910 | 0.1162 | 0.029994 | 63.07 | 13.7081 (5) | 0.0176 |
| **estimate_purity_M4** | 6 | 1755 | 0.151860 | -0.0692 to 0.3729 | -0.0200 to 0.3237 | 1.1640 | 0.9331 to 1.4520 | 0.1377 | 0.026228 | 59.63 | 12.5792 (5) | 0.0277 |
| **score_143_M1** | 6 | 1755 | 0.121748 | 0.0065 to 0.2370 | 0.0341 to 0.2094 | 1.1295 | 1.0065 to 1.2674 | 0.0420 | 0.001283 | 10.55 | 5.6114 (5) | 0.3459 |
| **score_143_M2** | 6 | 1755 | 0.122129 | 0.0042 to 0.2400 | 0.0326 to 0.2117 | 1.1299 | 1.0042 to 1.2713 | 0.0447 | 0.001384 | 10.93 | 5.6524 (5) | 0.3415 |
| **score_143_M3** | 6 | 1755 | 0.157238 | -0.0386 to 0.3531 | 0.0073 to 0.3072 | 1.1703 | 0.9622 to 1.4234 | 0.0939 | 0.017779 | 52.53 | 10.6241 (5) | 0.0594 |
| **score_143_M4** | 6 | 1755 | 0.153401 | -0.0595 to 0.3663 | -0.0130 to 0.3198 | 1.1658 | 0.9422 to 1.4424 | 0.1232 | 0.023303 | 56.61 | 11.7246 (5) | 0.0388 |
| **no_redacted_M1** | 6 | 1754 | 0.121394 | -0.0000 to 0.2428 | 0.0294 to 0.2134 | 1.1291 | 1.0000 to 1.2749 | 0.0500 | 0.002309 | 17.43 | 6.1164 (5) | 0.2951 |
| **no_redacted_M2** | 6 | 1754 | 0.122392 | 0.0002 to 0.2446 | 0.0301 to 0.2147 | 1.1302 | 1.0002 to 1.2771 | 0.0497 | 0.002019 | 15.08 | 5.9769 (5) | 0.3085 |
| **no_redacted_M3** | 6 | 1754 | 0.156609 | -0.0413 to 0.3545 | 0.0059 to 0.3073 | 1.1695 | 0.9595 to 1.4255 | 0.0976 | 0.018335 | 53.65 | 10.8845 (5) | 0.0537 |
| **no_redacted_M4** | 6 | 1754 | 0.151737 | -0.0618 to 0.3653 | -0.0142 to 0.3176 | 1.1639 | 0.9401 to 1.4409 | 0.1273 | 0.023399 | 57.15 | 11.8572 (5) | 0.0368 |

### 5.1 Leave-one-cohort-out — `output/meta_leave_one_cohort_out.csv`

| Analysis | k | Cohorts retained | n | Estimate | HK CI | HR | τ² | I² |
|---|---|---|---|---|---|---|---|---|
| LOCO_M1_minus_COAD | 5 | READ+STAD+ESCA+PAAD+LIHC | 1300 | 0.151411 | 0.0187 to 0.2842 | 1.1635 | 0.000003 | 0.03 |
| LOCO_M1_minus_READ | 5 | COAD+STAD+ESCA+PAAD+LIHC | 1590 | 0.133444 | -0.0046 to 0.2715 | 1.1428 | 0.002437 | 20.25 |
| LOCO_M1_minus_STAD | 5 | COAD+READ+ESCA+PAAD+LIHC | 1352 | 0.121086 | -0.0474 to 0.2895 | 1.1287 | 0.006616 | 35.01 |
| LOCO_M1_minus_ESCA | 5 | COAD+READ+STAD+PAAD+LIHC | 1571 | 0.131527 | -0.0225 to 0.2856 | 1.1406 | 0.004025 | 27.53 |
| LOCO_M1_minus_PAAD | 5 | COAD+READ+STAD+ESCA+LIHC | 1577 | 0.097889 | -0.0292 to 0.2250 | 1.1028 | 0.001351 | 11.38 |
| LOCO_M1_minus_LIHC | 5 | COAD+READ+STAD+ESCA+PAAD | 1385 | 0.090137 | -0.0434 to 0.2236 | 1.0943 | 0.000000 | 0.00 |
| LOCO_M2_minus_COAD | 5 | READ+STAD+ESCA+PAAD+LIHC | 1300 | 0.151111 | 0.0145 to 0.2877 | 1.1631 | 0.000217 | 1.79 |
| LOCO_M2_minus_READ | 5 | COAD+STAD+ESCA+PAAD+LIHC | 1590 | 0.138181 | 0.0019 to 0.2745 | 1.1482 | 0.002232 | 18.05 |
| LOCO_M2_minus_STAD | 5 | COAD+READ+ESCA+PAAD+LIHC | 1352 | 0.135932 | -0.0298 to 0.3017 | 1.1456 | 0.005170 | 28.92 |
| LOCO_M2_minus_ESCA | 5 | COAD+READ+STAD+PAAD+LIHC | 1571 | 0.118367 | -0.0433 to 0.2801 | 1.1257 | 0.005040 | 31.43 |
| LOCO_M2_minus_PAAD | 5 | COAD+READ+STAD+ESCA+LIHC | 1577 | 0.100081 | -0.0316 to 0.2317 | 1.1053 | 0.001312 | 10.74 |
| LOCO_M2_minus_LIHC | 5 | COAD+READ+STAD+ESCA+PAAD | 1385 | 0.090878 | -0.0461 to 0.2279 | 1.0951 | 0.000000 | 0.00 |
| LOCO_M3_minus_COAD | 5 | READ+STAD+ESCA+PAAD+LIHC | 1300 | 0.196439 | -0.0356 to 0.4284 | 1.2171 | 0.018064 | 52.28 |
| LOCO_M3_minus_READ | 5 | COAD+STAD+ESCA+PAAD+LIHC | 1590 | 0.176803 | -0.0676 to 0.4212 | 1.1934 | 0.022781 | 61.73 |
| LOCO_M3_minus_STAD | 5 | COAD+READ+ESCA+PAAD+LIHC | 1352 | 0.198122 | -0.0449 to 0.4411 | 1.2191 | 0.018611 | 49.81 |
| LOCO_M3_minus_ESCA | 5 | COAD+READ+STAD+PAAD+LIHC | 1571 | 0.148887 | -0.1123 to 0.4100 | 1.1605 | 0.026196 | 62.83 |
| LOCO_M3_minus_PAAD | 5 | COAD+READ+STAD+ESCA+LIHC | 1577 | 0.114248 | -0.0791 to 0.3076 | 1.1210 | 0.012275 | 45.88 |
| LOCO_M3_minus_LIHC | 5 | COAD+READ+STAD+ESCA+PAAD | 1385 | 0.112469 | -0.1206 to 0.3455 | 1.1190 | 0.013864 | 44.10 |
| LOCO_M4_minus_COAD | 5 | READ+STAD+ESCA+PAAD+LIHC | 1300 | 0.197343 | -0.0409 to 0.4356 | 1.2182 | 0.020292 | 53.63 |
| LOCO_M4_minus_READ | 5 | COAD+STAD+ESCA+PAAD+LIHC | 1590 | 0.172278 | -0.0872 to 0.4317 | 1.1880 | 0.027296 | 64.07 |
| LOCO_M4_minus_STAD | 5 | COAD+READ+ESCA+PAAD+LIHC | 1352 | 0.193535 | -0.0738 to 0.4608 | 1.2135 | 0.024773 | 53.53 |
| LOCO_M4_minus_ESCA | 5 | COAD+READ+STAD+PAAD+LIHC | 1571 | 0.141107 | -0.1410 to 0.4232 | 1.1515 | 0.033002 | 65.39 |
| LOCO_M4_minus_PAAD | 5 | COAD+READ+STAD+ESCA+LIHC | 1577 | 0.112661 | -0.1202 to 0.3456 | 1.1193 | 0.021364 | 56.90 |
| LOCO_M4_minus_LIHC | 5 | COAD+READ+STAD+ESCA+PAAD | 1385 | 0.092555 | -0.1393 to 0.3244 | 1.0970 | 0.012205 | 37.93 |
| attenuation_total_LOCO_minus_COAD | 5 | READ+STAD+ESCA+PAAD+LIHC | 1300 | -0.040349 | -0.1566 to 0.0759 | 0.9605 | 0.006121 | 50.51 |
| attenuation_total_LOCO_minus_READ | 5 | COAD+STAD+ESCA+PAAD+LIHC | 1590 | -0.025414 | -0.1501 to 0.0992 | 0.9749 | 0.005822 | 51.60 |
| attenuation_total_LOCO_minus_STAD | 5 | COAD+READ+ESCA+PAAD+LIHC | 1352 | -0.076268 | -0.1725 to 0.0200 | 0.9266 | 0.000000 | 0.00 |
| attenuation_total_LOCO_minus_ESCA | 5 | COAD+READ+STAD+PAAD+LIHC | 1571 | -0.017309 | -0.1501 to 0.1155 | 0.9828 | 0.007161 | 49.95 |
| attenuation_total_LOCO_minus_PAAD | 5 | COAD+READ+STAD+ESCA+LIHC | 1577 | -0.015984 | -0.1335 to 0.1016 | 0.9841 | 0.005621 | 50.07 |
| attenuation_total_LOCO_minus_LIHC | 5 | COAD+READ+STAD+ESCA+PAAD | 1385 | 0.018037 | -0.0796 to 0.1157 | 1.0182 | 0.001664 | 16.63 |

### 5.2 Attenuation under each sensitivity variant — `output/sensitivity_attenuation.csv`

| Variant | COAD | READ | STAD | ESCA | PAAD | LIHC |
|---|---|---|---|---|---|---|
| **alt_endpoint** | 0.071873 | — | 0.026247 | -0.009882 | -0.074892 | -0.109532 |
| **estimate_purity** | 0.101303 | -0.041398 | 0.067129 | -0.061392 | -0.135673 | -0.138481 |
| **no_redacted** | 0.075929 | -0.023638 | 0.067129 | -0.061392 | -0.135673 | -0.128449 |
| **score_143** | 0.077145 | -0.027170 | 0.070547 | -0.068341 | -0.138488 | -0.134885 |### 5.3 Where each requested sensitivity actually lives

All seven are present, but in **three different files**, and one is a *column* not
a row — worth stating so nobody reports it as missing:

| Requested sensitivity | File | Rows / columns |
|---|---|---|
| alternative endpoint | `sensitivity_meta.csv` | `alt_endpoint_M1`…`M4` |
| ESTIMATE purity | `sensitivity_meta.csv` | `estimate_purity_M1`…`M4` |
| 143-gene score | `sensitivity_meta.csv` | `score_143_M1`…`M4` |
| redaction-excluded | `sensitivity_meta.csv` | `no_redacted_M1`…`M4` |
| drop-low-EPV | **`meta_analysis.csv`** | `drop_lowEPV_M1`…`M4`, `attenuation_total_drop_lowEPV` |
| leave-one-out | **`meta_leave_one_cohort_out.csv`** | 30 rows |
| **fixed effect** | **every meta row** | the `fe_est`, `fe_ci_lo`, `fe_ci_hi` **columns** — not a separate analysis row |

⚠ **`alt_endpoint_M3` and `alt_endpoint_M4` pool k = 5, not 6** (COAD+STAD+ESCA+PAAD+LIHC —
READ drops out). Every other sensitivity pools all six. Do not compare those two
rows to the primary without noting the different cohort set.

---

## 6. Proportional-hazards diagnostics and collinearity

### 6.1 PH tests — `output/ph_tests.csv`

`p_score` is the Therneau–Grambsch test on the score term; `p_global` is the global
test. A row is flagged `violated_*` where p < 0.05.

| Cohort | Model | p (score) | p (global) | Violated (score) | Violated (global) |
|---|---|---|---|---|---|
| COAD | M1 | 0.3441 | 0.3441 | FALSE | FALSE |
| COAD | M2 | 0.4146 | **0.0423** | FALSE | **TRUE** |
| COAD | M3 | 0.4368 | 0.0547 | FALSE | FALSE |
| COAD | M4 | 0.4191 | **0.0073** | FALSE | **TRUE** |
| READ | M1 | 0.8315 | 0.8315 | FALSE | FALSE |
| READ | M2 | 0.7544 | 0.5881 | FALSE | FALSE |
| READ | M3 | 0.7838 | 0.6640 | FALSE | FALSE |
| READ | M4 | 0.7913 | 0.7276 | FALSE | FALSE |
| STAD | M1 | 0.0779 | 0.0779 | FALSE | FALSE |
| STAD | M2 | 0.0642 | 0.3984 | FALSE | FALSE |
| STAD | M3 | 0.0586 | 0.2086 | FALSE | FALSE |
| STAD | M4 | **0.0454** | 0.2121 | **TRUE** | FALSE |
| ESCA | M1 | 0.7565 | 0.7565 | FALSE | FALSE |
| ESCA | M2 | 0.3778 | 0.0975 | FALSE | FALSE |
| ESCA | M3 | 0.3558 | 0.0904 | FALSE | FALSE |
| ESCA | M4 | 0.3557 | 0.1216 | FALSE | FALSE |
| PAAD | M1 | 0.4826 | 0.4826 | FALSE | FALSE |
| PAAD | M2 | 0.6719 | 0.5817 | FALSE | FALSE |
| PAAD | M3 | 0.6128 | **0.0498** | FALSE | **TRUE** |
| PAAD | M4 | 0.5840 | 0.0825 | FALSE | FALSE |
| LIHC | M1 | **0.0489** | **0.0489** | **TRUE** | **TRUE** |
| LIHC | M2 | 0.1218 | **0.0021** | FALSE | **TRUE** |
| LIHC | M3 | 0.1155 | **0.0003** | FALSE | **TRUE** |
| LIHC | M4 | 0.1112 | **0.0001** | FALSE | **TRUE** |
| CHOL | M1 | 0.7158 | 0.7158 | FALSE | FALSE |

**Violating model–cohort pairs: 8 of 25.**
- **COAD M2** — global (p_score 0.4146, p_global 0.0423)
- **COAD M4** — global (p_score 0.4191, p_global 0.0073)
- **STAD M4** — score (p_score 0.0454, p_global 0.2121)
- **PAAD M3** — global (p_score 0.6128, p_global 0.0498)
- **LIHC M1** — score and global (p_score 0.0489, p_global 0.0489)
- **LIHC M2** — global (p_score 0.1218, p_global 0.0021)
- **LIHC M3** — global (p_score 0.1155, p_global 0.0003)
- **LIHC M4** — global (p_score 0.1112, p_global 0.0001)

### 6.2 Score × log-time interaction — `output/ph_sensitivity_score_logtime.csv`

| Cohort | Model | β score | SE | β (tt) | SE (tt) | p (tt) | Failed | Reason |
|---|---|---|---|---|---|---|---|---|
| COAD | M2 | -0.035138 | 0.375609 | 0.006735 | 0.061470 | 0.9127 | FALSE | — |
| COAD | M4 | -0.106595 | 0.381928 | 0.006038 | 0.061033 | 0.9212 | FALSE | — |
| STAD | M4 | -0.367264 | 0.435820 | 0.066900 | 0.076183 | 0.3799 | FALSE | — |
| PAAD | M3 | -0.436200 | 0.808697 | 0.153404 | 0.138400 | 0.2677 | FALSE | — |
| LIHC | M1 | 0.991316 | 0.383520 | -0.133102 | 0.066044 | 0.0439 | FALSE | — |
| LIHC | M2 | 0.853351 | 0.386704 | -0.109241 | 0.066769 | 0.1018 | FALSE | — |
| LIHC | M3 | 0.954503 | 0.396971 | -0.111491 | 0.067761 | 0.0999 | FALSE | — |
| LIHC | M4 | 1.013667 | 0.401794 | -0.114712 | 0.068524 | 0.0941 | FALSE | — |

### 6.3 M4 variance inflation — `output/vif_m4.csv`

| Cohort | age (GVIF) | purity (GVIF) | score (GVIF) | sex (GVIF) | stage_group (GVIF) | stromal_score (GVIF) |
|---|---|---|---|---|---|---|
| COAD | 1.076 | 2.002 | 2.024 | 1.066 | 1.073 | 2.610 |
| READ | 1.045 | 2.072 | 2.102 | 1.178 | 1.089 | 2.777 |
| STAD | 1.036 | 5.214 | 1.174 | 1.018 | 1.069 | 5.089 |
| ESCA | 1.068 | 7.957 | 1.514 | 1.067 | 1.158 | 8.026 |
| PAAD | 1.009 | 7.239 | 1.798 | 1.009 | 1.014 | 7.228 |
| LIHC | 1.159 | 1.774 | 1.392 | 1.096 | 1.096 | 1.681 |

**Maximum GVIF across all cohorts and terms: 8.026** (ESCA, `stromal_score`). Conventional concern thresholds are 5 or 10.

---

## 7. Null benchmark

### 7.1 p-values, all three configurations — `output/null_pvalues.csv`

The primary p is the **add-one two-sided empirical p**, denominator = `N_requested`
(so failed draws cannot inflate significance). The one-sided value is the
directional companion.

| Config | Role | Statistic | N req. | n pooled | Observed | **p two-sided (PRIMARY)** | hits | p one-sided | tail | granularity |
|---|---|---|---|---|---|---|---|---|---|---|
| `registered_152|registered` | **primary_registered** | p_crude_M1 | 10000 | 10000 | 0.121174 | **0.395260** | 3952 | 0.395260 | upper | 0.000100 |
| `registered_152|registered` | **primary_registered** | p_adjusted_M2 | 10000 | 10000 | 0.122869 | **0.473053** | 4730 | 0.473053 | upper | 0.000100 |
| `registered_152|registered` | **primary_registered** | p_atten | 10000 | 10000 | -0.024766 | **0.115688** | 1156 | 0.000100 | lower | 0.000100 |
| `tested_140|registered` | **size_sensitivity** | p_crude_M1 | 10000 | 10000 | 0.121174 | **0.265173** | 2651 | 0.265173 | upper | 0.000100 |
| `tested_140|registered` | **size_sensitivity** | p_adjusted_M2 | 10000 | 10000 | 0.122869 | **0.360664** | 3606 | 0.360664 | upper | 0.000100 |
| `tested_140|registered` | **size_sensitivity** | p_atten | 10000 | 10000 | -0.024766 | **0.126587** | 1265 | 0.000100 | lower | 0.000100 |
| `tested_140|EXPLORATORY_POSTHOC` | **exploratory_posthoc** | p_crude_M1 | 1000 | 1000 | 0.121174 | **0.721279** | 721 | 0.721279 | upper | 0.000999 |
| `tested_140|EXPLORATORY_POSTHOC` | **exploratory_posthoc** | p_adjusted_M2 | 1000 | 1000 | 0.122869 | **0.855145** | 855 | 0.855145 | upper | 0.000999 |
| `tested_140|EXPLORATORY_POSTHOC` | **exploratory_posthoc** | p_atten | 1000 | 1000 | -0.024766 | **0.022977** | 22 | 0.022977 | lower | 0.000999 |

### 7.2 Null distribution five-number summary — `output/null_distributions.csv`

| Config | Role | Statistic | N | Min | Q25 | Median | Mean | Q75 | Max | Observed | **Percentile** |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `registered_152|registered` | primary_registered | pooled_M1_logHR | 10000 | 0.064988 | 0.107915 | **0.117486** | 0.117481 | 0.126996 | 0.168142 | 0.121174 | **60.48** |
| `registered_152|registered` | primary_registered | pooled_M2_logHR | 10000 | 0.069253 | 0.112147 | **0.121734** | 0.121938 | 0.131689 | 0.178049 | 0.122869 | **52.70** |
| `registered_152|registered` | primary_registered | pooled_attenuation_total | 10000 | -0.008000 | 0.011152 | **0.015318** | 0.015244 | 0.019522 | 0.036920 | -0.024766 | **0.00** |
| `tested_140|registered` | size_sensitivity | pooled_M1_logHR | 10000 | 0.060477 | 0.102820 | **0.112376** | 0.112343 | 0.121797 | 0.163874 | 0.121174 | **73.49** |
| `tested_140|registered` | size_sensitivity | pooled_M2_logHR | 10000 | 0.060866 | 0.108127 | **0.117870** | 0.117945 | 0.127584 | 0.168054 | 0.122869 | **63.94** |
| `tested_140|registered` | size_sensitivity | pooled_attenuation_total | 10000 | -0.009747 | 0.011086 | **0.015430** | 0.015303 | 0.019747 | 0.037115 | -0.024766 | **0.00** |
| `tested_140|EXPLORATORY_POSTHOC` | exploratory_posthoc | pooled_M1_logHR | 1000 | 0.095350 | 0.119968 | **0.127263** | 0.127563 | 0.135024 | 0.158609 | 0.121174 | **27.90** |
| `tested_140|EXPLORATORY_POSTHOC` | exploratory_posthoc | pooled_M2_logHR | 1000 | 0.105109 | 0.126852 | **0.134670** | 0.134515 | 0.142203 | 0.165811 | 0.122869 | **14.50** |
| `tested_140|EXPLORATORY_POSTHOC` | exploratory_posthoc | pooled_attenuation_total | 1000 | -0.036531 | -0.014119 | **-0.009283** | -0.009601 | -0.005322 | 0.009675 | -0.024766 | **1.20** |

### 7.3 Proportion of null sets whose pooled CI excludes 1 — `output/null_ci_exclusion.csv`

| Config | Role | n | Proportion with CI excluding 1 |
|---|---|---|---|
| `registered_152|registered` | primary_registered | 10000 | **0.2808** |
| `tested_140|registered` | size_sensitivity | 10000 | **0.1745** |
| `tested_140|EXPLORATORY_POSTHOC` | exploratory_posthoc | 1000 | **0.9720** |

### 7.4 Per-cohort percentiles — `output/null_per_cohort_percentile.csv`

Primary registered configuration only (other configurations are in the file).

| Cohort | Model | Observed β | Null median | Null 2.5% | Null 97.5% | **Percentile** | N |
|---|---|---|---|---|---|---|---|
| COAD | M1 | 0.005159 | 0.033496 | -0.028680 | 0.097874 | **18.32** | 10000 |
| COAD | M2 | 0.004596 | 0.022290 | -0.034160 | 0.080470 | **26.97** | 10000 |
| COAD | M4 | -0.071333 | -0.004853 | -0.070797 | 0.065763 | **2.40** | 10000 |
| READ | M1 | -0.035099 | 0.039663 | -0.048939 | 0.138471 | **4.92** | 10000 |
| READ | M2 | -0.045714 | 0.022361 | -0.056969 | 0.110371 | **4.84** | 10000 |
| READ | M4 | -0.022076 | 0.030346 | -0.059750 | 0.130020 | **12.68** | 10000 |
| STAD | M1 | 0.113122 | 0.099449 | 0.043592 | 0.156138 | **67.86** | 10000 |
| STAD | M2 | 0.075526 | 0.065367 | 0.007421 | 0.125012 | **62.74** | 10000 |
| STAD | M4 | 0.008397 | 0.033873 | -0.006536 | 0.075537 | **11.07** | 10000 |
| ESCA | M1 | 0.054184 | -0.024980 | -0.079784 | 0.032426 | **99.75** | 10000 |
| ESCA | M2 | 0.146456 | 0.107623 | 0.040850 | 0.182645 | **85.20** | 10000 |
| ESCA | M4 | 0.207848 | 0.103045 | 0.027661 | 0.191161 | **98.89** | 10000 |
| PAAD | M1 | 0.276873 | 0.215003 | 0.129639 | 0.295215 | **93.17** | 10000 |
| PAAD | M2 | 0.275399 | 0.210664 | 0.134389 | 0.286793 | **95.13** | 10000 |
| PAAD | M4 | 0.411072 | 0.168354 | 0.073084 | 0.261945 | **100.00** | 10000 |
| LIHC | M1 | 0.233098 | 0.284737 | 0.212677 | 0.355266 | **7.66** | 10000 |
| LIHC | M2 | 0.233928 | 0.269006 | 0.195802 | 0.338969 | **16.41** | 10000 |
| LIHC | M4 | 0.362574 | 0.288065 | 0.210395 | 0.362116 | **97.56** | 10000 |
| CHOL | M1 | -0.160616 | -0.071697 | -0.146325 | 0.003184 | **0.99** | 10000 |---

## 8. CMS orthogonality (B.o) — colorectal only

### 8.1 Subtype distribution — `output/cms_distribution.csv`

| Cohort | CMS | n | % |
|---|---|---|---|
| COAD | CMS1 | 83 | 18.24 |
| READ | CMS1 | 16 | 9.70 |
| COAD | CMS2 | 157 | 34.51 |
| READ | CMS2 | 52 | 31.52 |
| COAD | CMS3 | 73 | 16.04 |
| READ | CMS3 | 33 | 20.00 |
| COAD | CMS4 | 125 | 27.47 |
| READ | CMS4 | 54 | 32.73 |
| COAD | unclassified | 17 | 3.74 |
| READ | unclassified | 10 | 6.06 |

### 8.2 Score-tertile × CMS cross-tabulation — `output/cms_tertile_crosstab.csv`


**COAD** — χ² = **122.9298**, df = **6**, p (asymptotic) = **&lt;0.0001**, p (simulated) = &lt;0.0001, min expected count = 24.167

| Tertile | CMS1 | CMS2 | CMS3 | CMS4 |
|---|---|---|---|---|
| T1 | 16 | 81 | 35 | 15 |
| T2 | 22 | 61 | 29 | 34 |
| T3 | 45 | 15 | 9 | 76 |

**READ** — χ² = **38.8027**, df = **6**, p (asymptotic) = **&lt;0.0001**, p (simulated) = &lt;0.0001, min expected count = 5.161

| Tertile | CMS1 | CMS2 | CMS3 | CMS4 |
|---|---|---|---|---|
| T1 | 3 | 27 | 16 | 4 |
| T2 | 6 | 19 | 9 | 19 |
| T3 | 7 | 6 | 8 | 31 |

### 8.3 CMS4 association and R² — `output/cms4_association.csv`

| Cohort | n | n CMS4 | r (point-biserial) | Mean score CMS4 | Mean score other | **R² of score on CMS** |
|---|---|---|---|---|---|---|
| COAD | 455 | 125 | **0.4054** | 0.6580 | -0.2492 | **0.2861** |
| READ | 165 | 54 | **0.4399** | 0.6288 | -0.3059 | **0.2381** |

### 8.4 Models with CMS added as a covariate — `output/cms_adjusted_models.csv`

| Cohort | Model | n | Events | β without CMS | SE | β with CMS | SE | Δβ | params | EPV | Band | Fitted | HR without | HR with | CMS levels |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| COAD | M2 | 455 | 102 | 0.004596 | 0.098124 | -0.111941 | 0.114982 | **-0.116537** | 9 | 11.33 | fit_and_pool | TRUE | 1.0046 | 0.8941 | CMS1/CMS2/CMS3/CMS4/unclassified |
| COAD | M4 | 455 | 102 | -0.071333 | 0.137779 | -0.137093 | 0.146600 | **-0.065760** | 11 | 9.27 | fit_pool_flag_LOO | TRUE | 0.9312 | 0.8719 | CMS1/CMS2/CMS3/CMS4/unclassified |
| READ | M2 | 165 | 36 | -0.045714 | 0.157689 | — | — | **—** | 9 | 4.00 | do_not_fit | FALSE | 0.9553 | — | CMS1/CMS2/CMS3/CMS4/unclassified |
| READ | M4 | 165 | 36 | -0.022076 | 0.221829 | — | — | **—** | 11 | 3.27 | do_not_fit | FALSE | 0.9782 | — | CMS1/CMS2/CMS3/CMS4/unclassified |

### 8.5 Stratified by CMS4 vs non-CMS4 — `output/cms_stratified_models.csv`

| Cohort | Stratum | Model | n | Events | EPV | Band | Fitted | β | SE | HR |
|---|---|---|---|---|---|---|---|---|---|---|
| COAD | CMS4 | M1 | 125 | 32 | 32.00 | fit_and_pool | TRUE | -0.231816 | 0.218280 | 0.7931 |
| COAD | CMS4 | M2 | 125 | 32 | 6.40 | fit_pool_flag_LOO | TRUE | -0.054235 | 0.230056 | 0.9472 |
| COAD | CMS4 | M4 | 125 | 32 | 4.57 | **do_not_fit** | FALSE | *not fitted* | | |
| COAD | non_CMS4 | M1 | 330 | 70 | 70.00 | fit_and_pool | TRUE | -0.046026 | 0.123030 | 0.9550 |
| COAD | non_CMS4 | M2 | 330 | 70 | 14.00 | fit_and_pool | TRUE | -0.066902 | 0.124057 | 0.9353 |
| COAD | non_CMS4 | M4 | 330 | 70 | 10.00 | fit_and_pool | TRUE | -0.067355 | 0.157762 | 0.9349 |
| READ | CMS4 | M1 | 54 | 14 | 14.00 | fit_and_pool | TRUE | -0.457968 | 0.351319 | 0.6326 |
| READ | CMS4 | M2 | 54 | 14 | 2.80 | **do_not_fit** | FALSE | *not fitted* | | |
| READ | CMS4 | M4 | 54 | 14 | 2.00 | **do_not_fit** | FALSE | *not fitted* | | |
| READ | non_CMS4 | M1 | 111 | 22 | 22.00 | fit_and_pool | TRUE | 0.034696 | 0.240933 | 1.0353 |
| READ | non_CMS4 | M2 | 111 | 22 | 4.40 | **do_not_fit** | FALSE | *not fitted* | | |
| READ | non_CMS4 | M4 | 111 | 22 | 3.14 | **do_not_fit** | FALSE | *not fitted* | | |

### 8.6 Attenuation with and without CMS — `output/cms_attenuation.csv`

| Cohort | Attenuation without CMS | Attenuation with CMS |
|---|---|---|
| COAD | 0.075929 | 0.025152 |
| READ | -0.023638 | — |

---

## 9. GSE39582 external validation — **never pooled with discovery**

### 9.1 Models — `output/validation_gse39582_models.csv`

| Model | n | Events | Params | EPV | Band | Fitted | β | SE | HR | 95% CI | p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 | 561 | 188 | 1 | 188.00 | fit_and_pool | TRUE | -0.039900 | 0.075759 | **0.9609** | 0.8283 to 1.1147 | 0.5984 |
| M2 | 561 | 188 | 5 | 37.60 | fit_and_pool | TRUE | -0.008495 | 0.073889 | **0.9915** | 0.8579 to 1.1461 | 0.9085 |
| M3 | 561 | 188 | 6 | 31.33 | fit_and_pool | TRUE | -0.133234 | 0.106959 | **0.8753** | 0.7097 to 1.0794 | 0.2129 |
| M4 | 561 | 188 | 7 | 26.86 | fit_and_pool | TRUE | -0.107981 | 0.111987 | **0.8976** | 0.7207 to 1.1180 | 0.3349 |

### 9.2 Attenuation — `output/validation_gse39582_attenuation.csv`

- **attenuation_total** = 0.099486 (95% CI -0.0764 to 0.2797, SE 0.089696)
- **attenuation_purity** = 0.124739 (95% CI -0.0334 to 0.2914, SE 0.083595)
- **attenuation_stroma** = -0.025254 (95% CI -0.0995 to 0.0478, SE 0.037545)
- β M2 -0.008495, β M3 -0.133234, β M4 -0.107981
- `prop_attenuated` = -11.710900
- n tumours 562, complete-case 561, purity NA 0
- bootstrap ok 2000, failed 0, nuisance unstable 32
- **pooled with discovery: `FALSE`**

### 9.3 PH and VIF

| Model | p (score) | p (global) | Violated score | Violated global |
|---|---|---|---|---|
| M1 | 0.5097 | 0.5097 | FALSE | FALSE |
| M2 | 0.5366 | 0.3460 | FALSE | FALSE |
| M3 | 0.5534 | 0.4347 | FALSE | FALSE |
| M4 | 0.5481 | 0.2112 | FALSE | FALSE |

| Term | GVIF | Df | GVIF^(1/2Df) |
|---|---|---|---|
| `score` | 2.329 | 1 | 1.526 |
| `age` | 1.024 | 1 | 1.012 |
| `sex` | 1.035 | 1 | 1.017 |
| `stage_group` | 1.036 | 2 | 1.009 |
| `purity` | 10.485 | 1 | 3.238 |
| `stromal_score` | 7.695 | 1 | 2.774 |

### 9.4 Sensitivities — `output/validation_gse39582_sensitivity.csv`

| Variant | Model | Fitted | β | SE | HR |
|---|---|---|---|---|---|
| **score_143_sensitivity** | M1 | TRUE | -0.041507 | 0.075733 | 0.9593 |
| **score_143_sensitivity** | M2 | TRUE | -0.009100 | 0.073955 | 0.9909 |
| **score_143_sensitivity** | M3 | TRUE | -0.137216 | 0.107754 | 0.8718 |
| **score_143_sensitivity** | M4 | TRUE | -0.112005 | 0.112831 | 0.8940 |
| **drop_KRT17** | M1 | TRUE | -0.044654 | 0.075771 | 0.9563 |
| **drop_KRT17** | M2 | TRUE | -0.012517 | 0.073935 | 0.9876 |
| **drop_KRT17** | M3 | TRUE | -0.142046 | 0.107017 | 0.8676 |
| **drop_KRT17** | M4 | TRUE | -0.117190 | 0.112304 | 0.8894 |

### 9.5 CIT subtype check — `output/validation_gse39582_cit_models.csv`

- registered: `FALSE` | CIT levels: 6
- **R² of score on CIT = 0.2503**
- β M4 without CIT -0.107981 → with CIT -0.041127 (SE 0.115495)
- params with CIT 12, EPV 15.67, band `fit_and_pool`, fitted TRUE

Cross-tabulation: χ² = **125.3572**, df = **10**, p (asymptotic) = **&lt;0.0001**, p (simulated) = &lt;0.0001, min expected = 19.667

| Tertile | C1 | C2 | C3 | C4 | C5 | C6 |
|---|---|---|---|---|---|---|
| T1 | 54 | 10 | 36 | 7 | 64 | 16 |
| T2 | 43 | 22 | 30 | 15 | 50 | 27 |
| T3 | 19 | 68 | 8 | 37 | 38 | 17 |---

## 10. Colorectal pool — **EXPLORATORY, POST-HOC**

Requested after the primary result was known. `output/exploratory_colorectal_pooled.csv`
and `output/exploratory_colorectal_inputs.csv`.

### 10.1 Inputs

| Cohort | n | Events | β | SE | Source |
|---|---|---|---|---|---|
| COAD | 455 | 102 | 0.005159 | 0.094338 | discovery |
| READ | 165 | 36 | -0.035099 | 0.166886 | discovery |
| GSE39582 | 561 | 188 | -0.039900 | 0.075759 | external validation |
| **sum** | **1181** | **326** | | | |

### 10.2 Pooled

| Quantity | Value |
|---|---|
| label | EXPLORATORY_POSTHOC |
| **k** | **3** (TCGA-COAD+TCGA-READ+GSE39582) |
| n total | **1181** |
| events total | **326** |
| estimate (log-HR) | **-0.023666** |
| **HK 95% CI** | -0.0879 to 0.0406 (p 0.2540) |
| **Wald 95% CI** | -0.1328 to 0.0855 (p 0.6708) |
| fixed effect | -0.023666 (-0.1328 to 0.0855) |
| HR | **0.9766** |
| HR, HK CI | 0.9158 to 1.0415 |
| HR, Wald CI | 0.8756 to 1.0892 |
| τ² | **0.000000** |
| I² | **0.00** |
| Q (df) | 0.1440 (2), p 0.9305 |
| prediction interval | -0.0879 to 0.0406 |
| **compatibility upper bound, HR per SD (Wald)** | **1.0892** |
| compatibility upper bound, HR per SD (HK) | 1.0415 |

**Note in the file:** EXPLORATORY, POST-HOC. NOT the registered discovery meta-analysis: GSE39582 is an external validation cohort and Amendment 16 states validation cohorts are not meta-analysed with discovery. With k = 3 and Q p = 0.931, the Hartung-Knapp interval is NARROWER than the Wald interval -- the known small-k behaviour when heterogeneity is near zero. The WALD/FE interval is the conservative one and is the bound quoted. Registered discovery pooled M1 for comparison: 0.121174 (HR 1.1288, HKSJ 1.0002-1.2740) over SIX cohorts; this row is two of those six plus one external cohort and is NOT a re-estimate of it. The three inputs are not exchangeable: COAD and READ are RNA-seq scored by expression_log2tpm, GSE39582 is GPL570 microarray scored through a probe-collapse path with a KRT17 rescue; and READ's registered endpoint is PFI (Amendment 7) while COAD and GSE39582 are OS. Pooling across differing endpoints and platforms is the substantive limitation, not the arithmetic.

⚠ **τ² = 0.000000 at k = 3.** In this regime the Hartung–Knapp interval is
**narrower** than Wald — the known small-k behaviour — so the **Wald bound is the
conservative one** and is what the figure quotes.

---

## 11. FU-iCCA — `output/fuicca_correlations.csv`

Amendment 16 PRIMARY validation, executed under Amendment 18.

| Comparison | Tier | Method | n | r | 95% CI | p | n-note |
|---|---|---|---|---|---|---|---|
| score_140 vs STAT3:Y705 | **PRIMARY** | pearson | 114 | **0.3491** | 0.1765 to 0.5009 | 0.0001 | 114 of 120 finite Y705 values are used: 6 of those patients have no mRNA column |
| score_140 vs STAT3:Y705 | **PRIMARY** | spearman | 114 | **0.3718** | 0.1957 to 0.5248 | &lt;0.0001 | 114 of 120 finite Y705 values are used: 6 of those patients have no mRNA column |
| score_140 vs STAT3:S727 | **secondary** | pearson | 132 | **0.1468** | -0.0247 to 0.3099 | 0.0931 | 132 of 135 finite S727 values are used: 3 of those patients have no mRNA column |
| score_140 vs STAT3:S727 | **secondary** | spearman | 132 | **0.0958** | -0.0768 to 0.2627 | 0.2748 | 132 of 135 finite S727 values are used: 3 of those patients have no mRNA column |
| score_140 vs STAT3 protein (S1D) | **secondary** | pearson | 208 | **0.2796** | 0.1492 to 0.4004 | &lt;0.0001 | — |
| score_140 vs STAT3 protein (S1D) | **secondary** | spearman | 208 | **0.3364** | 0.2063 to 0.4548 | &lt;0.0001 | — |
| score_139 vs STAT3:Y705 (CXCL8 dropped, A19 sensitivity) | **sensitivity** | pearson | 114 | **0.3501** | 0.1776 to 0.5017 | 0.0001 | — |
| score_139 vs STAT3:Y705 (CXCL8 dropped, A19 sensitivity) | **sensitivity** | spearman | 114 | **0.3707** | 0.1944 to 0.5238 | &lt;0.0001 | — |
| score_143 vs STAT3:Y705 | **sensitivity** | pearson | 114 | **0.3445** | 0.1714 to 0.4969 | 0.0002 | — |
| score_143 vs STAT3:Y705 | **sensitivity** | spearman | 114 | **0.3628** | 0.1858 to 0.5169 | &lt;0.0001 | — |
| score_140 vs STAT3 mRNA (S1C) | **secondary** | pearson | 255 | **0.5409** | 0.4478 to 0.6224 | &lt;0.0001 | — |
| score_140 vs STAT3 mRNA (S1C) | **secondary** | spearman | 255 | **0.5365** | 0.4359 to 0.6239 | &lt;0.0001 | — |
| score_140 vs ESTIMATE stromal | **secondary** | pearson | 255 | **0.5270** | 0.4321 to 0.6103 | &lt;0.0001 | — |
| score_140 vs ESTIMATE stromal | **secondary** | spearman | 255 | **0.5274** | 0.4258 to 0.6159 | &lt;0.0001 | — |
| ESTIMATE stromal vs STAT3:Y705 | **secondary** | pearson | 114 | **0.0853** | -0.1001 to 0.2651 | 0.3666 | — |
| ESTIMATE stromal vs STAT3:Y705 | **secondary** | spearman | 114 | **0.1372** | -0.0488 to 0.3140 | 0.1454 | — |

### 11.1 Per-gene correlation with pY705 — `output/fuicca_per_gene_y705.csv`

**140 genes** tested (the 140-gene primary list). **22 significant after
Benjamini–Hochberg** at 0.05. Ten largest \|r\|:

| Gene | n | r | 95% CI | p | p (BH) |
|---|---|---|---|---|---|
| `OSMR` | 114 | 0.3826 | 0.2138 to 0.5293 | &lt;0.0001 | 0.0031 |
| `TYK2` | 114 | 0.3729 | 0.2029 to 0.5211 | &lt;0.0001 | 0.0031 |
| `ETV6` | 114 | 0.3456 | 0.1726 to 0.4979 | 0.0002 | 0.0065 |
| `BCL6` | 114 | 0.3432 | 0.1700 to 0.4958 | 0.0002 | 0.0065 |
| `DNMT1` | 114 | 0.3166 | 0.1409 to 0.4730 | 0.0006 | 0.0168 |
| `REG1A` | 114 | 0.3090 | 0.1326 to 0.4664 | 0.0008 | 0.0176 |
| `IKBKE` | 114 | 0.3064 | 0.1298 to 0.4641 | 0.0009 | 0.0176 |
| `GRB2` | 114 | 0.2994 | 0.1222 to 0.4581 | 0.0012 | 0.0176 |
| `JAK2` | 114 | 0.2983 | 0.1210 to 0.4571 | 0.0013 | 0.0176 |
| `STAT3` | 114 | 0.2961 | 0.1187 to 0.4552 | 0.0014 | 0.0176 |

### 11.2 Missingness comparison — `output/fuicca_missingness.csv`

| Group | n | Mean score | SD | Median | Mean diff | t p | Wilcoxon p | Welch 95% CI | Tests |
|---|---|---|---|---|---|---|---|---|---|
| **Y705 measured** | 114 | -0.0511 | 0.8561 | 0.0002 | 0.0390 | 0.7719 | 0.9788 | -0.2263 to 0.3043 | association of Y705 missingness with the SCORE (predictor-side selection). Does NOT test Amendment 18's range-truncation claim, which concerns the unobserved Y705 values in the 94 excluded patients and is unfalsifiable from these data. |
| **Y705 missing** | 94 | -0.0901 | 1.0463 | 0.0059 | — | — | — | — | — |

**Gene coverage** — `output/fuicca_gene_coverage.csv`

| List | Genes | Present | Missing | Missing names | Present after rename | Missing after |
|---|---|---|---|---|---|---|
| primary_140 | 140 | 139 | 1 | CXCL8 | 140 | 0 |
| sensitivity_143 | 143 | 142 | 1 | CXCL8 | 143 | 0 |

---

## 12. TCGA RPPA — **EXPLORATORY, POST-HOC**

`output/exploratory_rppa_control_correlations.csv` (per cohort) and
`output/exploratory_rppa_control_pooled.csv` (pooled).

### 12.1 Per cohort, all four comparisons

| Cohort | Comparison | Method | n | r | 95% CI | p |
|---|---|---|---|---|---|---|
| COAD | 1. total STAT3 RPPA vs STAT3 mRNA | pearson | 345 | **0.3091** | 0.2104 to 0.4016 | &lt;0.0001 |
| READ | 1. total STAT3 RPPA vs STAT3 mRNA | pearson | 121 | **0.2608** | 0.0864 to 0.4198 | 0.0039 |
| STAD | 1. total STAT3 RPPA vs STAT3 mRNA | pearson | 327 | **0.4149** | 0.3209 to 0.5009 | &lt;0.0001 |
| ESCA | 1. total STAT3 RPPA vs STAT3 mRNA | pearson | 125 | **0.4470** | 0.2946 to 0.5773 | &lt;0.0001 |
| PAAD | 1. total STAT3 RPPA vs STAT3 mRNA | pearson | 111 | **0.0742** | -0.1137 to 0.2571 | 0.4387 |
| LIHC | 1. total STAT3 RPPA vs STAT3 mRNA | pearson | 181 | **0.3676** | 0.2343 to 0.4874 | &lt;0.0001 |
| COAD | 2. score_140 vs STAT3 mRNA | pearson | 356 | **0.6987** | 0.6413 to 0.7483 | &lt;0.0001 |
| READ | 2. score_140 vs STAT3 mRNA | pearson | 127 | **0.7951** | 0.7207 to 0.8514 | &lt;0.0001 |
| STAD | 2. score_140 vs STAT3 mRNA | pearson | 327 | **0.7433** | 0.6906 to 0.7882 | &lt;0.0001 |
| ESCA | 2. score_140 vs STAT3 mRNA | pearson | 125 | **0.5650** | 0.4323 to 0.6738 | &lt;0.0001 |
| PAAD | 2. score_140 vs STAT3 mRNA | pearson | 113 | **0.8144** | 0.7412 to 0.8685 | &lt;0.0001 |
| LIHC | 2. score_140 vs STAT3 mRNA | pearson | 181 | **0.6405** | 0.5456 to 0.7191 | &lt;0.0001 |
| COAD | 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | pearson | 356 | **0.6897** | 0.6309 to 0.7405 | &lt;0.0001 |
| READ | 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | pearson | 127 | **0.7889** | 0.7126 to 0.8467 | &lt;0.0001 |
| STAD | 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | pearson | 327 | **0.7353** | 0.6812 to 0.7815 | &lt;0.0001 |
| ESCA | 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | pearson | 125 | **0.5504** | 0.4149 to 0.6621 | &lt;0.0001 |
| PAAD | 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | pearson | 113 | **0.8097** | 0.7349 to 0.8650 | &lt;0.0001 |
| LIHC | 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | pearson | 181 | **0.6323** | 0.5359 to 0.7125 | &lt;0.0001 |
| COAD | 3. total STAT3 RPPA vs STAT3_pY705 RPPA | pearson | 345 | **-0.0358** | -0.1408 to 0.0701 | 0.5077 |
| READ | 3. total STAT3 RPPA vs STAT3_pY705 RPPA | pearson | 121 | **-0.0841** | -0.2587 to 0.0958 | 0.3591 |
| STAD | 3. total STAT3 RPPA vs STAT3_pY705 RPPA | pearson | 327 | **0.0502** | -0.0586 to 0.1578 | 0.3656 |
| ESCA | 3. total STAT3 RPPA vs STAT3_pY705 RPPA | pearson | 125 | **0.5433** | 0.4064 to 0.6563 | &lt;0.0001 |
| PAAD | 3. total STAT3 RPPA vs STAT3_pY705 RPPA | pearson | 111 | **0.0371** | -0.1503 to 0.2219 | 0.6992 |
| LIHC | 3. total STAT3 RPPA vs STAT3_pY705 RPPA | pearson | 181 | **0.2722** | 0.1316 to 0.4021 | 0.0002 |
| COAD | 4. score_140 vs STAT3_pY705 RPPA (from 13) | pearson | 356 | **0.1998** | 0.0979 to 0.2976 | 0.0001 |
| READ | 4. score_140 vs STAT3_pY705 RPPA (from 13) | pearson | 127 | **-0.0092** | -0.1831 to 0.1653 | 0.9183 |
| STAD | 4. score_140 vs STAT3_pY705 RPPA (from 13) | pearson | 327 | **0.0532** | -0.0556 to 0.1607 | 0.3379 |
| ESCA | 4. score_140 vs STAT3_pY705 RPPA (from 13) | pearson | 125 | **0.2073** | 0.0329 to 0.3695 | 0.0204 |
| PAAD | 4. score_140 vs STAT3_pY705 RPPA (from 13) | pearson | 113 | **0.1802** | -0.0047 to 0.3531 | 0.0562 |
| LIHC | 4. score_140 vs STAT3_pY705 RPPA (from 13) | pearson | 181 | **-0.0567** | -0.2009 to 0.0899 | 0.4483 |

### 12.2 Pooled — `output/exploratory_rppa_control_pooled.csv`

| Comparison | k | n | **r pooled** | Wald 95% CI | HK 95% CI | p (Wald) | τ² | **I²** | Q (df) | Q p | Can establish? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1. total STAT3 RPPA vs STAT3 mRNA | 6 | 1210 | **0.3247** | 0.2273 to 0.4157 | 0.1891 to 0.4481 | &lt;0.0001 | 0.011523 | **68.54** | 14.3673 (5) | 0.0134 | Protein-mRNA agreement for ONE gene across two independent assays. Typical TCGA protein-mRNA r is modest (~0.4-0.5 for well-behaved antibodies), so a value in that range is consistent with a working platform and a value near zero would indict the total-STAT3 antibody in that cohort. NO threshold is prespecified; this is descriptive. |
| 2. score_140 vs STAT3 mRNA | 6 | 1229 | **0.7189** | 0.6414 to 0.7818 | 0.6107 to 0.8007 | &lt;0.0001 | 0.026724 | **83.69** | 23.8934 (5) | 0.0002 | Panel internal coherence. PARTLY A SELF-CORRELATION: STAT3 is 1 of the 140 scoring genes (0.71% of the mean). Comparator: FU-iCCA 0.5409, which has the same property. See 2b for the leave-STAT3-out value. |
| 2b. score_139 (STAT3 dropped) vs STAT3 mRNA | 6 | 1229 | **0.7107** | 0.6308 to 0.7756 | 0.5992 to 0.7951 | &lt;0.0001 | 0.027230 | **83.95** | 24.1572 (5) | 0.0002 | As 2, with STAT3 excluded from the score. The difference from 2 is the self-correlation artefact. |
| 3. total STAT3 RPPA vs STAT3_pY705 RPPA | 6 | 1210 | **0.1393** | -0.0641 to 0.3315 | -0.1297 to 0.3891 | 0.1788 | 0.059064 | **91.78** | 48.8368 (5) | &lt;0.0001 | WITHIN-PLATFORM: both antibodies are measured on the SAME lysate, the same array, and are normalised together. It is NOT independent evidence that either antibody is valid -- shared loading and normalisation inflate it. A high value is uninformative about assay validity; a LOW value would be surprising and would suggest one antibody is noise. |
| 4. score_140 vs STAT3_pY705 RPPA (from 13) | 6 | 1229 | **0.0964** | 0.0041 to 0.1872 | -0.0242 to 0.2143 | 0.0407 | 0.007734 | **59.76** | 12.4967 (5) | 0.0286 | Reproduces 13's committed score-vs-pY705 correlations; the quantity under test. |

### 12.3 Coverage — `output/exploratory_rppa_coverage.csv`

| Cohort | RPPA files | RPPA patients | Dup aliquots dropped | Analysis-set patients | Overlap | Analysed | Note |
|---|---|---|---|---|---|---|---|
| COAD | 360 | 360 | 0 | 455 | **356** | TRUE | — |
| READ | 131 | 131 | 0 | 165 | **127** | TRUE | — |
| STAD | 357 | 357 | 0 | 406 | **327** | TRUE | — |
| ESCA | 126 | 126 | 0 | 184 | **125** | TRUE | — |
| PAAD | 120 | 120 | 0 | 178 | **113** | TRUE | — |
| LIHC | 184 | 184 | 0 | 370 | **181** | TRUE | — |
| CHOL | 0 | 0 | 0 | 35 | **0** | FALSE | descriptive only: not meta-eligible (Amendment 8) AND the author's prior ESMO Asia work used TCGA RPPA STAT3_pY705 in this cohort -- prior knowledge, disclosed |

### 12.4 Pooled heterogeneity, Fisher-z scale — `output/exploratory_rppa_pooled_heterogeneity.csv`

statistic `pooled Fisher z, random effects (REML)` | k 6 | n 1229 | z 0.096727 |
**r back-transformed 0.0964** | Wald CI 0.0041 to 0.1872 |
HK CI -0.0242 to 0.2143 | τ² 0.007734 | I² 59.76 |
Q 12.4967 (5) p 0.0286 | PI -0.1583 to 0.3391

Note: EXPLORATORY, POST-HOC. tanh of the pooled z is the back-transformed pooled correlation, not an unbiased mean of the per-cohort r. Wald interval quoted, as in the colorectal pooling. Cohorts differ in tissue and in n, and the RPPA antibody is one phosphosite assay.

---

## Nothing was computed for this document

Every number above is a transcription of the file named in its section heading. Two
derived *lists* are labelled as such (§1.4's 24 gene names, whose count reproduces
the committed `k_all3` = 24; and the row/column selections in the wide tables).
`output/null_pooled_draws.csv` does not exist — the per-draw values live in
`output/null_replicates.rds` and are exported to `figures/refined/` by
`11b_figures_refined.R`; §7 reports the committed summaries instead.
