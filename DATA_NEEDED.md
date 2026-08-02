# Data needed for script 10 (external validation, Amendment 16)

Everything below was verified against the live repositories on 2026-08-02 by
listing directories and reading series-matrix headers. Sizes and file names are
**as reported by the server**, not recalled. Where something could not be verified
because a repository is unreachable from this sandbox or requires an account, that
is stated as such rather than guessed.

**No account was created and no registration was attempted.** Two cohorts require
credentials that a human must obtain.

Target layout: `data/validation/<cohort>/`. Nothing here is downloaded yet.

---

## 1. FU-iCCA — PRIMARY validation analysis

Dong et al. 2022, *Cancer Cell* 40(1):70–87, `10.1016/j.ccell.2021.12.006`,
"Proteogenomic characterization identifies clinically relevant subgroups of
intrahepatic cholangiocarcinoma". Abstract confirms 262 patients, paired tumour
and adjacent liver.

| | |
|---|---|
| Repository | NODE (National Omics Data Encyclopedia), `https://www.biosino.org` |
| Accession | **OEP001105** |
| Project URL | `https://www.biosino.org/node/project/detail/OEP001105` |
| Needed | RNA-seq expression matrix (n up to 255); phosphoproteomic matrix including the **STAT3 pY705** site (n 214); sample-level clinical/linking table |
| **Registration** | **Likely required — a NODE account, and possibly author authorisation.** |
| **File names / sizes** | **NOT VERIFIED — see below.** |

**Why the file list is not here.** The project page is a JavaScript application
that returns a 906-byte shell with no file names in the HTML; there is no publicly
documented API (`/node/api/project/OEP001105`, `/api/project/detail/…`,
`/api/sample/list/…` all return `{"code":500,"msg":"404 NOT_FOUND"}`), and the
frontend bundle exposes only browser routes (`/project/detail/:id`,
`/sample/detail/:id`, `/download/node/data/:id`,
`/download/node/data/public/:id`). Listing the files therefore requires a browser
session. **Do not let anyone fill this table in from memory** — the file names
must be read off the project page.

The bundle also contains user-facing text describing **restricted data requiring
the submitter's authorisation** (`您必须向作者申请授权才能使用受限数据` — "you must
apply to the author for authorisation to use restricted data"), and separate
`/download/node/data/:id` versus `/download/node/data/public/:id` routes. Whether
OEP001105's files sit on the public or restricted path is the **first thing to
check**, because if they are restricted the primary validation analysis depends on
a request to the authors and has a lead time measured in weeks.

**Human steps:**
1. Open the project URL in a browser, register a NODE account if prompted.
2. Record, for each file: exact name, size, and whether it is on the public or the
   restricted download route.
3. If restricted, submit the authorisation request to the submitting author. The
   paper is closed-access (Unpaywall, Semantic Scholar and PMC all return no OA
   copy), so the data-availability statement must be read from the publisher PDF —
   it may name a contact or an alternative mirror.
4. Confirm the phosphoproteomics table actually reports **STAT3 pY705** as a
   distinct site. Amendment 16's primary analysis is defined on that phosphosite;
   if the deposit reports only protein-level STAT3, or a different site, the
   primary analysis as registered cannot be run and that must be raised **before**
   any substitute is considered.

---

## 2. GSE39582 — colorectal, SECONDARY

Marisa et al., "Gene expression classification of colon cancer into molecular
subtypes". **Verified on the GEO FTP server, 2026-08-02.**

| File | Size | URL |
|---|---|---|
| `GSE39582_series_matrix.txt.gz` | **164 MB** | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz` |
| `GSE39582_RAW.tar` (CEL files, only if re-normalising) | **4.4 GB** | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/suppl/GSE39582_RAW.tar` |
| `filelist.txt` | 39 KB | same `suppl/` directory |

- **Platform GPL570** (Affymetrix HG-U133 Plus 2.0), **585 samples**.
- **No registration.** Public FTP.
- **The series matrix is sufficient.** It carries both expression and clinical
  data; the 4.4 GB RAW tar is needed only to re-normalise from CEL files.

**Covariates present in the series matrix — verified by reading its 33
characteristic fields:** `Sex`, `age.at.diagnosis (year)`, `tnm.stage`, `tnm.t`,
`tnm.n`, `tnm.m`, `os.event`, `os.delay (months)`, `rfs.event`, `rfs.delay`,
plus `tumor.location`, `chemotherapy.adjuvant`, `mmr.status`, `cimp.status`,
`cin.status`, `tp53/kras/braf.mutation`, and `cit.molecularsubtype`.

**Consequence:** every covariate the B.j sequence needs (age, sex, stage, OS) is
present, so **M1–M4 can all be fitted** here, with purity ESTIMATE-derived per
Amendment 16. `cit.molecularsubtype` additionally allows a CMS-style orthogonality
check without re-running a classifier.

---

## 3. GSE66229 — gastric (ACRG), SECONDARY

**Verified on the GEO FTP server, 2026-08-02, and this cohort has a complication.**

| File | Size | URL |
|---|---|---|
| `GSE66229_series_matrix.txt.gz` | **108 MB** | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66229/matrix/GSE66229_series_matrix.txt.gz` |
| `GSE62254_series_matrix.txt.gz` (tumour subseries) | **85 MB** | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62254/matrix/GSE62254_series_matrix.txt.gz` |
| `GSE66229_RAW.tar` (CEL files, optional) | **1.8 GB** | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66229/suppl/GSE66229_RAW.tar` |

- **Platform GPL570**, **400 samples** in the SuperSeries.
- GSE66229 is a **SuperSeries**: `GSE62254` (300 tumours) + `GSE66222` (100
  normals). For validation use **GSE62254**, the tumour subseries — the 100
  normals are not part of the estimand.
- **No registration.** Public FTP.

**BLOCKER — survival data is not in GEO.** Both GSE66229 and GSE62254 carry
exactly **two** characteristic fields: `tissue` and `patient`. There is no
survival time, no event indicator, and no age, sex or stage. The `suppl/`
directory holds only `GSE62254_RAW.tar` and `filelist.txt` — no clinical table.

The ACRG clinical data must come from **the supplementary tables of the source
publication** (Cristescu et al. 2015, *Nat Med* 21:449–456, "Molecular analysis of
gastric cancer identifies subtypes associated with distinct clinical outcomes").
A human must download that supplement and confirm it contains a GSM-linkable
patient identifier.

**Consequence, to be decided before 10 is written:** if the supplement provides
survival but not age/sex/stage, then per Amendment 16 — "where a covariate is
unavailable, the model requiring it is not fitted and is reported as not fitted,
never approximated" — only **M1** is fittable in this cohort, and M2–M4 are
reported as not fitted. That is a legitimate outcome under the amendment, not a
reason to substitute proxies.

---

## 4. ICGC pancreatic (PACA-AU, PACA-CA) — SECONDARY

| | |
|---|---|
| Projects | **PACA-AU** (Australia), **PACA-CA** (Canada) |
| Needed | normalised expression matrix (`exp_seq` / `exp_array`), `donor` (survival), `specimen` (to restrict to primary tumours) |
| **Registration** | **Required for controlled-tier data — DACO / ICGC ARGO application.** Open-tier summary files do not include per-donor expression. |
| **File names / sizes** | **NOT VERIFIED — repository unreachable from this sandbox.** |

`dcc.icgc.org` and `docs.icgc-argo.org` are both blocked by this environment's
network allowlist, so no file listing could be obtained. Note also that the legacy
ICGC DCC portal has been retired and its holdings migrated; a human should confirm
the **current** host before downloading, rather than assume the legacy DCC paths
still resolve.

**Human steps:**
1. Determine the current home of PACA-AU / PACA-CA (ICGC ARGO platform, or the
   ICGC 25K legacy archive).
2. Establish whether the expression matrices are open-tier or controlled. If
   controlled, a DACO application is required and has a lead time of weeks.
3. Record exact file names and sizes once visible.
4. Confirm the expression unit. ICGC RNA-seq is typically distributed as
   normalised read counts rather than TPM; **B.h specifies log2(TPM+1)**, so any
   conversion is a deviation that needs a dated amendment before 10 is written.

---

## Feasibility already established

**The gene lists are fully measurable on the microarray cohorts.** Checked against
the GPL570 platform table (24,442 distinct gene symbols) fetched from GEO:

| List | Present on GPL570 | Missing |
|---|---|---|
| `final_gene_list_140.csv` (primary) | **140 / 140** | none |
| `final_gene_list_143.csv` (sensitivity) | **143 / 143** | none |

So GSE39582 and GSE62254 can both carry the full score. Probe-to-symbol collapse
will still be needed (GPL570 has multiple probes per gene), and per the fix
already applied in 07 that collapse must be **cohort-independent** — summing or
taking a fixed rule, never a within-cohort statistic such as maximum mean
expression, which would make the retained probe differ between cohorts.

## Summary of what blocks what

| Cohort | Role | Data reachable now | Blocking issue |
|---|---|---|---|
| FU-iCCA | **PRIMARY** | no | NODE account; possible author authorisation; file list must be read in a browser; pY705 presence unconfirmed |
| GSE39582 | secondary | **yes, fully** | none — 164 MB public download, all covariates present |
| GSE62254 | secondary | expression yes, clinical **no** | survival must come from the Cristescu 2015 supplement |
| ICGC PACA-AU/CA | secondary | no | likely controlled-access (DACO); host unconfirmed; expression unit may not be TPM |

**The primary validation analysis is the one with the hardest access path.** If
FU-iCCA proves restricted, that should be known before any secondary cohort is
downloaded, because it determines whether Amendment 16's primary analysis is
runnable at all.
