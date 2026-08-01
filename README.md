# stat3-purity-gi

Preregistered secondary analysis asking whether a STAT3 transcriptional activity
score carries prognostic information in gastrointestinal cancers independent of
tumour purity and stromal content.

**Preregistration:** [OSF registration tcvgb](https://osf.io/tcvgb/), registered
2026-08-01 (frozen snapshot). Parent project: [rka4f](https://osf.io/rka4f/),
DOI 10.17605/OSF.IO/RKA4F.
**Registered at commit:** `468c7b4` (recorded in `commit.txt` on the parent
project; the frozen registration snapshot predates that file and archives the
pre-debug versions of `01_download.R` and `02_panel.R`)

## Status

Panel locked. Analysis plan **v1.5** is the post-registration implementation,
carrying Amendments 9 and 10: the compartment atlas set is **three** atlases
(GSE125449 liver/biliary, GSE178341 colorectal, Peng pancreatic) after GSE155698
and GSE183904 were both found to deposit no cell-type annotation. No compartment
fraction or survival model has been computed.

## Key documents

- `panel_definition.md` — gene panel prespecification, locked, with dated amendments
- `analysis_plan.md` — the preregistered analysis plan
- `feasibility_assessment.md` — pilot work and prior-art review
- `NOTES_FOR_REVIEW.md` — observations surfaced during setup, not acted on

## Data

No patient-level data is included. All datasets are public: TCGA (GDC), GEO,
ICGC, and BioSino NODE. Retrieval is scripted in `01_download.R`.

## Licence

Code is MIT licensed. The underlying datasets carry their own terms of use.
