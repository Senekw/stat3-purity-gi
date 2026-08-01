# stat3-purity-gi

Preregistered secondary analysis asking whether a STAT3 transcriptional activity
score carries prognostic information in gastrointestinal cancers independent of
tumour purity and stromal content.

**Preregistration:** [OSF registration tcvgb](https://osf.io/tcvgb/)
**Registered at commit:** `468c7b4`

## Status

Panel locked. Analysis plan v1.4 is the post-registration Amendment 9
implementation. No compartment fraction or survival model has been computed.

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
