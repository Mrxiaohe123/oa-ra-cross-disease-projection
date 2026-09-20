# Stage 2C Independent Validation Report

Date: 2026-09-12 13:45:44

## Scope and locked inputs
The analysis uses the Stage2B technically frozen GSE283079 gene TPM matrix. The primary analysis contains 36 OA samples; 5 non-OA samples are descriptive only. No genes were reselected and no subtype, DEG, GSEA, clustering, WGCNA, machine learning, ROC, drug prediction, or new database analysis was performed.
The primary score is the previously used cohort-background score: log2(TPM+1), row-wise z score across all 41 samples, mean available genes, and RA projection = frozen_ra_up - frozen_ra_down. A rank score is retained only as sensitivity output.

## Technical scope
Stage2B retained 41/41 samples: 36 OA and 5 non-OA. Patient-level identity was not independently confirmed in the public run manifest, so no patient-level clinical claim is made.
Program coverage ranged from 6 to 200 genes; coverage is reported in STAGE2C_program_coverage.csv.

## Primary OA36 results
RA projection SD = 0.742; range = [-1.261, 1.635]. This supports patient-to-patient variation as a quantitative property, not a discrete subtype.
Correlation with general inflammation: rho = 0.328; correlation with 10-gene APC: rho = 0.091; mean absolute correlation across comparator programs = 0.18. Full bootstrap CIs are in GSE283079_program_correlations.csv.
The general-only model R2 = 0.07; adding the fixed APC score changes R2 by 0.01. These are association and variance-partitioning results, not evidence of a RA-specific mechanism.

## Cross-cohort interpretation
The cross-cohort files combine the new GSE283079 OA36 result with the previously generated OA-only cohort relationships. Random-effects pooling uses Fisher z with DerSimonian-Laird tau2; leave-one-cohort-out results are provided. The added cohort is not merged with prior cohorts and is not treated as an independent sample within those cohorts.

## Decision
**C_NOT_SUPPORTED**
This decision is limited to molecular independent validation. A positive RA-derived projection is not called RA-specific; interpretation remains conditional on general inflammation, APC/MHC-II, and tissue-state programs.

## Files
Frozen input manifest, OA36 patient-level scores, non-OA descriptive scores, coverage, correlations with bootstrap CIs, competing models, cross-cohort relationships, random-effects meta-analysis, leave-one-cohort-out analysis, and execution logs are stored in this Stage2C directory.
