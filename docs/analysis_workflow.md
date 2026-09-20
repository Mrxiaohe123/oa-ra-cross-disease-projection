# Analysis workflow

1. Define and freeze the 150/150 RA programme from GSE89408.
2. Reconstruct and QC GSE283079; use 36 OA samples for primary validation and retain 5 non-OA samples descriptively.
3. Score frozen programmes with log2(TPM+1), gene-wise z-scores and mean available genes; RA projection is RA-up minus RA-down.
4. Evaluate within-cohort coupling, cross-cohort portability, meta-analysis, leave-one-cohort-out stability, score sensitivities and competing models.
5. Rebuild figures from frozen small tables where possible.

Path-only portability substitutions are documented in CODE_PROVENANCE.md.
