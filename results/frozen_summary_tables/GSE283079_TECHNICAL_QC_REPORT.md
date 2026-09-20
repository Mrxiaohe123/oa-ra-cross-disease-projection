# GSE283079 Stage 2B Technical QC Report

## Scope

Only gene-level matrix construction and sample-level technical QC were performed. No RA score, DEG, GSEA, clustering for biology, or subtype analysis was performed.

## Matrix

Samples: 41; genes: 78277; OA: 36; non-OA: 5
Transcript-to-gene unmatched transcript rows: 0 in the matrix QC summary.

## Sample decisions

Retained: 41; excluded: 0
All samples passed paired FASTQ/FastQC/Salmon job completion checks.

## SRR31542944
Mapping rate: 71.2087%; detected genes: 17941; median sample correlation: 0.8073
Decision: QC-FLAG-BUT-KEEP. The sample is retained because the low mapping rate alone is not sufficient for exclusion and no independent pairing or FastQC failure was documented.

## Decision

TECHNICAL QC CONDITIONAL PASS

The matrix is suitable for the next pre-specified analysis stage, with SRR31542944 carried forward as a technical QC flag. Patient identity is not independently confirmed from the public manifest and no patient-level clinical association is claimed.
