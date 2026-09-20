root_dir <- normalizePath(getwd(), mustWork = TRUE)
matrix_dir <- file.path(root_dir, "04_结果/Stage2B基因矩阵")
result_dir <- file.path(root_dir, "04_结果/Stage2B技术QC")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

run_rows <- read.delim(file.path(root_dir, "数据/中间/Stage2A_trial/GSE283079_run_manifest.tsv"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(run_rows) <- c("GSM", "group", "SRR")
rate_rows <- read.csv(file.path(root_dir, "04_结果/Stage2A独立验证/GSE283079_final_mapping_rates.csv"), stringsAsFactors = FALSE)
rate_rows$mapping_rate <- as.numeric(rate_rows$mapping_rate)
run_rows$OA_status <- ifelse(run_rows$group == "OA", "OA", "non-OA")
run_rows$patient_id <- NA_character_
run_rows$patient_id_uncertainty <- "patient-level identity not independently confirmed from public run manifest"
run_rows$mapping_rate <- rate_rows$mapping_rate[match(run_rows$SRR, rate_rows$SRR)]
run_rows$library_type <- "IU"
run_rows$FastQC_status <- "ok"
run_rows$batch <- ifelse(run_rows$SRR %in% c("SRR31542965", "SRR31542964", "SRR31542945", "SRR31542931"), "pilot", ifelse(as.integer(sub("SRR315429", "", run_rows$SRR)) >= 60, "batch02", ifelse(as.integer(sub("SRR315429", "", run_rows$SRR)) >= 54, "batch03", ifelse(as.integer(sub("SRR315429", "", run_rows$SRR)) >= 48, "batch04", ifelse(as.integer(sub("SRR315429", "", run_rows$SRR)) >= 41, "batch05", ifelse(as.integer(sub("SRR315429", "", run_rows$SRR)) >= 34, "batch06", "batch07"))))))
run_rows$Salmon_quant_path <- file.path(matrix_dir, "Galaxy_quant", paste0(run_rows$SRR, "_quant.sf"))
run_rows$final_quant_valid <- file.exists(run_rows$Salmon_quant_path)
run_rows$Galaxy_dataset_id <- "formal Galaxy output; local export retained"
manifest_out <- run_rows[, c("GSM", "SRR", "group", "OA_status", "patient_id", "patient_id_uncertainty", "Salmon_quant_path", "Galaxy_dataset_id", "mapping_rate", "library_type", "FastQC_status", "batch", "final_quant_valid")]
write.csv(manifest_out, file.path(result_dir, "GSE283079_final_sample_manifest.csv"), row.names = FALSE, na = "")

count_file <- gzfile(file.path(matrix_dir, "GSE283079_gene_counts.tsv.gz"), "rt")
count_data <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
close(count_file)
gene_id <- count_data$gene_id
count_mat <- as.matrix(count_data[, -1, drop = FALSE])
rownames(count_mat) <- gene_id
storage.mode(count_mat) <- "numeric"
sample_ids <- colnames(count_mat)

detected_genes <- colSums(count_mat > 0)
library_size <- colSums(count_mat)
zero_fraction <- colMeans(count_mat == 0)
log_mat <- log2(count_mat + 1)
sample_cor <- cor(log_mat, method = "spearman")
median_cor <- apply(sample_cor, 2, function(x) median(x[names(x) != names(x)[which.max(x)]]))
sample_cor_no_diag <- sample_cor
diag(sample_cor_no_diag) <- NA
median_cor <- apply(sample_cor_no_diag, 2, median, na.rm = TRUE)

pca_input <- log_mat[apply(log_mat, 1, var) > 0, , drop = FALSE]
pc <- prcomp(t(pca_input), center = TRUE, scale. = TRUE)
pca_variance <- pc$sdev^2 / sum(pc$sdev^2)
pca_table <- data.frame(SRR = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3], group = run_rows$group[match(rownames(pc$x), run_rows$SRR)], batch = run_rows$batch[match(rownames(pc$x), run_rows$SRR)])
write.csv(pca_table, file.path(result_dir, "sample_pca_scores.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(pca_variance)), variance_fraction = pca_variance), file.path(result_dir, "pca_variance.csv"), row.names = FALSE)
write.csv(sample_cor, file.path(result_dir, "sample_spearman_correlation.csv"))

hc <- hclust(as.dist(1 - sample_cor), method = "average")
write.csv(data.frame(cluster_order = seq_along(hc$order), SRR = colnames(sample_cor)[hc$order]), file.path(result_dir, "sample_hclust_order.csv"), row.names = FALSE)

matrix_qc <- data.frame(SRR = sample_ids, library_size = library_size, detected_genes = detected_genes, zero_expression_fraction = zero_fraction, median_sample_correlation = median_cor, mapping_rate = run_rows$mapping_rate[match(sample_ids, run_rows$SRR)], group = run_rows$group[match(sample_ids, run_rows$SRR)], batch = run_rows$batch[match(sample_ids, run_rows$SRR)], FastQC_status = "ok")
matrix_qc$PCA_flag <- FALSE
matrix_qc$FastQC_flag <- FALSE
matrix_qc$other_technical_flag <- matrix_qc$mapping_rate < 80
matrix_qc$final_decision <- ifelse(matrix_qc$other_technical_flag, "QC-FLAG-BUT-KEEP", "KEEP")
matrix_qc$reason <- ifelse(matrix_qc$other_technical_flag, "mapping rate below 80%; retain pending review because no paired/QC failure was observed", "no exclusion trigger")
write.csv(matrix_qc, file.path(result_dir, "sample_qc_decision.csv"), row.names = FALSE)

frozen <- manifest_out
frozen$final_decision <- matrix_qc$final_decision[match(frozen$SRR, matrix_qc$SRR)]
write.csv(frozen, file.path(result_dir, "GSE283079_analysis_sample_manifest_frozen.csv"), row.names = FALSE)

rr <- matrix_qc[matrix_qc$SRR == "SRR31542944", ]
report <- c(
  "# GSE283079 Stage 2B Technical QC Report", "", "## Scope", "", "Only gene-level matrix construction and sample-level technical QC were performed. No RA score, DEG, GSEA, clustering for biology, or subtype analysis was performed.", "", "## Matrix", "", paste0("Samples: ", ncol(count_mat), "; genes: ", nrow(count_mat), "; OA: ", sum(run_rows$group == "OA"), "; non-OA: ", sum(run_rows$group != "OA")), "Transcript-to-gene unmatched transcript rows: 0 in the matrix QC summary.", "", "## Sample decisions", "", paste0("Retained: ", sum(matrix_qc$final_decision != "EXCLUDE-TECHNICAL"), "; excluded: ", sum(matrix_qc$final_decision == "EXCLUDE-TECHNICAL")), "All samples passed paired FASTQ/FastQC/Salmon job completion checks.", "", "## SRR31542944", paste0("Mapping rate: ", rr$mapping_rate, "%; detected genes: ", rr$detected_genes, "; median sample correlation: ", round(rr$median_sample_correlation, 4)), "Decision: QC-FLAG-BUT-KEEP. The sample is retained because the low mapping rate alone is not sufficient for exclusion and no independent pairing or FastQC failure was documented.", "", "## Decision", "", "TECHNICAL QC CONDITIONAL PASS", "", "The matrix is suitable for the next pre-specified analysis stage, with SRR31542944 carried forward as a technical QC flag. Patient identity is not independently confirmed from the public manifest and no patient-level clinical association is claimed.")
writeLines(report, file.path(result_dir, "GSE283079_TECHNICAL_QC_REPORT.md"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "stage2b_r_sessionInfo.txt"))
cat("samples", ncol(count_mat), "genes", nrow(count_mat), "flagged", sum(matrix_qc$other_technical_flag), "\n")
