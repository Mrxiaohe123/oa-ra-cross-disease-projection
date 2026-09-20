options(stringsAsFactors = FALSE)

root_dir <- normalizePath(getwd())
result_dir <- file.path(root_dir, "04_结果", "Stage2C独立验证")
log_dir <- file.path(root_dir, "05_日志", "Stage2C独立验证")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(root_dir, "04_结果", "Stage2B基因矩阵", "GSE283079_gene_tpm.tsv.gz")
tx2gene_file <- file.path(root_dir, "04_结果", "Stage2B基因矩阵", "GENCODE_v47_tx2gene.tsv.gz")
manifest_file <- file.path(root_dir, "04_结果", "Stage2B技术QC", "GSE283079_analysis_sample_manifest_frozen.csv")
decision_file <- file.path(root_dir, "04_结果", "Stage2B技术QC", "sample_qc_decision.csv")
frozen_manifest_file <- file.path(root_dir, "结果", "Stage1B增量验证", "frozen_program_manifest.csv")
prior_corr_file <- file.path(root_dir, "04_结果", "OA内部连续谱分析", "oa_ra_program_correlations.csv")

message("Reading Stage2B matrix")
tpm <- read.delim(gzfile(matrix_file), check.names = FALSE, row.names = 1)
tpm <- as.matrix(tpm)
storage.mode(tpm) <- "numeric"
sample_manifest <- read.csv(manifest_file, check.names = FALSE)
qc_decision <- read.csv(decision_file, check.names = FALSE)
stopifnot(all(colnames(tpm) %in% sample_manifest$SRR))
sample_manifest <- sample_manifest[match(colnames(tpm), sample_manifest$SRR), ]
sample_manifest$group <- ifelse(sample_manifest$group %in% c("non_OA", "non-OA"), "non-OA", "OA")
sample_manifest$analysis_keep <- sample_manifest$final_quant_valid & sample_manifest$final_decision %in% c("KEEP", "QC-FLAG-BUT-KEEP")
oa_idx <- sample_manifest$group == "OA" & sample_manifest$analysis_keep
non_oa_idx <- sample_manifest$group == "non-OA" & sample_manifest$analysis_keep
if (sum(oa_idx) != 36L || sum(non_oa_idx) != 5L) stop("Unexpected Stage2B analysis groups")

tx2gene <- read.delim(gzfile(tx2gene_file), check.names = FALSE)
tx2gene$gene_id_clean <- sub("\\..*$", "", tx2gene$gene_id)
tx2gene$gene_symbol <- trimws(tx2gene$gene_symbol)
tx2gene <- tx2gene[!duplicated(tx2gene$gene_id_clean), c("gene_id_clean", "gene_symbol")]
gene_id_clean <- sub("\\..*$", "", rownames(tpm))
gene_symbol <- tx2gene$gene_symbol[match(gene_id_clean, tx2gene$gene_id_clean)]
gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- gene_id_clean[is.na(gene_symbol) | gene_symbol == ""]

frozen_manifest <- read.csv(frozen_manifest_file, check.names = FALSE)
program_names <- frozen_manifest$program
program_genes <- setNames(lapply(frozen_manifest$genes, function(x) unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))), program_names)

# Reproduce the established cohort-background scoring rule: log2(TPM+1),
# row-wise z score across all 41 samples, then mean of available genes.
log_tpm <- log2(tpm + 1)
row_mean <- rowMeans(log_tpm, na.rm = TRUE)
row_sd <- apply(log_tpm, 1, sd, na.rm = TRUE)
row_sd[!is.finite(row_sd) | row_sd == 0] <- NA_real_
zmat <- (log_tpm - row_mean) / row_sd
rownames(zmat) <- gene_symbol

score_program <- function(genes) {
  idx <- which(rownames(zmat) %in% genes)
  if (!length(idx)) return(rep(NA_real_, ncol(zmat)))
  colMeans(zmat[idx, , drop = FALSE], na.rm = TRUE)
}
scores_all <- as.data.frame(matrix(nrow = ncol(zmat), ncol = length(program_names)))
names(scores_all) <- program_names
for (p in program_names) scores_all[[p]] <- score_program(program_genes[[p]])
scores_all$SRR <- colnames(zmat)
scores_all$GSM <- sample_manifest$GSM
scores_all$group <- sample_manifest$group
scores_all$batch <- sample_manifest$batch
scores_all$ra_projection <- scores_all$frozen_ra_up - scores_all$frozen_ra_down
scores_all$direction_control <- scores_all$frozen_ra_up + scores_all$frozen_ra_down
scores_all$rank_score <- rank(scores_all$frozen_ra_up, ties.method = "average") / nrow(scores_all) - rank(scores_all$frozen_ra_down, ties.method = "average") / nrow(scores_all)
scores_all <- scores_all[, c("GSM", "SRR", "group", "batch", program_names, "ra_projection", "direction_control", "rank_score")]
scores_oa <- scores_all[oa_idx, ]
scores_non_oa <- scores_all[non_oa_idx, ]
write.csv(scores_oa, file.path(result_dir, "GSE283079_OA36_frozen_program_scores.csv"), row.names = FALSE, quote = TRUE)
write.csv(scores_non_oa, file.path(result_dir, "GSE283079_non_OA5_descriptive_program_scores.csv"), row.names = FALSE, quote = TRUE)
rank_sensitivity <- data.frame(
  n = nrow(scores_oa),
  spearman_ra_vs_rank_score = suppressWarnings(cor(scores_oa$ra_projection, scores_oa$rank_score, method = "spearman")),
  spearman_ra_vs_direction_control = suppressWarnings(cor(scores_oa$ra_projection, scores_oa$direction_control, method = "spearman")),
  ra_projection_sd = sd(scores_oa$ra_projection),
  rank_score_sd = sd(scores_oa$rank_score),
  stringsAsFactors = FALSE
)
write.csv(rank_sensitivity, file.path(result_dir, "GSE283079_score_sensitivity.csv"), row.names = FALSE)

coverage <- data.frame(
  program = program_names,
  frozen_gene_count = vapply(program_genes, length, integer(1)),
  covered_gene_count = vapply(program_genes, function(g) sum(g %in% rownames(zmat)), integer(1)),
  coverage_fraction = vapply(program_genes, function(g) mean(g %in% rownames(zmat)), numeric(1)),
  source = frozen_manifest$source,
  use = frozen_manifest$use,
  stringsAsFactors = FALSE
)
write.csv(coverage, file.path(result_dir, "STAGE2C_program_coverage.csv"), row.names = FALSE)

manifest_out <- data.frame(
  input_type = c("expression_matrix", "stage2b_sample_manifest", "stage2b_qc_decision", "frozen_program_manifest", "scoring_rule", "ra_reference"),
  path = c(matrix_file, manifest_file, decision_file, frozen_manifest_file, "log2(TPM+1), cohort-wide row-wise z, mean available genes; RA = up - down", "Stage1 frozen GSE89408 RA reference; 150 up and 150 down"),
  status = c("read_ok", "read_ok", "read_ok", "read_ok", "locked", "locked"),
  notes = c("Stage2B gene TPM; 41 samples", "36 OA and 5 non-OA; patient IDs unresolved", "all 41 retained; SRR31542944 QC-FLAG-BUT-KEEP", "no genes reselected", "rank_score retained as prespecified sensitivity", "not interpreted as RA-specific mechanism"),
  stringsAsFactors = FALSE
)
write.csv(manifest_out, file.path(result_dir, "STAGE2C_FROZEN_INPUT_MANIFEST.csv"), row.names = FALSE)

bootstrap_spearman <- function(x, y, n_boot = 2000, seed = 20260912) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 5) return(c(n = n, rho = NA, p = NA, ci_low = NA, ci_high = NA))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  p <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  set.seed(seed + n)
  boots <- replicate(n_boot, {
    ii <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[ii], y[ii], method = "spearman"))
  })
  c(n = n, rho = rho, p = p, ci_low = unname(quantile(boots, .025, na.rm = TRUE)), ci_high = unname(quantile(boots, .975, na.rm = TRUE)))
}

comparison_programs <- setdiff(program_names, c("frozen_ra_up", "frozen_ra_down"))
cor_rows <- lapply(comparison_programs, function(p) {
  z <- bootstrap_spearman(scores_oa$ra_projection, scores_oa[[p]], seed = 20260912 + match(p, comparison_programs))
  data.frame(cohort = "GSE283079", program = p, n = z[["n"]], spearman_rho = z[["rho"]], p = z[["p"]], ci_low = z[["ci_low"]], ci_high = z[["ci_high"]])
})
cor_gse <- do.call(rbind, cor_rows)
write.csv(cor_gse, file.path(result_dir, "GSE283079_program_correlations.csv"), row.names = FALSE)

fit_model <- function(formula_text, data, model_name) {
  fit <- lm(as.formula(formula_text), data = data)
  sm <- summary(fit)
  cf <- coef(summary(fit))
  ci <- suppressWarnings(confint(fit))
  out <- data.frame(cohort = "GSE283079", model = model_name, n = nrow(data), r_squared = sm$r.squared, adjusted_r_squared = sm$adj.r.squared, residual_sd = sd(residuals(fit)), stringsAsFactors = FALSE)
  for (term in rownames(cf)) {
    nm <- gsub("[^A-Za-z0-9]+", "_", term)
    out[[paste0("beta_", nm)]] <- cf[term, "Estimate"]
    out[[paste0("ci_low_", nm)]] <- ci[term, 1]
    out[[paste0("ci_high_", nm)]] <- ci[term, 2]
    out[[paste0("p_", nm)]] <- cf[term, "Pr(>|t|)"]
  }
  out
}
model_rows <- list()
model_rows[[1]] <- fit_model("ra_projection ~ general_inflammation", scores_oa, "general_only")
model_rows[[2]] <- fit_model("ra_projection ~ general_inflammation + mhc_ii_apc_primary", scores_oa, "general_plus_apc")
model_rows[[3]] <- fit_model("ra_projection ~ general_inflammation + mhc_ii_apc_primary + myeloid_context + interferon_state + t_cell_context + b_cell_context + fibroblast_ecm + osteoclast_context", scores_oa, "general_apc_cell_states")
all_model_columns <- unique(unlist(lapply(model_rows, names)))
model_rows <- lapply(model_rows, function(x) {
  miss <- setdiff(all_model_columns, names(x))
  for (m in miss) x[[m]] <- NA_real_
  x[, all_model_columns, drop = FALSE]
})
models <- do.call(rbind, model_rows)
models$delta_r_squared_vs_general <- models$r_squared - models$r_squared[models$model == "general_only"]
write.csv(models, file.path(result_dir, "GSE283079_competing_models.csv"), row.names = FALSE)

prior <- read.csv(prior_corr_file, check.names = FALSE)
prior <- prior[prior$program %in% comparison_programs, c("accession", "program", "n", "spearman_rho", "p", "ci_low", "ci_high")]
names(prior)[names(prior) == "accession"] <- "cohort"
relationships <- rbind(prior, cor_gse)
relationships$source <- ifelse(relationships$cohort == "GSE283079", "new_independent_validation", "previous_frozen_OA_analysis")
write.csv(relationships, file.path(result_dir, "OA_crosscohort_program_relationships.csv"), row.names = FALSE)

meta_re <- function(d) {
  d <- d[is.finite(d$spearman_rho) & is.finite(d$n) & d$n > 3, ]
  if (nrow(d) < 2) return(data.frame(k = nrow(d), pooled_rho = NA, ci_low = NA, ci_high = NA, tau2 = NA, i2 = NA))
  z <- atanh(pmax(pmin(d$spearman_rho, .999999), -.999999)); v <- 1 / (d$n - 3)
  w <- 1 / v; z_fixed <- sum(w * z) / sum(w); q <- sum(w * (z - z_fixed)^2)
  c0 <- sum(w) - sum(w^2) / sum(w); tau2 <- max(0, (q - (nrow(d) - 1)) / c0)
  wr <- 1 / (v + tau2); z_re <- sum(wr * z) / sum(wr); se <- sqrt(1 / sum(wr))
  ci <- z_re + c(-1, 1) * 1.96 * se
  data.frame(k = nrow(d), pooled_rho = tanh(z_re), ci_low = tanh(ci[1]), ci_high = tanh(ci[2]), tau2 = tau2, i2 = ifelse(q > 0, max(0, (q - (nrow(d)-1))/q) * 100, 0))
}
meta_rows <- lapply(comparison_programs, function(p) {
  d <- relationships[relationships$program == p, ]
  m <- meta_re(d)
  cbind(data.frame(program = p, stringsAsFactors = FALSE), m)
})
meta_table <- do.call(rbind, meta_rows)
write.csv(meta_table, file.path(result_dir, "OA_crosscohort_meta_analysis.csv"), row.names = FALSE)

loo_rows <- list()
for (p in comparison_programs) {
  d <- relationships[relationships$program == p, ]
  cohorts <- unique(d$cohort)
  for (omit in c("none", cohorts)) {
    dd <- if (omit == "none") d else d[d$cohort != omit, ]
    m <- meta_re(dd)
    loo_rows[[paste(p, omit, sep = "__")]] <- cbind(data.frame(program = p, omitted_cohort = omit, stringsAsFactors = FALSE), m)
  }
}
write.csv(do.call(rbind, loo_rows), file.path(result_dir, "OA_crosscohort_leave_one_cohort_out.csv"), row.names = FALSE)

mean_abs_r <- mean(abs(cor_gse$spearman_rho), na.rm = TRUE)
general_r <- cor_gse$spearman_rho[cor_gse$program == "general_inflammation"]
apc_r <- cor_gse$spearman_rho[cor_gse$program == "mhc_ii_apc_primary"]
general_model <- models[models$model == "general_only", ]
apc_model <- models[models$model == "general_plus_apc", ]
delta_apc <- apc_model$delta_r_squared_vs_general
if (is.finite(general_r) && is.finite(delta_apc) && abs(general_r) >= 0.5 && delta_apc < 0.10) {
  decision <- "B_CONTEXT-DEPENDENT SUPPORT"
} else if (is.finite(general_r) && is.finite(apc_r) && abs(general_r) >= 0.4 && abs(apc_r) >= 0.4) {
  decision <- "B_CONTEXT-DEPENDENT SUPPORT"
} else {
  decision <- "C_NOT_SUPPORTED"
}

report_lines <- c(
  "# Stage 2C Independent Validation Report",
  "",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Scope and locked inputs",
  "The analysis uses the Stage2B technically frozen GSE283079 gene TPM matrix. The primary analysis contains 36 OA samples; 5 non-OA samples are descriptive only. No genes were reselected and no subtype, DEG, GSEA, clustering, WGCNA, machine learning, ROC, drug prediction, or new database analysis was performed.",
  "The primary score is the previously used cohort-background score: log2(TPM+1), row-wise z score across all 41 samples, mean available genes, and RA projection = frozen_ra_up - frozen_ra_down. A rank score is retained only as sensitivity output.",
  "",
  "## Technical scope",
  paste0("Stage2B retained ", sum(sample_manifest$analysis_keep), "/41 samples: ", sum(oa_idx), " OA and ", sum(non_oa_idx), " non-OA. Patient-level identity was not independently confirmed in the public run manifest, so no patient-level clinical claim is made."),
  paste0("Program coverage ranged from ", min(coverage$covered_gene_count), " to ", max(coverage$covered_gene_count), " genes; coverage is reported in STAGE2C_program_coverage.csv."),
  "",
  "## Primary OA36 results",
  paste0("RA projection SD = ", round(sd(scores_oa$ra_projection), 3), "; range = [", round(min(scores_oa$ra_projection), 3), ", ", round(max(scores_oa$ra_projection), 3), "]. This supports patient-to-patient variation as a quantitative property, not a discrete subtype."),
  paste0("Correlation with general inflammation: rho = ", round(general_r, 3), "; correlation with 10-gene APC: rho = ", round(apc_r, 3), "; mean absolute correlation across comparator programs = ", round(mean_abs_r, 3), ". Full bootstrap CIs are in GSE283079_program_correlations.csv."),
  paste0("The general-only model R2 = ", round(general_model$r_squared, 3), "; adding the fixed APC score changes R2 by ", round(delta_apc, 3), ". These are association and variance-partitioning results, not evidence of a RA-specific mechanism."),
  "",
  "## Cross-cohort interpretation",
  "The cross-cohort files combine the new GSE283079 OA36 result with the previously generated OA-only cohort relationships. Random-effects pooling uses Fisher z with DerSimonian-Laird tau2; leave-one-cohort-out results are provided. The added cohort is not merged with prior cohorts and is not treated as an independent sample within those cohorts.",
  "",
  "## Decision",
  paste0("**", decision, "**"),
  "This decision is limited to molecular independent validation. A positive RA-derived projection is not called RA-specific; interpretation remains conditional on general inflammation, APC/MHC-II, and tissue-state programs.",
  "",
  "## Files",
  "Frozen input manifest, OA36 patient-level scores, non-OA descriptive scores, coverage, correlations with bootstrap CIs, competing models, cross-cohort relationships, random-effects meta-analysis, leave-one-cohort-out analysis, and execution logs are stored in this Stage2C directory."
)
writeLines(report_lines, file.path(result_dir, "STAGE2C_INDEPENDENT_VALIDATION_REPORT.md"))

writeLines(c(paste("Started", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "Stage2C frozen input validation executed", paste("Decision", decision), capture.output(sessionInfo())), file.path(log_dir, "stage2c_session_info.txt"))
writeLines(c("Stage2C executed", paste("OA n", sum(oa_idx)), paste("non-OA n", sum(non_oa_idx)), paste("Decision", decision)), file.path(log_dir, "stage2c_execution_status.txt"))
message("Stage2C complete: ", decision)
