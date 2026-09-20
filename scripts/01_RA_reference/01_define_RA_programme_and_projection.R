#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(AnnotationDbi)
  library(hgu133a.db)
  library(hgu133plus2.db)
})

root_dir <- normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(root_dir, "data")
raw_dir <- file.path(root_dir, "数据", "原始")
meta_dir <- file.path(root_dir, "数据", "元数据")
result_dir <- file.path(root_dir, "结果", "阶段1")
log_dir <- file.path(root_dir, "日志", "阶段1")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  "Stage 1 scope: cross-disease transcriptomic projection",
  "No WGCNA, LASSO, ROC, docking, drug prediction, or subtype cutoff was used.",
  "RA reference: GSE89408 RA versus healthy, frozen before projection.",
  "General inflammation control: HALLMARK_INFLAMMATORY_RESPONSE gene set.",
  "Scores are continuous and cohort-standardized; effect estimates are descriptive."
), file.path(result_dir, "analysis_scope.txt"))

make_gene_matrix <- function(expr, platform) {
  if (platform == "GPL96") {
    annotation_package <- hgu133a.db
  } else if (platform == "GPL570") {
    annotation_package <- hgu133plus2.db
  } else {
    stop("Unsupported platform: ", platform)
  }
  key_type <- if ("PROBEID" %in% keytypes(annotation_package)) "PROBEID" else keytypes(annotation_package)[1]
  annotation_table <- AnnotationDbi::select(annotation_package,
                                             keys = rownames(expr),
                                             columns = "SYMBOL",
                                             keytype = key_type)
  annotation_table <- annotation_table[!is.na(annotation_table$SYMBOL) & annotation_table$SYMBOL != "", ]
  annotation_table <- annotation_table[!duplicated(annotation_table[[key_type]]), ]
  gene_symbol <- annotation_table$SYMBOL[match(rownames(expr), annotation_table[[key_type]])]
  keep <- !is.na(gene_symbol)
  expr <- expr[keep, , drop = FALSE]
  gene_symbol <- gene_symbol[keep]
  expr <- rowsum(expr, group = gene_symbol, reorder = FALSE) /
    as.vector(table(gene_symbol)[unique(gene_symbol)])
  list(expr = expr, mapped_features = sum(keep), total_features = length(keep), annotation = platform)
}

to_log_expression <- function(expr) {
  if (max(expr, na.rm = TRUE) > 50) log2(expr + 1) else expr
}

row_z <- function(expr) {
  centered <- expr - rowMeans(expr, na.rm = TRUE)
  spread <- apply(expr, 1, sd, na.rm = TRUE)
  spread[spread == 0 | is.na(spread)] <- 1
  centered / spread
}

score_set <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  colMeans(row_z(expr)[genes, , drop = FALSE], na.rm = TRUE)
}

read_rna_reference <- function() {
  meta_file <- file.path(root_dir, "结果", "阶段0.7", "GSE89408_patient_group_mapping.csv")
  if (!file.exists(meta_file)) meta_file <- file.path(root_dir, "结果", "Stage0.7", "GSE89408_patient_group_mapping.csv")
  meta <- fread(meta_file)
  raw <- fread(cmd = paste("gzip -cd", shQuote(file.path(raw_dir, "GSE89408_GEO_count_matrix_rename.txt.gz"))),
               data.table = FALSE, check.names = FALSE)
  gene_id <- raw[[1]]
  raw[[1]] <- NULL
  expr <- as.matrix(raw)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- sub(":.*$", "", gene_id)
  gene_ok <- grepl("^[A-Za-z][A-Za-z0-9.-]*$", rownames(expr))
  expr <- expr[gene_ok, , drop = FALSE]
  expr <- rowsum(expr, group = rownames(expr), reorder = FALSE) /
    as.vector(table(rownames(expr))[unique(rownames(expr))])
  normalize_label <- function(value) {
    value <- tolower(value)
    value <- gsub("healthy", "normal", value)
    value <- gsub("osteoarthritis", "oa", value)
    value <- gsub("rheumatoid arthritis", "ra", value)
    value <- gsub("arthralgia", "ag", value)
    value <- gsub("undifferentiated arthritis", "undiff", value)
    value <- gsub("tissue|synovial biopsy|synovial", "", value)
    gsub("[^a-z0-9]+", "_", value)
  }
  group <- meta$group[match(normalize_label(colnames(expr)), normalize_label(meta$title))]
  if (anyNA(group)) stop("GSE89408 RNA-seq group mapping failed")
  list(expr = expr, group = group, sample = colnames(expr), platform = "RNA-seq")
}

read_array_cohort <- function(accession) {
  object <- readRDS(file.path(data_dir, paste0(accession, "_parsed.rds")))
  expr <- to_log_expression(object$expr)
  mapped <- make_gene_matrix(expr, object$platform)
  list(expr = mapped$expr, group = object$pheno$group, sample = object$pheno$sample,
       title = object$pheno$title, platform = object$platform,
       mapped_features = mapped$mapped_features, total_features = mapped$total_features)
}

rna_reference <- read_rna_reference()
reference_keep <- rna_reference$group %in% c("RA", "healthy")
reference_expr <- rna_reference$expr[, reference_keep, drop = FALSE]
reference_group <- factor(rna_reference$group[reference_keep], levels = c("healthy", "RA"))
reference_dge <- DGEList(counts = reference_expr)
reference_dge <- calcNormFactors(reference_dge)
reference_voom <- voom(reference_dge, model.matrix(~ reference_group), plot = FALSE)
reference_fit <- eBayes(lmFit(reference_voom, model.matrix(~ reference_group)))
reference_table <- topTable(reference_fit, coef = 2, number = Inf, sort.by = "P")
reference_table$gene_symbol <- rownames(reference_table)
fwrite(reference_table, file.path(result_dir, "GSE89408_RA_reference_statistics.csv"))

select_reference <- function(table_data, direction, n_genes = 150) {
  if (direction == "up") {
    selected <- rownames(table_data)[table_data$adj.P.Val < 0.05 & table_data$logFC > 0]
    if (length(selected) < n_genes) selected <- rownames(table_data)[order(table_data$t, decreasing = TRUE)]
  } else {
    selected <- rownames(table_data)[table_data$adj.P.Val < 0.05 & table_data$logFC < 0]
    if (length(selected) < n_genes) selected <- rownames(table_data)[order(table_data$t)]
  }
  head(selected, n_genes)
}

ra_up <- select_reference(reference_table, "up")
ra_down <- select_reference(reference_table, "down")
generic_genes <- if (file.exists(file.path(meta_dir, "frozen_HALLMARK_INFLAMMATORY_RESPONSE_genes.txt"))) {
  readLines(file.path(meta_dir, "frozen_HALLMARK_INFLAMMATORY_RESPONSE_genes.txt"), warn = FALSE)
} else {
  stop("Frozen generic inflammation program is missing")
}
writeLines(ra_up, file.path(result_dir, "frozen_ra_reference_up_genes.txt"))
writeLines(ra_down, file.path(result_dir, "frozen_ra_reference_down_genes.txt"))
writeLines(sort(unique(generic_genes)), file.path(result_dir, "frozen_general_inflammation_genes.txt"))

program_list <- list(ra_up = ra_up, ra_down = ra_down)
program_overlap <- rbindlist(lapply(names(program_list), function(set_name) {
  set_data <- program_list[[set_name]]
  data.table(program = set_name, n_genes = length(set_data),
             overlap_general_inflammation = length(intersect(set_data, generic_genes)),
             overlap_fraction = length(intersect(set_data, generic_genes)) / length(set_data))
}), fill = TRUE)
fwrite(program_overlap, file.path(result_dir, "program_overlap_audit.csv"))

cell_sets <- list(
  immune_content = c("PTPRC", "CD3D", "CD3E", "MS4A1", "CD79A", "LST1", "TYROBP", "FCER1G"),
  myeloid_content = c("LST1", "TYROBP", "FCER1G", "CTSS", "FCGR3A", "AIF1"),
  t_cell_content = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R"),
  b_cell_content = c("MS4A1", "CD79A", "CD37", "CD74", "HLA-DRA"),
  fibroblast_content = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDPN"),
  endothelial_content = c("VWF", "KDR", "EMCN", "PECAM1", "ENG"),
  inflammatory_state = c("IL1B", "TNF", "CXCL8", "CCL2", "NFKBIA", "JUN", "FOS"),
  interferon_state = c("ISG15", "IFI6", "IFIT1", "IFIT3", "OAS1", "MX1")
)

cohort_names <- c("GSE55235", "GSE55457", "GSE55584", "GSE206848", "GSE82107")
cohorts <- setNames(lapply(cohort_names, read_array_cohort), cohort_names)
cohorts[["GSE89408"]] <- list(expr = to_log_expression(rna_reference$expr), group = rna_reference$group,
                               sample = rna_reference$sample, title = rna_reference$sample,
                               platform = "RNA-seq", mapped_features = nrow(rna_reference$expr),
                               total_features = nrow(rna_reference$expr))

platform_audit <- rbindlist(lapply(names(cohorts), function(name) {
  item <- cohorts[[name]]
  data.table(accession = name, platform = item$platform, total_features = item$total_features,
             mapped_features = item$mapped_features, mapped_fraction = item$mapped_features / item$total_features,
             samples = ncol(item$expr), group_summary = paste(names(table(item$group)), as.integer(table(item$group)), sep = "=", collapse = ";"))
}), fill = TRUE)
fwrite(platform_audit, file.path(result_dir, "platform_annotation_audit.csv"))

extract_patient_token <- function(title) {
  token <- sub(".*patient[ ]+", "", title, ignore.case = TRUE)
  token <- sub("[ ]*\\(.*$", "", token)
  token <- sub("[ ]+$", "", token)
  ifelse(grepl("patient", title, ignore.case = TRUE), token, title)
}
independence_audit <- rbindlist(lapply(names(cohorts), function(name) {
  item <- cohorts[[name]]
  data.table(accession = name, n_samples = length(item$sample), unique_sample_ids = uniqueN(item$sample),
             duplicate_sample_ids = length(item$sample) - uniqueN(item$sample),
             unique_patient_tokens = uniqueN(extract_patient_token(item$title)),
             duplicate_patient_tokens = length(item$sample) - uniqueN(extract_patient_token(item$title)))
}), fill = TRUE)
fwrite(independence_audit, file.path(result_dir, "sample_independence_audit.csv"))

all_tokens <- lapply(names(cohorts), function(name) extract_patient_token(cohorts[[name]]$title))
names(all_tokens) <- names(cohorts)
token_overlap <- rbindlist(lapply(seq_along(all_tokens), function(i) {
  if (i == length(all_tokens)) return(NULL)
  rbindlist(lapply((i + 1):length(all_tokens), function(j) {
    data.table(cohort_a = names(all_tokens)[i], cohort_b = names(all_tokens)[j],
               shared_patient_tokens = paste(intersect(all_tokens[[i]], all_tokens[[j]]), collapse = ";"),
               n_shared = length(intersect(all_tokens[[i]], all_tokens[[j]])))
  }))
}), fill = TRUE)
fwrite(token_overlap, file.path(result_dir, "cross_cohort_patient_token_overlap.csv"))

score_rows <- list()
gsea_rows <- list()
effect_rows <- list()
sensitivity_rows <- list()
specificity_rows <- list()

effect_size <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(c(n_x = length(x), n_y = length(y), difference = NA, d = NA))
  pooled <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
  c(n_x = length(x), n_y = length(y), difference = mean(x) - mean(y), d = (mean(x) - mean(y)) / pooled)
}

for (name in names(cohorts)) {
  item <- cohorts[[name]]
  expr <- item$expr
  score_table <- data.table(accession = name, sample = item$sample, group = item$group,
                            ra_convergence_score = score_set(expr, ra_up) - score_set(expr, ra_down),
                            general_inflammation_score = score_set(expr, generic_genes))
  for (set_name in names(cell_sets)) score_table[[set_name]] <- score_set(expr, cell_sets[[set_name]])
  score_rows[[name]] <- score_table

  for (contrast_name in c("RA_vs_OA", "OA_vs_Control")) {
    groups <- strsplit(contrast_name, "_vs_", fixed = TRUE)[[1]]
    if (!all(groups %in% unique(item$group))) next
    keep <- item$group %in% groups
    y <- factor(item$group[keep], levels = c(groups[2], groups[1]))
    fit <- eBayes(lmFit(expr[, keep, drop = FALSE], model.matrix(~ y)))
    stats <- fit$t[, 2]
    names(stats) <- rownames(expr)
    pathways <- list(ra_up = intersect(ra_up, names(stats)), ra_down = intersect(ra_down, names(stats)),
                     general_inflammation = intersect(generic_genes, names(stats)))
    pathways <- pathways[vapply(pathways, length, integer(1)) >= 10]
    if (length(pathways)) {
      gsea <- as.data.table(fgsea(pathways = pathways, stats = sort(stats, decreasing = TRUE),
                                  minSize = 10, maxSize = 1000, eps = 0))
      gsea[, `:=`(accession = name, contrast = contrast_name)]
      gsea_rows[[paste(name, contrast_name)]] <- gsea
    }
  }

  for (score_name in c("ra_convergence_score", "general_inflammation_score")) {
    for (group_name in c("RA_vs_OA", "OA_vs_Control")) {
      groups <- strsplit(group_name, "_vs_", fixed = TRUE)[[1]]
      if (!all(groups %in% unique(item$group))) next
      values <- score_table[[score_name]]
      effect <- effect_size(values[item$group == groups[1]], values[item$group == groups[2]])
      effect_rows[[paste(name, score_name, group_name)]] <- data.table(
        accession = name, score = score_name, contrast = group_name,
        n_first = effect[["n_x"]], n_second = effect[["n_y"]],
        mean_difference = effect[["difference"]], standardized_difference = effect[["d"]]
      )
    }
  }

  if (all(c("ra_convergence_score", "general_inflammation_score") %in% names(score_table))) {
    cell_names <- names(cell_sets)
    design_names <- c("general_inflammation_score", cell_names)
    complete <- complete.cases(score_table[, c("ra_convergence_score", design_names), with = FALSE])
    if (sum(complete) >= length(design_names) + 5) {
      base_model <- lm(ra_convergence_score ~ general_inflammation_score, data = score_table[complete])
      full_model <- lm(as.formula(paste("ra_convergence_score ~", paste(design_names, collapse = " + "))), data = score_table[complete])
      sensitivity_rows[[name]] <- data.table(accession = name, n = sum(complete),
                                              generic_only_r2 = summary(base_model)$r.squared,
                                              generic_plus_cell_state_r2 = summary(full_model)$r.squared,
                                              residual_ra_sd = sd(residuals(full_model)))
    }
  }

  if (all(c("RA", "OA") %in% unique(item$group))) {
    keep <- item$group %in% c("RA", "OA")
    adjusted_data <- score_table[keep]
    adjusted_data[, disease := factor(group, levels = c("OA", "RA"))]
    raw_model <- lm(ra_convergence_score ~ disease, data = adjusted_data)
    generic_model <- lm(ra_convergence_score ~ disease + general_inflammation_score, data = adjusted_data)
    cell_model_data <- adjusted_data[complete.cases(adjusted_data[, c("ra_convergence_score", "disease", "general_inflammation_score", names(cell_sets)), with = FALSE])]
    if (nrow(cell_model_data) >= 10) {
      cell_formula <- as.formula(paste("ra_convergence_score ~ disease + general_inflammation_score +", paste(names(cell_sets), collapse = " + ")))
      cell_model <- lm(cell_formula, data = cell_model_data)
      cell_beta <- coef(summary(cell_model))["diseaseRA", "Estimate"]
      cell_p <- coef(summary(cell_model))["diseaseRA", "Pr(>|t|)"]
    } else {
      cell_beta <- NA_real_; cell_p <- NA_real_
    }
    specificity_rows[[name]] <- data.table(
      accession = name, n = nrow(adjusted_data),
      raw_ra_vs_oa_difference = coef(summary(raw_model))["diseaseRA", "Estimate"],
      raw_ra_vs_oa_p = coef(summary(raw_model))["diseaseRA", "Pr(>|t|)"],
      adjusted_for_generic_difference = coef(summary(generic_model))["diseaseRA", "Estimate"],
      adjusted_for_generic_p = coef(summary(generic_model))["diseaseRA", "Pr(>|t|)"],
      adjusted_for_generic_and_cells_difference = cell_beta,
      adjusted_for_generic_and_cells_p = cell_p
    )
  }
}

score_table_all <- rbindlist(score_rows, fill = TRUE)
fwrite(score_table_all, file.path(result_dir, "continuous_projection_scores.csv"))
if (length(gsea_rows)) fwrite(rbindlist(gsea_rows, fill = TRUE), file.path(result_dir, "ranked_gsea_results.csv"))
fwrite(rbindlist(effect_rows, fill = TRUE), file.path(result_dir, "cross_cohort_effect_matrix.csv"))
if (length(sensitivity_rows)) fwrite(rbindlist(sensitivity_rows, fill = TRUE), file.path(result_dir, "cell_state_sensitivity.csv"))
if (length(specificity_rows)) fwrite(rbindlist(specificity_rows, fill = TRUE), file.path(result_dir, "ra_specificity_adjusted_effects.csv"))

capture.output(sessionInfo(), file = file.path(log_dir, "stage1_sessionInfo.txt"))
writeLines(c(
  paste("Completed", Sys.time()),
  paste("Cohorts", paste(names(cohorts), collapse = ",")),
  paste("Score rows", nrow(score_table_all)),
  "No subtype cutoff or formal clinical outcome model was run."
), file.path(log_dir, "stage1_execution_status.txt"))
print(platform_audit)
print(independence_audit)
