#!/usr/bin/env Rscript

project_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(project_dir, "results/frozen_summary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s2c <- file.path(project_dir, "04_结果/Stage2C独立验证")
rels <- read.csv(file.path(s2c, "OA_crosscohort_program_relationships.csv"), stringsAsFactors = FALSE)
meta <- read.csv(file.path(s2c, "OA_crosscohort_meta_analysis.csv"), stringsAsFactors = FALSE)
keys <- c("general_inflammation", "mhc_ii_apc_primary", "mhc_ii_apc_generic", "myeloid_context", "interferon_state")
labels <- c(general_inflammation = "RA–inflammation", mhc_ii_apc_primary = "RA–APC",
            mhc_ii_apc_generic = "RA–MHC-II", myeloid_context = "RA–myeloid",
            interferon_state = "RA–IFN")

reml_tau2 <- function(y, v) {
  nll <- function(tau2) {
    wi <- 1 / (v + tau2)
    mu <- sum(wi * y) / sum(wi)
    0.5 * (sum(log(v + tau2)) + log(sum(wi)) + sum(wi * (y - mu)^2))
  }
  upper <- max(10, max(y^2) + 1)
  opt <- optimize(nll, interval = c(0, upper), tol = 1e-12)
  max(0, opt$minimum)
}

rows <- lapply(keys, function(key) {
  x <- rels[rels$program == key & is.finite(rels$spearman_rho) & rels$n > 3, , drop = FALSE]
  y <- atanh(x$spearman_rho)
  v <- 1 / (x$n - 3)
  k <- length(y)
  tau2 <- reml_tau2(y, v)
  w <- 1 / (v + tau2)
  mu <- sum(w * y) / sum(w)
  q <- sum(w * (y - mu)^2)
  hk_var <- q / ((k - 1) * sum(w))
  crit <- qt(0.975, df = k - 1)
  lo_z <- mu - crit * sqrt(hk_var)
  hi_z <- mu + crit * sqrt(hk_var)
  pooled <- tanh(mu)
  lo <- tanh(lo_z)
  hi <- tanh(hi_z)
  main <- meta[meta$program == key, , drop = FALSE]
  data.frame(
    relationship = unname(labels[key]), program = key, k = k,
    main_DL_pooled_rho = main$pooled_rho, main_DL_ci_low = main$ci_low,
    main_DL_ci_high = main$ci_high, REML_tau2 = tau2,
    REML_HK_pooled_rho = pooled, REML_HK_ci_low = lo, REML_HK_ci_high = hi,
    HK_df = k - 1, HK_q = q,
    stringsAsFactors = FALSE
  )
})
tab <- do.call(rbind, rows)
write.table(tab, file.path(out_dir, "Table_S9_REML_HK_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

qual <- paste(
  "# Supplementary Methods update", "",
  "As a sensitivity analysis, cross-cohort Spearman correlations were transformed to Fisher z values and pooled using restricted maximum likelihood (REML) estimation of between-cohort variance. Hartung–Knapp confidence intervals were calculated on the Fisher-z scale with k−1 degrees of freedom and back-transformed with tanh. The DerSimonian–Laird analysis remained the prespecified main analysis.", "",
  "# Supplementary Results update", "",
  "The REML/Hartung–Knapp sensitivity analysis is reported in Table S9 and was not used to replace the main DerSimonian–Laird estimates. The general-inflammation relationship remained the most comparatively portable among the five relationships by point estimate and low estimated heterogeneity, whereas APC/MHC-II, myeloid and interferon estimates remained context dependent or imprecise. Hartung–Knapp intervals were wider than the main-analysis intervals for several relationships, as expected with seven cohorts; the qualitative evidence boundary was unchanged.", sep = "\n")
writeLines(qual, file.path(out_dir, "Supplementary_REML_HK_update.md"), useBytes = TRUE)
print(tab)
