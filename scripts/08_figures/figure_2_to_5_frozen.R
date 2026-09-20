#!/usr/bin/env Rscript

# User-local library path intentionally omitted.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

project_dir <- normalizePath(getwd(), mustWork = TRUE)
submission_dir <- project_dir
figure_dir <- file.path(project_dir, "results/figure_exports")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

s2c <- file.path(project_dir, "04_结果/Stage2C独立验证")
s2b <- file.path(project_dir, "04_结果/Stage2B技术QC")

scores <- read.csv(file.path(s2c, "GSE283079_OA36_frozen_program_scores.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
corr <- read.csv(file.path(s2c, "GSE283079_program_correlations.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
rels <- read.csv(file.path(s2c, "OA_crosscohort_program_relationships.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(file.path(s2c, "OA_crosscohort_meta_analysis.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
models <- read.csv(file.path(s2c, "GSE283079_competing_models.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
qc <- read.csv(file.path(s2b, "GSE283079_analysis_sample_manifest_frozen.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(nrow(scores) == 36, nrow(models) == 3, length(unique(rels$cohort)) == 7)

pal <- c(
  ra = "#C44E3F", oa = "#2C7FB8", inflammation = "#C98720",
  apc = "#6F4AA8", myeloid = "#2A9D8F", neutral_dark = "#455A64",
  neutral_light = "#B8C2CC", zero = "#9AA5AE",
  positive = "#C44E3F", negative = "#3B6EA8"
)

theme_ir <- theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 9, face = "bold", colour = pal[["neutral_dark"]],
                              hjust = 0, margin = margin(b = 5)),
    axis.title = element_text(size = 8, colour = pal[["neutral_dark"]]),
    axis.text = element_text(size = 7, colour = pal[["neutral_dark"]]),
    axis.line = element_line(linewidth = 0.35, colour = pal[["neutral_dark"]]),
    axis.ticks = element_line(linewidth = 0.3, colour = pal[["neutral_dark"]]),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 7, 5, 5)
  )

panel_title <- function(letter, title) paste0(letter, "  ", title)

save_figure <- function(plot, stem, width_mm, height_mm) {
  base <- file.path(figure_dir, stem)
  ggsave(paste0(base, ".pdf"), plot, device = grDevices::pdf,
         width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
         bg = "white", useDingbats = FALSE)
  ggsave(paste0(base, ".svg"), plot, device = svglite::svglite,
         width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
         bg = "white")
  ggsave(paste0(base, ".png"), plot, device = ragg::agg_png,
         width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
         dpi = 150, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, device = ragg::agg_tiff,
         width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
}

## Figure 2: independent reconstruction and continuous projection.
qc$mapping_rate <- as.numeric(qc$mapping_rate)
qc <- qc[order(qc$mapping_rate), , drop = FALSE]
s2 <- scores[order(scores$ra_projection), , drop = FALSE]
s2$rank <- seq_len(nrow(s2))

p2a <- ggplot(qc, aes(reorder(SRR, mapping_rate), mapping_rate)) +
  geom_hline(yintercept = median(qc$mapping_rate), linetype = "dashed",
             linewidth = 0.4, colour = pal[["zero"]]) +
  geom_point(aes(colour = SRR == "SRR31542944"), size = 1.8) +
  geom_text(data = qc[qc$SRR == "SRR31542944", , drop = FALSE],
            aes(label = "SRR31542944 | 71.2%"), hjust = -0.05, vjust = -0.6,
            size = 2.35, colour = pal[["ra"]], family = "Helvetica") +
  scale_colour_manual(values = c(`TRUE` = pal[["ra"]], `FALSE` = pal[["neutral_dark"]]), guide = "none") +
  scale_x_discrete(labels = function(x) ifelse(x == "SRR31542944", "", "")) +
  coord_flip() + labs(title = panel_title("a", "GSE283079 reconstruction quality"),
                      x = NULL, y = "Salmon mapping rate (%)") + theme_ir +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p2b <- ggplot(s2, aes(rank, ra_projection)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_segment(aes(xend = rank, y = 0, yend = ra_projection),
               colour = pal[["ra"]], linewidth = 0.35) +
  geom_point(colour = pal[["ra"]], size = 1.8) +
  labs(title = panel_title("b", "Continuous RA-derived projection in OA"),
       x = "OA samples ordered by RA-derived projection", y = "RA-derived projection") + theme_ir

p2c <- ggplot(scores, aes(rank_score, ra_projection)) +
  geom_point(colour = pal[["neutral_dark"]], size = 1.7) +
  annotate("text", x = -Inf, y = Inf, label = "rho = 0.944\nn = 36",
           hjust = -0.08, vjust = 1.15, size = 2.55, family = "Helvetica") +
  labs(title = panel_title("c", "Scoring robustness"),
       x = "Rank-based sensitivity score", y = "Primary RA projection") + theme_ir

dprog <- c("general_inflammation", "mhc_ii_apc_primary", "mhc_ii_apc_generic",
            "myeloid_context", "interferon_state", "t_cell_context", "b_cell_context", "fibroblast_ecm")
d <- corr[match(dprog, corr$program), , drop = FALSE]
d$label <- c("General inflammation", "APC", "Generic MHC-II", "Myeloid", "IFN", "T-cell", "B-cell", "Fibroblast/ECM")
d$key <- c("inflammation", "apc", "apc", "myeloid", "neutral_dark", "neutral_dark", "neutral_dark", "neutral_dark")
d$label <- factor(d$label, levels = d$label[order(d$spearman_rho)])
d$y <- as.numeric(d$label)
p2d <- ggplot(d, aes(spearman_rho, y, colour = key)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, y = y, yend = y), linewidth = 0.5) +
  geom_segment(aes(x = ci_low, xend = ci_low, y = y - 0.09, yend = y + 0.09), linewidth = 0.5) +
  geom_segment(aes(x = ci_high, xend = ci_high, y = y - 0.09, yend = y + 0.09), linewidth = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(inflammation = pal[["inflammation"]], apc = pal[["apc"]],
                                 myeloid = pal[["myeloid"]], neutral_dark = pal[["neutral_dark"]]), guide = "none") +
  scale_y_continuous(breaks = seq_len(nrow(d)), labels = levels(d$label), expand = expansion(add = 0.35)) +
  labs(title = panel_title("d", "Programme relationships"),
       x = "Spearman rho (95% bootstrap CI)", y = NULL) + theme_ir

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_figure(fig2, "Figure2_IR_FINAL", 180, 120)

## Figure 3: genuine vector redraw from frozen cross-cohort relationships.
cohort_order <- c("GSE89408", "GSE55235", "GSE55457", "GSE55584", "GSE82107", "GSE206848", "GSE283079")
program_order <- c("general_inflammation", "mhc_ii_apc_primary", "mhc_ii_apc_generic",
                   "myeloid_context", "interferon_state", "t_cell_context", "b_cell_context", "fibroblast_ecm")
program_labels <- c(general_inflammation = "Inflammation", mhc_ii_apc_primary = "APC",
                    mhc_ii_apc_generic = "MHC-II", myeloid_context = "Myeloid",
                    interferon_state = "IFN", t_cell_context = "T-cell",
                    b_cell_context = "B-cell", fibroblast_ecm = "Fibroblast/ECM")
rr <- rels[rels$program %in% program_order, , drop = FALSE]
rr$cohort <- factor(rr$cohort, levels = rev(cohort_order))
rr$program <- factor(rr$program, levels = program_order)

p3a <- ggplot(rr, aes(program, cohort, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(is.na(spearman_rho), "", sprintf("%.2f", spearman_rho))),
            size = 2.25, family = "Helvetica", na.rm = FALSE) +
  scale_x_discrete(labels = program_labels) +
  scale_fill_gradient2(low = pal[["negative"]], mid = "white", high = pal[["positive"]],
                       midpoint = 0, limits = c(-1, 1), na.value = pal[["neutral_light"]],
                       name = "Spearman rho") +
  labs(title = panel_title("a", "Programme relationship matrix"), x = NULL, y = NULL) +
  theme_ir + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                   axis.text.y = element_text(size = 7), legend.position = "right",
                   legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
                   legend.key.height = grid::unit(22, "pt"))

general <- rr[rr$program == "general_inflammation", , drop = FALSE]
p3b <- ggplot(general, aes(spearman_rho, cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_point(aes(colour = cohort == "GSE283079"), size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = pal[["oa"]], `FALSE` = pal[["inflammation"]]), guide = "none") +
  scale_x_continuous(limits = c(-0.9, 0.9), breaks = c(-0.8, -0.4, 0, 0.4, 0.8)) +
  labs(title = panel_title("b", "Study-level general-inflammation\nrelationships"),
       x = "RA–inflammation rho", y = NULL) + theme_ir

pair <- rr[rr$program %in% c("mhc_ii_apc_primary", "mhc_ii_apc_generic"), , drop = FALSE]
pair$short <- ifelse(pair$program == "mhc_ii_apc_primary", "APC", "MHC-II")
p3c <- ggplot(pair, aes(spearman_rho, cohort, colour = short, shape = short)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(APC = pal[["apc"]], `MHC-II` = "#8A6BB1"), guide = "none") +
  scale_shape_manual(values = c(APC = 16, `MHC-II` = 15), name = NULL) +
  scale_x_continuous(limits = c(-0.9, 0.9), breaks = c(-0.8, -0.4, 0, 0.4, 0.8)) +
  labs(title = panel_title("c", "Within-cohort APC and MHC-II relationships"),
       x = "Spearman rho", y = NULL) + theme_ir +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.key = element_blank())

myeloid <- rr[rr$program == "myeloid_context", , drop = FALSE]
p3d <- ggplot(myeloid, aes(spearman_rho, cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_point(aes(colour = cohort == "GSE283079"), size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = pal[["oa"]], `FALSE` = pal[["myeloid"]]), guide = "none") +
  scale_x_continuous(limits = c(-0.9, 0.9), breaks = c(-0.8, -0.4, 0, 0.4, 0.8)) +
  labs(title = panel_title("d", "Study-level myeloid relationships"),
       x = "RA–myeloid rho", y = NULL) + theme_ir

fig3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(widths = c(1.15, 1), heights = c(1.15, 1)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_figure(fig3, "Figure3_IR_FINAL", 180, 135)

## Figure 4: cross-cohort effects with frozen estimates and palette semantics.
meta_key <- c("general_inflammation", "mhc_ii_apc_primary", "mhc_ii_apc_generic", "myeloid_context", "interferon_state")
meta <- meta[meta$program %in% meta_key, , drop = FALSE]
meta$programme <- factor(meta$program, levels = rev(meta_key))
meta$label <- c(general_inflammation = "RA–inflammation", mhc_ii_apc_primary = "RA–APC",
                mhc_ii_apc_generic = "RA–MHC-II", myeloid_context = "RA–myeloid",
                interferon_state = "RA–IFN")[meta$program]
meta$programme_label <- factor(meta$label, levels = rev(unname(c(
  general_inflammation = "RA–inflammation", mhc_ii_apc_primary = "RA–APC",
  mhc_ii_apc_generic = "RA–MHC-II", myeloid_context = "RA–myeloid", interferon_state = "RA–IFN"))))
meta$colour <- c(general_inflammation = pal[["inflammation"]], mhc_ii_apc_primary = pal[["apc"]],
                 mhc_ii_apc_generic = pal[["apc"]], myeloid_context = pal[["myeloid"]],
                 interferon_state = "#4E9BB5")[meta$program]
st <- rels[rels$program %in% meta_key, , drop = FALSE]
st$programme_label <- factor(meta$label[match(st$program, meta$program)], levels = levels(meta$programme_label))

p4a <- ggplot(st, aes(spearman_rho, programme_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_point(aes(colour = program), size = 1.75, alpha = 0.8) +
  scale_colour_manual(values = setNames(meta$colour, meta$program), guide = "none") +
  scale_x_continuous(limits = c(-0.9, 1), breaks = c(-0.5, 0, 0.5, 1)) +
  labs(title = panel_title("a", "Study-level effect estimates"), x = "Spearman rho", y = NULL) + theme_ir

meta$y <- as.numeric(meta$programme_label)
p4b <- ggplot(meta, aes(pooled_rho, y)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_segment(aes(x = ci_low, xend = ci_high, y = y, yend = y, colour = program), linewidth = 0.7) +
  geom_segment(aes(x = ci_low, xend = ci_low, y = y - 0.07, yend = y + 0.07, colour = program), linewidth = 0.7) +
  geom_segment(aes(x = ci_high, xend = ci_high, y = y - 0.07, yend = y + 0.07, colour = program), linewidth = 0.7) +
  geom_point(aes(colour = program), size = 2.8) +
  geom_text(aes(x = 0.76, label = sprintf("%.3f [%.3f, %.3f]", pooled_rho, ci_low, ci_high)),
            hjust = 0, size = 2.0, colour = pal[["neutral_dark"]], family = "Helvetica") +
  scale_colour_manual(values = setNames(meta$colour, meta$program), guide = "none") +
  scale_y_continuous(breaks = seq_len(nrow(meta)), labels = levels(meta$programme_label), expand = expansion(add = 0.35)) +
  scale_x_continuous(limits = c(-0.65, 1.55), breaks = c(-0.5, 0, 0.5, 1)) +
  labs(title = panel_title("b", "Random-effects pooled estimates"), x = "Pooled Spearman rho", y = NULL) + theme_ir +
  theme(plot.margin = margin(5, 55, 5, 5))

 p4c <- ggplot(meta, aes(i2, y)) +
  geom_segment(aes(x = 0, xend = i2, yend = y), linewidth = 0.8,
               colour = pal[["neutral_light"]]) +
  geom_point(aes(colour = program), size = 2.8) +
  geom_text(aes(x = pmin(i2 + 3, 64), label = sprintf("%.1f%%", i2)),
            hjust = 0, size = 2.25, colour = pal[["neutral_dark"]], family = "Helvetica") +
  scale_colour_manual(values = setNames(meta$colour, meta$program), guide = "none") +
  scale_x_continuous(limits = c(0, 78), breaks = c(0, 20, 40, 60)) +
  scale_y_continuous(breaks = seq_len(nrow(meta)), labels = levels(meta$programme_label), expand = expansion(add = 0.35)) +
  labs(title = panel_title("c", "Cross-cohort stability\nand heterogeneity"), x = "I² (%)", y = NULL) + theme_ir +
  theme(plot.margin = margin(5, 12, 5, 5))

fig4 <- (p4a / (p4b | p4c)) +
  plot_layout(heights = c(1.08, 1), widths = c(1.25, 0.75)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_figure(fig4, "Figure4_IR_FINAL", 180, 125)

## Figure 5: evidence-boundary figure; no new claim layer.
model_general <- models[models$model == "general_only", , drop = FALSE]
model_apc <- models[models$model == "general_plus_apc", , drop = FALSE]
model_full <- models[models$model == "general_apc_cell_states", , drop = FALSE]
model_labels <- c(general_only = "Inflammation only", general_plus_apc = "Inflammation + APC",
                  general_apc_cell_states = "Inflammation + APC + cell states")
explained <- data.frame(label = factor(unname(model_labels), levels = rev(unname(model_labels))),
                        estimate = c(model_general$adjusted_r_squared, model_apc$adjusted_r_squared,
                                     model_full$adjusted_r_squared))
residual <- data.frame(label = factor(c("Original RA score", "After inflammation", "After inflammation + APC"),
                                      levels = rev(c("Original RA score", "After inflammation", "After inflammation + APC"))),
                       estimate = c(sd(scores$ra_projection), model_general$residual_sd, model_apc$residual_sd))
coefficients <- data.frame(
  predictor = factor(c("Inflammation | inflammation-only model",
                       "Inflammation | inflammation + APC model",
                       "APC | inflammation + APC model",
                       "Myeloid | full model"),
                     levels = rev(c("Inflammation | inflammation-only model",
                                    "Inflammation | inflammation + APC model",
                                    "APC | inflammation + APC model",
                                    "Myeloid | full model"))),
  estimate = c(model_full$beta_general_inflammation, model_apc$beta_general_inflammation,
               model_apc$beta_mhc_ii_apc_primary, model_full$beta_myeloid_context),
  lower = c(model_full$ci_low_general_inflammation, model_apc$ci_low_general_inflammation,
            model_apc$ci_low_mhc_ii_apc_primary, model_full$ci_low_myeloid_context),
  upper = c(model_full$ci_high_general_inflammation, model_apc$ci_high_general_inflammation,
            model_apc$ci_high_mhc_ii_apc_primary, model_full$ci_high_myeloid_context))

p5a <- ggplot(explained, aes(estimate, label)) +
  geom_segment(aes(x = 0, xend = estimate, yend = label), colour = pal[["neutral_light"]], linewidth = 0.8) +
  geom_point(size = 2.8, colour = pal[["oa"]]) +
  geom_text(aes(x = estimate + 0.018, label = sprintf("%.3f", estimate)), hjust = 0,
            size = 2.45, colour = pal[["neutral_dark"]], family = "Helvetica") +
  scale_x_continuous(limits = c(0, 0.44), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  labs(title = panel_title("a", "Explained variance"), x = "Adjusted R²", y = NULL) + theme_ir

p5b <- ggplot(residual, aes(estimate, label)) +
  geom_segment(aes(x = 0, xend = estimate, yend = label), colour = pal[["neutral_light"]], linewidth = 0.8) +
  geom_point(size = 2.8, colour = pal[["neutral_dark"]]) +
  geom_text(aes(x = estimate + 0.035, label = sprintf("%.3f", estimate)), hjust = 0,
            size = 2.45, colour = pal[["neutral_dark"]], family = "Helvetica") +
  scale_x_continuous(limits = c(0, 0.98), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  labs(title = panel_title("b", "Residual variation"), x = "Residual standard deviation", y = NULL) + theme_ir

coefficients$y <- as.numeric(coefficients$predictor)
p5c <- ggplot(coefficients, aes(estimate, y)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = pal[["zero"]]) +
  geom_segment(aes(x = lower, xend = upper, y = y, yend = y), linewidth = 0.65) +
  geom_segment(aes(x = lower, xend = lower, y = y - 0.08, yend = y + 0.08), linewidth = 0.65) +
  geom_segment(aes(x = upper, xend = upper, y = y - 0.08, yend = y + 0.08), linewidth = 0.65) +
  geom_point(size = 2.7, colour = pal[["neutral_dark"]]) +
  scale_y_continuous(breaks = seq_len(nrow(coefficients)), labels = levels(coefficients$predictor), expand = expansion(add = 0.35)) +
  scale_x_continuous(limits = c(-1.7, 4.8), breaks = c(-1, 0, 1, 2), expand = c(0, 0)) +
  labs(title = panel_title("c", "Frozen model coefficients"), x = "Coefficient (95% CI)", y = NULL) + theme_ir

ladder <- data.frame(y = c(4, 3, 2, 1),
                     claim = c("Measurable RA-derived resemblance", "Partial inflammatory convergence",
                               "Stable multicomponent immune architecture", "Stable RA-like subtype identity"),
                     status = c("SUPPORTED", "PARTIAL / SUPPORTED", "NOT ESTABLISHED", "NOT ESTABLISHED"),
                     colour = c(pal[["ra"]], pal[["inflammation"]], pal[["neutral_dark"]], pal[["neutral_dark"]]))
p5d <- ggplot(ladder, aes(x = 1, y = y)) +
  geom_segment(data = ladder[1:3, ], aes(x = 1, xend = 1, y = y - 0.22, yend = y - 0.78),
               colour = pal[["neutral_light"]], linewidth = 0.65,
               arrow = arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
  geom_point(aes(colour = claim), size = 3.2) +
  geom_text(aes(x = 1.10, label = claim), hjust = 0, size = 2.4,
            colour = pal[["neutral_dark"]], family = "Helvetica") +
  geom_text(aes(x = 1.10, y = y - 0.26, label = status), hjust = 0, size = 1.95,
            colour = ladder$colour, family = "Helvetica", fontface = "bold") +
  scale_colour_manual(values = setNames(ladder$colour, ladder$claim), guide = "none") +
  scale_x_continuous(limits = c(0.85, 2.6), breaks = NULL) +
  scale_y_continuous(limits = c(0.05, 4.5), breaks = NULL) +
  labs(title = panel_title("d", "Evidence boundary"), x = NULL, y = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 9, face = "bold", colour = pal[["neutral_dark"]], hjust = 0,
                                  margin = margin(b = 5)), plot.margin = margin(5, 5, 5, 5))

fig5 <- (p5a | p5b) / (p5c | p5d) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_figure(fig5, "Figure5_IR_FINAL", 180, 135)

writeLines(capture.output(sessionInfo()), file.path(submission_dir, "figure_build_session_info.txt"))
cat("Built Figures 2–5 in", normalizePath(figure_dir), "\n")
