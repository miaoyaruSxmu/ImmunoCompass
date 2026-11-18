#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

root <- Sys.getenv("IBI_TOOL_ROOT", unset = getwd())
analysis_dir <- file.path(
  root, "Result", "revision_781_BAS_CAS_IAS_sensitivity_top390_20261028"
)
expression_file <- file.path(root, "Result", "validation", "GSE186143_bulk_TPM.txt")
sample_file <- file.path(analysis_dir, "GSE186143_top390_IBI_and_grade.csv")
coefficient_file <- file.path(analysis_dir, "conventional_model_coefficients.csv")
fit_file <- file.path(analysis_dir, "conventional_model_fit_metrics.csv")
comparison_file <- file.path(analysis_dir, "conventional_model_comparison_summary.csv")
nested_file <- file.path(analysis_dir, "nested_model_comparison_summary.csv")
output_dir <- file.path(
  root, "Result", "revision_comment16_existing_methods_comparison_20261113"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  expression_file, sample_file, coefficient_file, fit_file,
  comparison_file, nested_file
)
if (!all(file.exists(required_files))) {
  stop("Missing input file(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

expression_all <- read.delim(
  expression_file, row.names = 1, check.names = FALSE,
  quote = "", comment.char = ""
)
expression_all <- as.matrix(expression_all)
storage.mode(expression_all) <- "double"
sample_data <- read.csv(sample_file, check.names = FALSE)
coefficients <- read.csv(coefficient_file, check.names = FALSE)
fit_metrics <- read.csv(fit_file, check.names = FALSE)
comparisons <- read.csv(comparison_file, check.names = FALSE)
nested <- read.csv(nested_file, check.names = FALSE)

stopifnot(
  nrow(sample_data) == 60L,
  all(sample_data$sample_id %in% colnames(expression_all)),
  nrow(coefficients) == 20647L,
  nrow(fit_metrics) == 20647L,
  nrow(comparisons) == 3L,
  nrow(nested) == 2L
)

expression <- expression_all[, sample_data$sample_id, drop = FALSE]
expression_filter <- rowMeans(expression > 0.5, na.rm = TRUE) > 0.20
expression_log <- log2(expression[expression_filter, , drop = FALSE] + 1)
gene_sd <- apply(expression_log, 1L, sd, na.rm = TRUE)
expression_log <- expression_log[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
expression_z <- t(scale(t(expression_log)))

if (!identical(rownames(expression_z), coefficients$gene)) {
  index <- match(coefficients$gene, rownames(expression_z))
  if (anyNA(index)) stop("Modeled genes could not be aligned to the expression matrix")
  expression_z <- expression_z[index, , drop = FALSE]
}

z_score <- function(x) as.numeric(scale(as.numeric(x)))
ibi_z <- z_score(sample_data$IBI)
grade_z <- z_score(sample_data$grade_number)
interaction_z <- z_score(ibi_z * grade_z)
design_full <- cbind(
  Intercept = 1, IBI = ibi_z, Grade = grade_z,
  Interaction = interaction_z
)

rss_for_design <- function(y_gene_by_sample, design) {
  response <- t(y_gene_by_sample)
  fit <- lm.fit(x = design, y = response)
  colSums(fit$residuals^2)
}

rss_full <- rss_for_design(expression_z, design_full)
rss_without_ibi <- rss_for_design(
  expression_z, design_full[, c("Intercept", "Grade", "Interaction"), drop = FALSE]
)
rss_without_grade <- rss_for_design(
  expression_z, design_full[, c("Intercept", "IBI", "Interaction"), drop = FALSE]
)
rss_without_interaction <- rss_for_design(
  expression_z, design_full[, c("Intercept", "IBI", "Grade"), drop = FALSE]
)

partial_r2 <- function(rss_reduced, rss_complete) {
  value <- (rss_reduced - rss_complete) / rss_reduced
  pmax(0, pmin(1, value))
}

gene_level <- data.frame(
  gene = coefficients$gene,
  BAS = abs(coefficients$beta_ibi_full),
  CAS = abs(coefficients$beta_grade_full),
  IAS = abs(coefficients$beta_interaction_full),
  partial_R2_BAS = partial_r2(rss_without_ibi, rss_full),
  partial_R2_CAS = partial_r2(rss_without_grade, rss_full),
  partial_R2_IAS = partial_r2(rss_without_interaction, rss_full),
  stringsAsFactors = FALSE
)

max_component <- function(x) {
  c("BAS", "CAS", "IAS")[max.col(as.matrix(x), ties.method = "first")]
}

dominance_component <- function(x, margin = 1.5) {
  matrix_value <- as.matrix(x)
  best <- max.col(matrix_value, ties.method = "first")
  ordered <- t(apply(matrix_value, 1L, sort, decreasing = TRUE))
  output <- rep("None", nrow(matrix_value))
  keep <- ordered[, 1L] >= margin * ordered[, 2L]
  output[keep] <- c("BAS", "CAS", "IAS")[best[keep]]
  output
}

gene_level$coefficient_max <- max_component(gene_level[, c("BAS", "CAS", "IAS")])
gene_level$partial_R2_max <- max_component(
  gene_level[, c("partial_R2_BAS", "partial_R2_CAS", "partial_R2_IAS")]
)
gene_level$coefficient_dominance_1_5 <- dominance_component(
  gene_level[, c("BAS", "CAS", "IAS")]
)
gene_level$partial_R2_dominance_1_5 <- dominance_component(
  gene_level[, c("partial_R2_BAS", "partial_R2_CAS", "partial_R2_IAS")]
)

cohen_kappa <- function(a, b, levels) {
  tab <- table(factor(a, levels = levels), factor(b, levels = levels))
  observed <- sum(diag(tab)) / sum(tab)
  expected <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  (observed - expected) / (1 - expected)
}

agreement <- data.frame(
  classification = c("Maximum component", "1.5-fold dominance"),
  agreement = c(
    mean(gene_level$coefficient_max == gene_level$partial_R2_max),
    mean(gene_level$coefficient_dominance_1_5 == gene_level$partial_R2_dominance_1_5)
  ),
  cohen_kappa = c(
    cohen_kappa(
      gene_level$coefficient_max, gene_level$partial_R2_max,
      c("BAS", "CAS", "IAS")
    ),
    cohen_kappa(
      gene_level$coefficient_dominance_1_5,
      gene_level$partial_R2_dominance_1_5,
      c("BAS", "CAS", "IAS", "None")
    )
  )
)

partial_long <- rbind(
  data.frame(component = "BAS", partial_R2 = gene_level$partial_R2_BAS),
  data.frame(component = "CAS", partial_R2 = gene_level$partial_R2_CAS),
  data.frame(component = "IAS", partial_R2 = gene_level$partial_R2_IAS)
)
partial_long$component <- factor(partial_long$component, levels = c("BAS", "CAS", "IAS"))

partial_summary <- do.call(rbind, lapply(split(partial_long$partial_R2, partial_long$component), function(x) {
  data.frame(
    median = median(x),
    q1 = unname(quantile(x, 0.25)),
    q3 = unname(quantile(x, 0.75))
  )
}))
partial_summary$component <- rownames(partial_summary)
rownames(partial_summary) <- NULL

write.csv(gene_level, file.path(output_dir, "source_data_gene_level_partial_R2.csv"), row.names = FALSE)
write.csv(agreement, file.path(output_dir, "coefficient_partial_R2_agreement.csv"), row.names = FALSE)
write.csv(partial_summary, file.path(output_dir, "partial_R2_summary.csv"), row.names = FALSE)

comparison_map <- c(
  "Condition-only Grade vs IBI-adjusted Grade" = "Condition only vs IBI adjusted",
  "Condition-only Grade vs full-model CAS" = "Condition only vs full model",
  "IBI-adjusted Grade vs full-model CAS" = "IBI adjusted vs full model"
)

scatter_specs <- list(
  list(
    x = "beta_grade_condition_only", y = "beta_grade_ibi_adjusted",
    xlab = "Grade coefficient: condition-only model",
    ylab = "Grade coefficient: IBI-adjusted model",
    key = "Condition-only Grade vs IBI-adjusted Grade"
  ),
  list(
    x = "beta_grade_condition_only", y = "beta_grade_full",
    xlab = "Grade coefficient: condition-only model",
    ylab = "Grade coefficient: full interaction model",
    key = "Condition-only Grade vs full-model CAS"
  ),
  list(
    x = "beta_grade_ibi_adjusted", y = "beta_grade_full",
    xlab = "Grade coefficient: IBI-adjusted model",
    ylab = "Grade coefficient: full interaction model",
    key = "IBI-adjusted Grade vs full-model CAS"
  )
)

palette <- c(
  neutral = "#777777", blue = "#3B6FB6", teal = "#2A9D8F",
  orange = "#D9822B", red = "#C94C4C", green = "#4C956C"
)

theme_response <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = 6.2, colour = "black"),
      axis.title = element_text(size = 7),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(0, 1),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

make_scatter <- function(spec) {
  summary_row <- comparisons[comparisons$comparison == spec$key, , drop = FALSE]
  annotation <- sprintf(
    "Spearman rho = %.3f\nTop-100 Jaccard = %.3f",
    summary_row$spearman_signed_beta,
    summary_row$top100_jaccard
  )
  x <- coefficients[[spec$x]]
  y <- coefficients[[spec$y]]
  bounds <- range(c(x, y), finite = TRUE)
  padding <- diff(bounds) * 0.04
  ggplot(coefficients, aes(x = .data[[spec$x]], y = .data[[spec$y]])) +
    geom_point(size = 0.35, alpha = 0.12, colour = palette[["blue"]], stroke = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#555555") +
    annotate(
      "text", x = bounds[1] + padding, y = bounds[2] - padding,
      label = annotation, hjust = 0, vjust = 1, size = 2.35,
      family = "Arial"
    ) +
    coord_equal(xlim = bounds, ylim = bounds, expand = FALSE) +
    labs(x = spec$xlab, y = spec$ylab) +
    theme_response()
}

p_a <- make_scatter(scatter_specs[[1]])
p_b <- make_scatter(scatter_specs[[2]])
p_c <- make_scatter(scatter_specs[[3]])

fit_long <- rbind(
  data.frame(model = "Condition only", adjusted_R2 = fit_metrics$adjusted_R2_condition_only),
  data.frame(model = "IBI adjusted", adjusted_R2 = fit_metrics$adjusted_R2_IBI_adjusted),
  data.frame(model = "Full interaction", adjusted_R2 = fit_metrics$adjusted_R2_full_interaction)
)
fit_long$model <- factor(
  fit_long$model,
  levels = c("Condition only", "IBI adjusted", "Full interaction")
)

p_d <- ggplot(fit_long, aes(model, adjusted_R2, fill = model)) +
  geom_violin(width = 0.82, linewidth = 0.25, trim = TRUE, colour = NA, alpha = 0.78) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  scale_fill_manual(
    values = c(
      "Condition only" = "#B8B8B8",
      "IBI adjusted" = palette[["blue"]],
      "Full interaction" = palette[["teal"]]
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = "Adjusted R-squared") +
  theme_response() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

delta_long <- rbind(
  data.frame(term = "Add IBI", delta_R2 = fit_metrics$delta_R2_add_IBI),
  data.frame(term = "Add interaction", delta_R2 = fit_metrics$delta_R2_add_interaction)
)
delta_long$term <- factor(delta_long$term, levels = c("Add IBI", "Add interaction"))
nested$term <- factor(
  ifelse(nested$added_term == "IBI", "Add IBI", "Add interaction"),
  levels = levels(delta_long$term)
)
delta_max <- quantile(delta_long$delta_R2, 0.995, na.rm = TRUE)
nested$label_y <- delta_max * c(0.98, 0.82)
nested$label_hjust <- c(0, 1)
nested$label <- sprintf(
  "Median change in R-squared: %.4f\nBH FDR < 0.05: %s genes",
  nested$median_delta_R2,
  comma(nested$n_genes_FDR_lt_0_05)
)

p_e <- ggplot(delta_long, aes(term, delta_R2, fill = term)) +
  geom_violin(width = 0.82, linewidth = 0.25, trim = TRUE, colour = NA, alpha = 0.78) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_text(
    data = nested,
    aes(x = term, y = label_y, label = label, hjust = label_hjust),
    inherit.aes = FALSE, size = 1.9, lineheight = 0.95, family = "Arial"
  ) +
  scale_fill_manual(
    values = c("Add IBI" = palette[["blue"]], "Add interaction" = palette[["teal"]]),
    guide = "none"
  ) +
  coord_cartesian(ylim = c(0, delta_max)) +
  labs(x = NULL, y = "Incremental R-squared") +
  theme_response()

partial_summary$label <- sprintf("Median %.4f", partial_summary$median)
partial_summary$component <- factor(partial_summary$component, levels = c("BAS", "CAS", "IAS"))
partial_label_y <- quantile(partial_long$partial_R2, 0.995, na.rm = TRUE)
partial_summary$label_y <- partial_label_y * c(0.96, 0.81, 0.66)
agreement_label <- sprintf(
  "Max-component agreement: %.1f%%\nCohen's kappa: %.3f",
  100 * agreement$agreement[agreement$classification == "Maximum component"],
  agreement$cohen_kappa[agreement$classification == "Maximum component"]
)

p_f <- ggplot(partial_long, aes(component, partial_R2, fill = component)) +
  geom_violin(width = 0.82, linewidth = 0.25, trim = TRUE, colour = NA, alpha = 0.78) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_text(
    data = partial_summary,
    aes(x = component, y = label_y, label = label),
    inherit.aes = FALSE, size = 2.1, family = "Arial"
  ) +
  annotate(
    "text", x = 2, y = partial_label_y * 1.08,
    label = agreement_label, size = 1.9, lineheight = 0.95, family = "Arial"
  ) +
  scale_fill_manual(
    values = c("BAS" = palette[["blue"]], "CAS" = palette[["orange"]], "IAS" = palette[["green"]]),
    guide = "none"
  ) +
  coord_cartesian(ylim = c(0, partial_label_y * 1.14)) +
  labs(x = NULL, y = "Drop-one partial R-squared") +
  theme_response()

figure <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_annotation(tag_levels = "A")

base <- file.path(output_dir, "Response_Figure_R2_16_existing_methods_comparison")
width_mm = 183
height_mm = 128
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()

grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(figure)
dev.off()

temp_tiff <- tempfile(pattern = "response_r2_16_", fileext = ".tiff")
temp_png <- tempfile(pattern = "response_r2_16_", fileext = ".png")

ragg::agg_tiff(
  temp_tiff, width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure)
dev.off()

ragg::agg_png(
  temp_png, width = width_in, height = height_in,
  units = "in", res = 300
)
print(figure)
dev.off()

if (!file.copy(temp_tiff, paste0(base, ".tiff"), overwrite = TRUE)) {
  stop("Failed to copy the TIFF export into the output directory")
}
if (!file.copy(temp_png, paste0(base, ".png"), overwrite = TRUE)) {
  stop("Failed to copy the PNG export into the output directory")
}
unlink(c(temp_tiff, temp_png))

qa_notes <- c(
  "Core conclusion: BAS/CAS/IAS is an interpretation and prioritization layer over conventional limma regression, and its effect-source attribution can be expressed on a unique-variance scale.",
  "Archetype: quantitative grid.",
  "Backend: R only (ggplot2 and patchwork).",
  "Analysis unit: 20,647 genes from 60 patients with complete Grade data.",
  "No modeled genes or patients were excluded relative to the archived standardized analysis.",
  "Panels A-C use signed standardized Grade coefficients; annotations report Spearman correlation and Top-100 Jaccard overlap.",
  "Panels D-E use archived OLS fit metrics and nested-model F tests with BH correction across genes.",
  "Panel F recomputes ordinary least-squares drop-one partial R-squared from the standardized expression and design matrices; it does not use limma moderated t statistics.",
  "Exports: editable SVG and PDF, 600-dpi TIFF, and 300-dpi PNG preview."
)
writeLines(qa_notes, file.path(output_dir, "figure_qa_notes.txt"), useBytes = TRUE)

cat("Created response figure bundle in: ", output_dir, "\n", sep = "")
print(partial_summary[, c("component", "median", "q1", "q3")])
print(agreement)
