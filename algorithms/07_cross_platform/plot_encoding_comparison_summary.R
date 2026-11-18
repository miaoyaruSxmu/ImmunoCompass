options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

root <- gsub("\\\\", "/", Sys.getenv("IBI_TOOL_ROOT", unset = "C:/temp/ibi_tool"))
output_dir <- file.path(root, "Result", "revision_711_encoding_comparison_fixed390_20261110")
data <- read.csv(file.path(output_dir, "external_validation_primary_metrics.csv"))

method_labels <- c(
  absolute_expression = "Absolute expression",
  continuous_log_ratio = "Continuous log-ratio",
  binary_gene_pair_01 = "Binary gene-pair"
)
data$method_label <- factor(
  method_labels[data$method],
  levels = unname(method_labels)
)
data$cohort <- factor(data$cohort, levels = rev(unique(data$cohort)))

palette <- c(
  "Absolute expression" = "#737A80",
  "Continuous log-ratio" = "#C77C5D",
  "Binary gene-pair" = "#176B6B"
)

plot <- ggplot(data, aes(AUC, cohort, colour = method_label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#A9A9A9", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = AUC_CI_low, xmax = AUC_CI_high),
    orientation = "y", position = position_dodge(width = 0.58),
    width = 0.14, linewidth = 0.45
  ) +
  geom_point(
    position = position_dodge(width = 0.58), shape = 21,
    fill = "white", size = 2.2, stroke = 0.65
  ) +
  scale_colour_manual(values = palette) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Matched comparison of input representations",
    subtitle = "Final 711-sample development set and fixed 390-pair panel",
    x = "Direction-aware AUC (95% CI)", y = NULL, colour = NULL
  ) +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 8.5),
    plot.subtitle = element_text(size = 7.2, colour = "#454545"),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 6.8),
    legend.key.width = grid::unit(6, "mm"),
    plot.margin = margin(6, 8, 5, 6)
  )

prefix <- file.path(output_dir, "final711_three_input_external_AUC_summary")
grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = 150 / 25.4, height = 105 / 25.4, family = "Arial")
print(plot)
dev.off()
grDevices::png(paste0(prefix, ".png"), width = 150 / 25.4, height = 105 / 25.4,
               units = "in", res = 600, type = "cairo-png", family = "Arial")
print(plot)
dev.off()
