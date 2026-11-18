args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript comment16_partial_R2_functional_validation.R <input_dir> <output_dir>")
}

input_dir <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_file <- file.path(input_dir, "standardized_BAS_CAS_IAS_gene_level.csv")
go_file <- file.path(input_dir, "GO_BP_threshold_enrichment.csv")
dat <- read.csv(gene_file, check.names = FALSE, stringsAsFactors = FALSE)
go <- read.csv(go_file, check.names = FALSE, stringsAsFactors = FALSE)

# Full model: standardized expression ~ IBI_z + Grade_z + Interaction_z.
# With n = 60 and four fitted coefficients (including the intercept), df = 56.
# For a one-df term, drop-one partial R2 is exactly t^2 / (t^2 + residual df).
n <- 60L
residual_df <- n - 4L
dat$full_model_R2 <- NA_real_
fit_file <- file.path(input_dir, "conventional_model_fit_metrics.csv")
if (file.exists(fit_file)) {
  fit <- read.csv(fit_file, check.names = FALSE, stringsAsFactors = FALSE)
  dat$full_model_R2 <- fit$R2_full_interaction[match(dat$gene, fit$gene)]
}
dat$partial_R2_IBI <- dat$t_BAS^2 / (dat$t_BAS^2 + residual_df)
dat$partial_R2_Grade <- dat$t_CAS^2 / (dat$t_CAS^2 + residual_df)
dat$partial_R2_Interaction <- dat$t_IAS^2 / (dat$t_IAS^2 + residual_df)

max_class <- function(a, b, c) {
  labs <- c("BAS", "CAS", "IAS")
  m <- cbind(a, b, c)
  labs[max.col(m, ties.method = "first")]
}

threshold_class <- function(a, b, c, fold = 1.5) {
  m <- cbind(a, b, c)
  best <- max.col(m, ties.method = "first")
  sm <- t(apply(m, 1, sort, decreasing = TRUE))
  out <- rep("None", nrow(m))
  keep <- sm[, 1] >= fold * sm[, 2]
  out[keep] <- c("BAS", "CAS", "IAS")[best[keep]]
  out
}

dat$coefficient_max_component <- max_class(dat$BAS, dat$CAS, dat$IAS)
dat$partial_R2_max_component <- max_class(dat$partial_R2_IBI, dat$partial_R2_Grade,
                                           dat$partial_R2_Interaction)
dat$coefficient_dominance_1_5 <- dat$dominance_1_5
dat$partial_R2_dominance_1_5 <- threshold_class(dat$partial_R2_IBI,
                                                 dat$partial_R2_Grade,
                                                 dat$partial_R2_Interaction)

top_stats <- function(beta, pr2, component) {
  rho <- suppressWarnings(cor(beta, pr2, method = "spearman", use = "complete.obs"))
  rows <- lapply(c(50L, 100L), function(k) {
    a <- head(order(beta, decreasing = TRUE), k)
    b <- head(order(pr2, decreasing = TRUE), k)
    inter <- length(intersect(a, b))
    data.frame(component = component, statistic = paste0("Top", k),
               spearman_rho = rho, overlap_n = inter,
               overlap_fraction = inter / k,
               jaccard = inter / length(union(a, b)))
  })
  do.call(rbind, rows)
}

summary_top <- rbind(
  top_stats(dat$BAS, dat$partial_R2_IBI, "BAS"),
  top_stats(dat$CAS, dat$partial_R2_Grade, "CAS"),
  top_stats(dat$IAS, dat$partial_R2_Interaction, "IAS")
)

cohen_kappa <- function(a, b, levels) {
  a <- factor(a, levels = levels)
  b <- factor(b, levels = levels)
  tab <- table(a, b)
  po <- sum(diag(tab)) / sum(tab)
  pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  (po - pe) / (1 - pe)
}

agreement <- rbind(
  data.frame(classification = "maximum component",
             agreement = mean(dat$coefficient_max_component == dat$partial_R2_max_component),
             cohen_kappa = cohen_kappa(dat$coefficient_max_component,
                                       dat$partial_R2_max_component,
                                       c("BAS", "CAS", "IAS"))),
  data.frame(classification = "1.5-fold dominance",
             agreement = mean(dat$coefficient_dominance_1_5 == dat$partial_R2_dominance_1_5),
             cohen_kappa = cohen_kappa(dat$coefficient_dominance_1_5,
                                       dat$partial_R2_dominance_1_5,
                                       c("BAS", "CAS", "IAS", "None")))
)

write.csv(dat, file.path(output_dir, "Comment16_gene_level_partial_R2.csv"), row.names = FALSE)
write.csv(summary_top, file.path(output_dir, "Comment16_partial_R2_concordance_summary.csv"), row.names = FALSE)
write.csv(agreement, file.path(output_dir, "Comment16_dominance_agreement_summary.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(table(dat$coefficient_max_component,
                                     dat$partial_R2_max_component)),
          file.path(output_dir, "Comment16_max_component_confusion_matrix.csv"))
write.csv(as.data.frame.matrix(table(dat$coefficient_dominance_1_5,
                                     dat$partial_R2_dominance_1_5)),
          file.path(output_dir, "Comment16_1_5fold_confusion_matrix.csv"))

# Existing functional enrichment is reorganized for the response figure.
needed <- c("threshold", "component", "Description", "p.adjust")
if (!all(needed %in% names(go))) {
  stop("GO enrichment file does not contain the expected columns: ",
       paste(setdiff(needed, names(go)), collapse = ", "))
}
go15 <- go[go$threshold == 1.5 & go$p.adjust < 0.05 & go$component %in% c("BAS", "CAS", "IAS"), ]
go15 <- go15[order(go15$component, go15$p.adjust), ]
go_plot <- do.call(rbind, lapply(split(go15, go15$component), function(x) head(x, 5)))
write.csv(go_plot, file.path(output_dir, "Comment16_existing_GO_enrichment_reorganized.csv"), row.names = FALSE)

component_cols <- c(BAS = "#377EB8", CAS = "#E41A1C", IAS = "#4DAF4A")
draw_scatter <- function(x, y, component, letter) {
  smoothScatter(x, y, nrpoints = 500, colramp = colorRampPalette(c("white", component_cols[component])),
                xlab = paste0(component, ": |standardized beta|"),
                ylab = "Drop-one partial R-squared", main = paste0(letter, "  ", component),
                cex.lab = 1.25, cex.axis = 1.10, cex.main = 1.35,
                font.main = 2)
  rho <- cor(x, y, method = "spearman", use = "complete.obs")
  legend("topleft", legend = sprintf("Spearman rho = %.3f", rho), bty = "n", cex = 1.05)
}

draw_figure <- function(device_fun, filename) {
  device_fun(filename)
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE), heights = c(1, 1.05))
  par(mar = c(5.0, 5.2, 3.0, 1.2), mgp = c(3.0, 0.9, 0))
  draw_scatter(dat$BAS, dat$partial_R2_IBI, "BAS", "A")
  draw_scatter(dat$CAS, dat$partial_R2_Grade, "CAS", "B")
  draw_scatter(dat$IAS, dat$partial_R2_Interaction, "IAS", "C")
  par(mar = c(6.3, 11.0, 3.0, 1.2), mgp = c(3.0, 0.9, 0))
  if (nrow(go_plot) > 0) {
    concise_labels <- c(
      "cellular defense response" = "Cellular defense response",
      "cristae formation" = "Cristae formation",
      "inner mitochondrial membrane organization" = "Inner membrane organization",
      "proton motive force-driven mitochondrial ATP synthesis" = "Mitochondrial ATP synthesis",
      "energy coupled proton transport, down electrochemical gradient" = "Proton transport down gradient",
      "proton motive force-driven ATP synthesis" = "Proton-driven ATP synthesis"
    )
    go_plot$label <- unname(concise_labels[go_plot$Description])
    go_plot$label[is.na(go_plot$label)] <- substr(go_plot$Description[is.na(go_plot$label)], 1, 34)
    go_plot <- go_plot[order(go_plot$component, -log10(go_plot$p.adjust)), ]
    vals <- -log10(go_plot$p.adjust)
    yy <- seq_along(vals)
    barplot(vals, names.arg = go_plot$label, horiz = TRUE, las = 1,
            col = component_cols[go_plot$component], border = NA,
            xlab = expression(-log[10](BH~FDR)), main = "D  Existing GO-BP enrichment",
            cex.names = 0.92, cex.axis = 1.08, cex.lab = 1.20, cex.main = 1.30,
            font.main = 2)
    legend("bottomright", legend = unique(go_plot$component),
           fill = component_cols[unique(go_plot$component)], bty = "n", cex = 1.0)
  } else {
    plot.new(); title("D  Existing GO-BP enrichment")
    text(0.5, 0.5, "No terms at BH FDR < 0.05")
  }
}

pdf_fun <- function(filename) pdf(filename, width = 12, height = 9, useDingbats = FALSE)
png_fun <- function(filename) png(filename, width = 3600, height = 2700, res = 400, type = "cairo")
draw_figure(pdf_fun, file.path(output_dir, "Response_Figure_R2_16_partial_R2_and_GO.pdf"))
draw_figure(png_fun, file.path(output_dir, "Response_Figure_R2_16_partial_R2_and_GO.png"))

if (requireNamespace("svglite", quietly = TRUE)) {
  draw_figure(function(filename) svglite::svglite(filename, width = 12, height = 9),
              file.path(output_dir, "Response_Figure_R2_16_partial_R2_and_GO.svg"))
}
draw_figure(function(filename) tiff(filename, width = 12, height = 9, units = "in",
                                    res = 600, compression = "lzw", type = "cairo"),
            file.path(output_dir, "Response_Figure_R2_16_partial_R2_and_GO.tiff"))

cat("Completed Comment 16 partial-R2 and functional-validation analysis.\n")
cat("Residual df:", residual_df, "\n")
print(summary_top)
print(agreement)
cat("GO terms plotted:", nrow(go_plot), "\n")
