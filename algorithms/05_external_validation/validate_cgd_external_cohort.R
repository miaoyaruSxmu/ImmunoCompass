options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(ggplot2)
  library(limma)
  library(e1071)
})

root <- "C:/temp/ibi_tool"
input_dir <- file.path(root, "revision_supp_data")
output_dir <- file.path(root, "Result", "revision_GSE44969_validation_pm1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_series <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  stopifnot(length(begin) == 1L, length(end) == 1L)
  expression <- read.delim(
    path, skip = begin, nrows = end - begin - 2L, row.names = 1,
    check.names = FALSE, quote = "\"", comment.char = ""
  )
  expression <- as.matrix(expression)
  storage.mode(expression) <- "double"
  header <- lines[seq_len(min(100L, length(lines)))]
  extract_field <- function(prefix) {
    values <- strsplit(header[grep(prefix, header)[1L]], "\t", fixed = TRUE)[[1L]][-1L]
    sub('^"(.*)"$', "\\1", values)
  }
  characteristic_lines <- header[grep("^!Sample_characteristics_ch1", header)]
  extract_characteristic <- function(key) {
    line <- characteristic_lines[
      grep(paste0('^!Sample_characteristics_ch1\\t"', key, ':'),
           characteristic_lines, ignore.case = TRUE)
    ]
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    trimws(sub(paste0("^", key, ":"), "",
               sub('^"(.*)"$', "\\1", values), ignore.case = TRUE))
  }
  metadata <- data.frame(
    sample_id = extract_field("^!Sample_geo_accession"),
    sample_title = extract_field("^!Sample_title"),
    clinical_group = extract_characteristic("group"),
    age_group = extract_characteristic("age"),
    sex = extract_characteristic("sex")
  )
  stopifnot(identical(colnames(expression), metadata$sample_id))
  list(expression = expression, metadata = metadata)
}

map_probes_to_genes <- function(expression, annotation) {
  probe <- trimws(annotation$Probe_Id)
  gene <- toupper(trimws(annotation$Symbol))
  mapped_gene <- gene[match(rownames(expression), probe)]
  keep <- !is.na(mapped_gene) & mapped_gene != ""
  mapped <- avereps(expression[keep, , drop = FALSE], ID = mapped_gene[keep])
  mapped[apply(mapped, 1L, function(x) all(is.finite(x))), , drop = FALSE]
}

build_gene_pairs <- function(expression, pairs) {
  feature_matrix <- matrix(
    -1L, nrow = ncol(expression), ncol = length(pairs),
    dimnames = list(colnames(expression), pairs)
  )
  present <- logical(length(pairs))
  for (i in seq_along(pairs)) {
    genes <- strsplit(pairs[i], "_vs_", fixed = TRUE)[[1L]]
    present[i] <- all(genes %in% rownames(expression))
    if (present[i]) {
      feature_matrix[, i] <- ifelse(
        expression[genes[1L], ] > expression[genes[2L], ], 1L, -1L
      )
    }
  }
  list(x = feature_matrix, present = present)
}

load(file.path(root, "Result", "revision_GSE51404_validation_pm1",
               "updated_svr_with_GSE51405_pm1.rda"))
stopifnot(length(features) == 390L)
series <- read_series(file.path(input_dir, "GSE44969_series_matrix.txt"))
annotation <- read.delim(
  file.path(input_dir, "GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx"),
  skip = 8L, quote = "", check.names = FALSE, fill = TRUE
)
expression <- map_probes_to_genes(series$expression, annotation)
pair_data <- build_gene_pairs(expression, features)

coverage <- data.frame(
  dataset = "GSE44969 independent validation",
  samples = ncol(expression), mapped_genes = nrow(expression),
  gene_pairs_present = sum(pair_data$present),
  gene_pairs_missing = sum(!pair_data$present)
)
write.csv(coverage, file.path(output_dir, "platform_and_feature_coverage.csv"), row.names = FALSE)
if (any(!pair_data$present)) {
  warning("Some model gene pairs are absent; see platform_and_feature_coverage.csv")
}

metadata <- series$metadata
metadata$clinical_group <- factor(
  metadata$clinical_group,
  levels = c("Chronic granulomatous patient", "Carrier", "Healthy control"),
  labels = c("CGD", "Carrier", "Healthy")
)
metadata$age_group <- factor(metadata$age_group, levels = c("Child", "Adult"))
metadata$sex <- factor(metadata$sex, levels = c("Female", "Male"))
stopifnot(!anyNA(metadata$clinical_group), !anyNA(metadata$age_group), !anyNA(metadata$sex))
metadata$IBI <- as.numeric(predict(updated_svr, pair_data$x[, features, drop = FALSE]))
write.csv(metadata, file.path(output_dir, "GSE44969_sample_level_IBI.csv"), row.names = FALSE)

group_summary <- do.call(rbind, lapply(split(metadata, metadata$clinical_group), function(d) {
  data.frame(
    clinical_group = as.character(d$clinical_group[1L]), n = nrow(d),
    mean = mean(d$IBI), sd = sd(d$IBI), median = median(d$IBI),
    q1 = unname(quantile(d$IBI, 0.25)), q3 = unname(quantile(d$IBI, 0.75)),
    minimum = min(d$IBI), maximum = max(d$IBI)
  )
}))
rownames(group_summary) <- NULL
write.csv(group_summary, file.path(output_dir, "GSE44969_IBI_group_summary.csv"), row.names = FALSE)

global_test <- kruskal.test(IBI ~ clinical_group, data = metadata)
comparisons <- list(c("CGD", "Healthy"), c("Carrier", "Healthy"), c("CGD", "Carrier"))
pairwise_tests <- do.call(rbind, lapply(comparisons, function(groups) {
  first <- metadata$IBI[metadata$clinical_group == groups[1L]]
  second <- metadata$IBI[metadata$clinical_group == groups[2L]]
  test <- suppressWarnings(wilcox.test(first, second, exact = FALSE))
  data.frame(
    comparison = paste(groups, collapse = " vs "), n_1 = length(first), n_2 = length(second),
    median_1 = median(first), median_2 = median(second),
    median_difference = median(first) - median(second),
    statistic = unname(test$statistic), p_value = test$p.value
  )
}))
group_tests <- rbind(
  data.frame(
    comparison = "Overall Kruskal-Wallis", n_1 = NA, n_2 = NA,
    median_1 = NA, median_2 = NA, median_difference = NA,
    statistic = unname(global_test$statistic), p_value = global_test$p.value
  ),
  pairwise_tests
)
write.csv(group_tests, file.path(output_dir, "GSE44969_IBI_group_tests.csv"), row.names = FALSE)

palette <- c(CGD = "#315A7D", Carrier = "#D79053", Healthy = "#777777")
distribution_plot <- ggplot(metadata, aes(clinical_group, IBI, colour = clinical_group)) +
  geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = 0.35, linetype = "dashed") +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white",
               colour = "black", linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, size = 1.7, alpha = 0.82) +
  scale_colour_manual(values = palette, guide = "none") +
  labs(
    title = "IBI distribution in the independent GSE44969 cohort",
    subtitle = sprintf("Kruskal-Wallis P = %.3g; n = 24 independent participants",
                       global_test$p.value),
    x = NULL, y = "Immune Balance Index (IBI)"
  ) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 9)
  )
cairo_pdf(file.path(output_dir, "GSE44969_IBI_distribution.pdf"),
          width = 110 / 25.4, height = 90 / 25.4, family = "Arial")
print(distribution_plot)
dev.off()

selected_signatures <- c(
  "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
signature_labels <- c(
  HALLMARK_ALLOGRAFT_REJECTION = "Allograft rejection",
  HALLMARK_COMPLEMENT = "Complement",
  HALLMARK_IL2_STAT5_SIGNALING = "IL2-STAT5 signaling",
  HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6-JAK-STAT3 signaling",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-alpha response",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma response",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF-alpha/NF-kB signaling"
)
model_genes <- unique(toupper(unlist(strsplit(features, "_vs_", fixed = TRUE))))
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- hallmark[hallmark$gs_name %in% selected_signatures, ]
signature_sets <- lapply(selected_signatures, function(signature) {
  setdiff(unique(toupper(hallmark$gene_symbol[hallmark$gs_name == signature])), model_genes)
})
names(signature_sets) <- signature_labels[selected_signatures]
signature_scores <- gsva(
  ssgseaParam(expression, signature_sets, minSize = 10, maxSize = 500, normalize = TRUE),
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)

partial_spearman <- function(x, y, covariates) {
  complete <- complete.cases(x, y, covariates)
  x <- x[complete]
  y <- y[complete]
  covariates <- droplevels(covariates[complete, , drop = FALSE])
  design <- model.matrix(~ ., data = covariates)
  residual_x <- residuals(lm.fit(design, rank(x)))
  residual_y <- residuals(lm.fit(design, rank(y)))
  rho <- cor(residual_x, residual_y)
  covariate_df <- qr(design)$rank - 1L
  degrees_freedom <- length(x) - covariate_df - 2L
  statistic <- rho * sqrt(degrees_freedom / pmax(1 - rho^2, .Machine$double.eps))
  p_value <- 2 * pt(-abs(statistic), df = degrees_freedom)
  fisher_se <- 1 / sqrt(length(x) - covariate_df - 3L)
  fisher_z <- atanh(pmax(pmin(rho, 0.999999), -0.999999))
  c(
    n = length(x), rho = rho,
    ci_low = tanh(fisher_z - qnorm(0.975) * fisher_se),
    ci_high = tanh(fisher_z + qnorm(0.975) * fisher_se),
    p_value = p_value
  )
}

correlation_rows <- list()
for (signature in rownames(signature_scores)) {
  score <- as.numeric(signature_scores[signature, metadata$sample_id])
  overall_test <- suppressWarnings(cor.test(metadata$IBI, score, method = "spearman", exact = FALSE))
  overall_rho <- unname(overall_test$estimate)
  overall_se <- 1 / sqrt(nrow(metadata) - 3L)
  overall_z <- atanh(pmax(pmin(overall_rho, 0.999999), -0.999999))
  correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
    signature = signature, analysis = "Overall Spearman", subgroup = "All",
    n = nrow(metadata), rho = overall_rho,
    ci_low = tanh(overall_z - qnorm(0.975) * overall_se),
    ci_high = tanh(overall_z + qnorm(0.975) * overall_se),
    p_value = overall_test$p.value
  )
  group_adjusted <- partial_spearman(
    metadata$IBI, score, data.frame(clinical_group = metadata$clinical_group)
  )
  correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
    signature = signature, analysis = "Clinical-group adjusted", subgroup = "All",
    n = group_adjusted["n"], rho = group_adjusted["rho"],
    ci_low = group_adjusted["ci_low"], ci_high = group_adjusted["ci_high"],
    p_value = group_adjusted["p_value"]
  )
  fully_adjusted <- partial_spearman(
    metadata$IBI, score,
    data.frame(
      clinical_group = metadata$clinical_group,
      age_group = metadata$age_group,
      sex = metadata$sex
    )
  )
  correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
    signature = signature, analysis = "Group, age and sex adjusted", subgroup = "All",
    n = fully_adjusted["n"], rho = fully_adjusted["rho"],
    ci_low = fully_adjusted["ci_low"], ci_high = fully_adjusted["ci_high"],
    p_value = fully_adjusted["p_value"]
  )
  for (group_name in levels(metadata$clinical_group)) {
    index <- metadata$clinical_group == group_name
    if (sum(index) < 6L) next
    within_test <- suppressWarnings(cor.test(
      metadata$IBI[index], score[index], method = "spearman", exact = FALSE
    ))
    within_rho <- unname(within_test$estimate)
    within_se <- 1 / sqrt(sum(index) - 3L)
    within_z <- atanh(pmax(pmin(within_rho, 0.999999), -0.999999))
    correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
      signature = signature, analysis = "Within clinical group", subgroup = group_name,
      n = sum(index), rho = within_rho,
      ci_low = tanh(within_z - qnorm(0.975) * within_se),
      ci_high = tanh(within_z + qnorm(0.975) * within_se),
      p_value = within_test$p.value
    )
  }
}
correlations <- do.call(rbind, correlation_rows)
correlations$FDR <- ave(
  correlations$p_value,
  interaction(correlations$analysis, correlations$subgroup, drop = TRUE),
  FUN = function(p) p.adjust(p, method = "BH")
)
write.csv(correlations, file.path(output_dir, "GSE44969_IBI_signature_correlations.csv"),
          row.names = FALSE)

signature_coverage <- data.frame(
  signature = names(signature_sets),
  genes_after_model_gene_exclusion = lengths(signature_sets),
  genes_present_in_GSE44969 = vapply(
    signature_sets, function(genes) sum(genes %in% rownames(expression)), integer(1L)
  )
)
write.csv(signature_coverage, file.path(output_dir, "GSE44969_signature_gene_coverage.csv"),
          row.names = FALSE)

plot_correlations <- correlations[
  correlations$analysis %in% c(
    "Overall Spearman", "Clinical-group adjusted", "Group, age and sex adjusted"
  ),
]
plot_correlations$analysis <- factor(
  plot_correlations$analysis,
  levels = c("Overall Spearman", "Clinical-group adjusted", "Group, age and sex adjusted")
)
plot_correlations$signature <- factor(
  plot_correlations$signature, levels = rev(unname(signature_labels[selected_signatures]))
)
plot_correlations$label <- sprintf(
  "%.2f%s", plot_correlations$rho, ifelse(plot_correlations$p_value < 0.05, "*", "")
)
correlation_plot <- ggplot(plot_correlations, aes(analysis, signature, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#D65F5F",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(
    title = "Association of IBI with immune signatures in GSE44969",
    subtitle = "Signature scores exclude all model genes; *P < 0.05",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold", size = 9),
    legend.position = "right"
  )
cairo_pdf(file.path(output_dir, "GSE44969_IBI_immune_signature_correlations.pdf"),
          width = 145 / 25.4, height = 100 / 25.4, family = "Arial")
print(correlation_plot)
dev.off()

provenance <- data.frame(
  item = c(
    "Model", "Feature selection", "Gene-pair encoding", "Independent validation",
    "Primary correlation analysis", "Sensitivity correlation analysis",
    "Within-group correlation minimum n", "Statistical threshold"
  ),
  value = c(
    "SVR updated with GSE51405; GSE44969 not used for training or tuning",
    "Existing 390 gene pairs retained without re-screening",
    "1 when gene A expression > gene B expression; otherwise -1",
    "All 24 GSE44969 participants",
    "Partial Spearman correlation controlling clinical group",
    "Partial Spearman correlation controlling clinical group, age category and sex",
    "6 independent participants",
    "Two-sided P < 0.05"
  )
)
write.csv(provenance, file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("GSE44969 validation completed:", output_dir, "\n")
