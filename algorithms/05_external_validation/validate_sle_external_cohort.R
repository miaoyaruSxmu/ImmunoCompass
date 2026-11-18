options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(e1071)
  library(ggplot2)
  library(limma)
  library(pROC)
})

root <- "C:/temp/ibi_tool"
input_dir <- file.path(root, "revision_supp_data")
model_dir <- file.path(root, "Result", "revision_60pct_training_fullcohort_pm1")
output_dir <- file.path(root, "Result", "revision_GSE65391_external_validation_pm1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qa_dir <- Sys.getenv("IBI_QA_PNG_DIR", unset = "")
skip_formal_figures <- identical(Sys.getenv("IBI_SKIP_FORMAL_FIGURES", unset = "0"), "1")
save_qa_preview <- function(plot, filename, width, height) {
  if (!nzchar(qa_dir)) return(invisible(NULL))
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(qa_dir, filename), plot = plot,
    width = width, height = height, units = "in", dpi = 180
  )
}

series_path <- file.path(input_dir, "GSE65391_series_matrix.txt")
annotation_path <- file.path(input_dir, "GPL10558-50081.txt")
model_path <- file.path(model_dir, "updated_svr_60pct_training_fullcohort_pm1.rda")

read_header <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con))
  lines <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("Series matrix table marker was not found")
    lines <- c(lines, line)
    if (grepl("^!series_matrix_table_begin", line)) break
  }
  lines
}

extract_field <- function(header, prefix) {
  line <- header[grep(prefix, header)[1L]]
  values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
  sub('^"(.*)"$', "\\1", values)
}

extract_characteristic <- function(header, key) {
  lines <- header[grep("^!Sample_characteristics_ch1", header)]
  line <- lines[grep(
    paste0('^!Sample_characteristics_ch1\\t"', key, ':'),
    lines, ignore.case = TRUE
  )]
  if (!length(line)) stop("Missing characteristic: ", key)
  values <- strsplit(line[1L], "\t", fixed = TRUE)[[1L]][-1L]
  trimws(sub(
    paste0("^", key, ":"), "",
    sub('^"(.*)"$', "\\1", values), ignore.case = TRUE
  ))
}

as_numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(ifelse(
    x %in% c("Data Not Available", "Not Applicable", "", "NA"), NA, x
  )))
}

header <- read_header(series_path)
keys <- c(
  "batch", "batch_replicate", "subject", "visit", "set",
  "disease state", "gender", "race", "age", "sledai",
  "disease_activity", "wbc", "neutrophil_count", "lymphocyte_count",
  "monocyte_count", "esr", "c3", "c4"
)
metadata <- data.frame(
  sample_id = extract_field(header, "^!Sample_geo_accession"),
  sample_title = extract_field(header, "^!Sample_title"),
  stringsAsFactors = FALSE
)
for (key in keys) metadata[[gsub(" ", "_", key, fixed = TRUE)]] <-
  extract_characteristic(header, key)
stopifnot(nrow(metadata) == 996L)

numeric_fields <- c(
  "visit", "age", "sledai", "disease_activity", "wbc",
  "neutrophil_count", "lymphocyte_count", "monocyte_count", "esr", "c3", "c4"
)
for (field in numeric_fields) metadata[[field]] <- as_numeric_or_na(metadata[[field]])
metadata$technical_replicate <- metadata$set == "Technical_Replicate"
metadata$clinical_group <- factor(metadata$disease_state, levels = c("Healthy", "SLE"))

load(model_path)
stopifnot(exists("updated_svr"), exists("features"), length(features) == 390L)
load(file.path(root, "X0723.rda"))
stopifnot(exists("X"), all(features %in% colnames(X)))
development_x <- as.matrix(X[, features, drop = FALSE])
feature_modes <- ifelse(colMeans(development_x == 1L) >= 0.5, 1L, -1L)
names(feature_modes) <- features
model_genes <- unique(toupper(unlist(strsplit(features, "_vs_", fixed = TRUE))))

annotation <- read.delim(
  annotation_path, comment.char = "#", quote = "",
  check.names = FALSE, fill = TRUE
)
annotation$ID <- trimws(annotation$ID)
annotation$Symbol <- toupper(trimws(annotation$Symbol))
annotation <- annotation[
  annotation$Symbol %in% model_genes & annotation$ID != "",
  c("ID", "Symbol"), drop = FALSE
]
target_probes <- unique(annotation$ID)

table_begin <- grep("^!series_matrix_table_begin", header)
raw_expression <- fread(
  series_path, skip = table_begin, sep = "\t", header = TRUE,
  quote = "\"", data.table = FALSE, showProgress = TRUE
)
names(raw_expression)[1L] <- "probe_id"
raw_expression$probe_id <- trimws(raw_expression$probe_id)
raw_expression <- raw_expression[raw_expression$probe_id %in% target_probes, , drop = FALSE]
probe_ids <- raw_expression$probe_id
raw_expression$probe_id <- NULL
expression <- as.matrix(raw_expression)
storage.mode(expression) <- "double"
rownames(expression) <- probe_ids
rm(raw_expression)
gc()

stopifnot(identical(colnames(expression), metadata$sample_id))
mapped_gene <- annotation$Symbol[match(rownames(expression), annotation$ID)]
keep <- !is.na(mapped_gene) & mapped_gene != ""
gene_expression <- avereps(expression[keep, , drop = FALSE], ID = mapped_gene[keep])
rm(expression)
gc()

pair_matrix <- matrix(
  rep(feature_modes[features], each = ncol(gene_expression)),
  nrow = ncol(gene_expression), ncol = length(features),
  dimnames = list(colnames(gene_expression), features)
)
pair_present <- logical(length(features))
for (i in seq_along(features)) {
  genes <- strsplit(features[i], "_vs_", fixed = TRUE)[[1L]]
  pair_present[i] <- all(genes %in% rownames(gene_expression))
  if (pair_present[i]) {
    pair_matrix[, i] <- ifelse(
      gene_expression[genes[1L], ] > gene_expression[genes[2L], ], 1L, -1L
    )
  }
}
missing_pairs <- features[!pair_present]
missing_genes <- unique(unlist(strsplit(missing_pairs, "_vs_", fixed = TRUE)))
missing_genes <- missing_genes[!missing_genes %in% rownames(gene_expression)]
write.csv(
  data.frame(
    gene_pair = missing_pairs,
    development_set_modal_encoding = feature_modes[missing_pairs],
    reason = "At least one gene probe absent from processed GSE65391 matrix"
  ),
  file.path(output_dir, "GSE65391_imputed_gene_pairs.csv"), row.names = FALSE
)

metadata$IBI <- as.numeric(predict(updated_svr, pair_matrix[, features, drop = FALSE]))
pair_matrix_all_minus <- pair_matrix
pair_matrix_all_plus <- pair_matrix
pair_matrix_all_minus[, !pair_present] <- -1L
pair_matrix_all_plus[, !pair_present] <- 1L
metadata$IBI_missing_all_minus <- as.numeric(predict(
  updated_svr, pair_matrix_all_minus[, features, drop = FALSE]
))
metadata$IBI_missing_all_plus <- as.numeric(predict(
  updated_svr, pair_matrix_all_plus[, features, drop = FALSE]
))
load(file.path(root, "best_svr0731.rda"))
metadata$IBI_original_model <- as.numeric(predict(
  best_svr, pair_matrix[, features, drop = FALSE]
))
nontechnical <- metadata[!metadata$technical_replicate, , drop = FALSE]
nontechnical <- nontechnical[order(
  nontechnical$subject,
  ifelse(is.na(nontechnical$visit), Inf, nontechnical$visit),
  nontechnical$sample_id
), , drop = FALSE]
primary <- nontechnical[!duplicated(nontechnical$subject), , drop = FALSE]
primary$analysis_set <- "Earliest non-technical sample per participant"

stopifnot(
  nrow(nontechnical) == 972L,
  nrow(primary) == 204L,
  sum(primary$clinical_group == "SLE") == 158L,
  sum(primary$clinical_group == "Healthy") == 46L
)

write.csv(metadata, file.path(output_dir, "GSE65391_all_sample_predictions.csv"), row.names = FALSE)
write.csv(
  nontechnical,
  file.path(output_dir, "GSE65391_all_nontechnical_predictions.csv"),
  row.names = FALSE
)
write.csv(
  primary,
  file.path(output_dir, "GSE65391_primary_independent_subject_predictions.csv"),
  row.names = FALSE
)

response <- as.integer(primary$clinical_group == "SLE")
roc_object <- roc(
  response = response, predictor = primary$IBI,
  levels = c(0, 1), direction = "<", quiet = TRUE
)
auc_value <- as.numeric(auc(roc_object))
auc_ci <- as.numeric(ci.auc(roc_object, method = "delong"))
wilcox_result <- wilcox.test(
  IBI ~ clinical_group, data = primary, exact = FALSE, conf.int = TRUE
)
auc_results <- data.frame(
  dataset = "GSE65391",
  evaluation_scope = "Independent earliest-visit subject-level validation",
  comparison = "SLE vs Healthy",
  n_case = sum(response == 1L),
  n_healthy = sum(response == 0L),
  auc = auc_value,
  ci_low = auc_ci[1L],
  ci_high = auc_ci[3L],
  score_direction = "Higher IBI indicates SLE",
  independent_from_model_training = TRUE
)
write.csv(auc_results, file.path(output_dir, "GSE65391_primary_AUC.csv"), row.names = FALSE)

calculate_directional_auc <- function(score, label) {
  object <- roc(
    response = label, predictor = score,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  interval <- as.numeric(ci.auc(object, method = "delong"))
  data.frame(
    auc_expected_direction = as.numeric(auc(object)),
    ci_low = interval[1L], ci_high = interval[3L],
    auc_if_direction_reversed = 1 - as.numeric(auc(object))
  )
}
sensitivity_auc <- do.call(rbind, lapply(
  c("IBI", "IBI_missing_all_minus", "IBI_missing_all_plus", "IBI_original_model"),
  function(score_column) {
    cbind(
      data.frame(score_definition = score_column),
      calculate_directional_auc(primary[[score_column]], response)
    )
  }
))
write.csv(
  sensitivity_auc,
  file.path(output_dir, "GSE65391_missing_feature_and_model_sensitivity_AUC.csv"),
  row.names = FALSE
)

summary_table <- do.call(rbind, lapply(
  split(primary, primary$clinical_group),
  function(d) data.frame(
    clinical_group = as.character(d$clinical_group[1L]),
    n = nrow(d), mean = mean(d$IBI), sd = sd(d$IBI), median = median(d$IBI),
    q1 = unname(quantile(d$IBI, 0.25)), q3 = unname(quantile(d$IBI, 0.75)),
    minimum = min(d$IBI), maximum = max(d$IBI)
  )
))
rownames(summary_table) <- NULL
write.csv(summary_table, file.path(output_dir, "GSE65391_primary_group_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(
    comparison = "SLE vs Healthy", statistic = unname(wilcox_result$statistic),
    p_value = wilcox_result$p.value,
    median_SLE = median(primary$IBI[primary$clinical_group == "SLE"]),
    median_Healthy = median(primary$IBI[primary$clinical_group == "Healthy"])
  ),
  file.path(output_dir, "GSE65391_primary_group_test.csv"), row.names = FALSE
)

palette <- c("Healthy" = "#777777", "SLE" = "#C44E52")
distribution_plot <- ggplot(primary, aes(clinical_group, IBI, colour = clinical_group)) +
  geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = 0.35, linetype = "dashed") +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, fill = "white",
    colour = "black", linewidth = 0.4
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 20261018L),
    size = 1.25, alpha = 0.68
  ) +
  scale_colour_manual(values = palette, guide = "none") +
  labs(
    title = "IBI distribution in the independent GSE65391 cohort",
    subtitle = sprintf(
      "Earliest non-technical sample per participant; P = %.3g",
      wilcox_result$p.value
    ),
    x = NULL, y = "Immune Balance Index (IBI)"
  ) +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 9)
  )
if (!skip_formal_figures) {
  cairo_pdf(
    file.path(output_dir, "GSE65391_primary_IBI_distribution.pdf"),
    width = 100 / 25.4, height = 88 / 25.4, family = "Arial"
  )
  print(distribution_plot)
  dev.off()
}
save_qa_preview(distribution_plot, "GSE65391_primary_IBI_distribution.png", 100 / 25.4, 88 / 25.4)

roc_coordinates <- data.frame(
  dataset = "GSE65391", comparison = "SLE vs Healthy",
  specificity = roc_object$specificities,
  sensitivity = roc_object$sensitivities,
  false_positive_rate = 1 - roc_object$specificities,
  true_positive_rate = roc_object$sensitivities
)
write.csv(
  roc_coordinates,
  file.path(output_dir, "GSE65391_primary_ROC_coordinates.csv"), row.names = FALSE
)
roc_plot <- ggplot(roc_coordinates, aes(false_positive_rate, true_positive_rate)) +
  geom_abline(
    intercept = 0, slope = 1, linewidth = 0.45,
    linetype = "dashed", colour = "#777777"
  ) +
  geom_path(linewidth = 0.9, colour = "#C44E52", lineend = "round") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(
    title = "SLE classification in GSE65391",
    subtitle = sprintf(
      "AUC %.3f (95%% CI %.3f-%.3f); cases/HC %d/%d",
      auc_value, auc_ci[1L], auc_ci[3L], sum(response == 1L), sum(response == 0L)
    ),
    x = "1 - Specificity", y = "Sensitivity"
  ) +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 9)
  )
if (!skip_formal_figures) {
  cairo_pdf(
    file.path(output_dir, "GSE65391_primary_ROC_curve.pdf"),
    width = 100 / 25.4, height = 92 / 25.4, family = "Arial"
  )
  print(roc_plot)
  dev.off()
}
save_qa_preview(roc_plot, "GSE65391_primary_ROC_curve.png", 100 / 25.4, 92 / 25.4)

clinical_variables <- c(
  sledai = "SLEDAI",
  disease_activity = "Disease activity category",
  c3 = "Complement C3",
  c4 = "Complement C4",
  esr = "Erythrocyte sedimentation rate",
  wbc = "White blood cell count",
  neutrophil_count = "Neutrophil count",
  lymphocyte_count = "Lymphocyte count",
  monocyte_count = "Monocyte count"
)
sle_primary <- primary[primary$clinical_group == "SLE", , drop = FALSE]
correlation_rows <- lapply(names(clinical_variables), function(variable) {
  keep <- is.finite(sle_primary$IBI) & is.finite(sle_primary[[variable]])
  test <- suppressWarnings(cor.test(
    sle_primary$IBI[keep], sle_primary[[variable]][keep],
    method = "spearman", exact = FALSE
  ))
  data.frame(
    variable = variable,
    indicator = unname(clinical_variables[variable]),
    n = sum(keep), rho = unname(test$estimate), p_value = test$p.value
  )
})
clinical_correlations <- do.call(rbind, correlation_rows)
clinical_correlations$FDR_BH <- p.adjust(clinical_correlations$p_value, method = "BH")
write.csv(
  clinical_correlations,
  file.path(output_dir, "GSE65391_primary_IBI_clinical_correlations.csv"),
  row.names = FALSE
)

correlation_plot_data <- clinical_correlations
correlation_plot_data$indicator <- factor(
  correlation_plot_data$indicator,
  levels = rev(correlation_plot_data$indicator)
)
correlation_plot_data$label <- sprintf(
  "rho = %.2f%s",
  correlation_plot_data$rho,
  ifelse(correlation_plot_data$p_value < 0.05, "*", "")
)
correlation_plot <- ggplot(
  correlation_plot_data,
  aes(x = rho, y = indicator, colour = rho)
) +
  geom_vline(xintercept = 0, colour = "#A8A8A8", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rho, yend = indicator), linewidth = 0.65) +
  geom_point(size = 2.4) +
  geom_text(
    aes(label = label), colour = "black", size = 2.4,
    hjust = ifelse(correlation_plot_data$rho >= 0, -0.08, 1.08)
  ) +
  scale_colour_gradient2(
    low = "#4C78A8", mid = "#777777", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1), guide = "none"
  ) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    title = "Association of IBI with clinical immune indicators",
    subtitle = "SLE participants at earliest available visit; *P < 0.05",
    x = "Spearman correlation coefficient", y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 9),
    plot.margin = margin(5, 18, 5, 5, unit = "pt")
  )
if (!skip_formal_figures) {
  cairo_pdf(
    file.path(output_dir, "GSE65391_primary_IBI_clinical_correlations.pdf"),
    width = 135 / 25.4, height = 95 / 25.4, family = "Arial"
  )
  print(correlation_plot)
  dev.off()
}
save_qa_preview(
  correlation_plot, "GSE65391_primary_IBI_clinical_correlations.png",
  135 / 25.4, 95 / 25.4
)

write.csv(
  data.frame(
    item = c(
      "Dataset role", "Model", "Feature policy", "Gene-pair encoding",
      "Primary analysis unit", "Technical replicate policy",
      "Longitudinal sample policy", "Missing-feature policy",
      "AUC orientation", "Statistical threshold",
      "Output figures"
    ),
    value = c(
      "Independent external validation; no GSE65391 samples used for training or tuning",
      "SVR augmented with GSE51405 and 60% subsets of GSE51404 and GSE44969",
      "Existing 390 gene pairs retained without re-screening",
      "1 when gene A expression > gene B expression; otherwise -1",
      "Earliest non-technical sample per unique participant (158 SLE and 46 healthy)",
      "Samples labelled Technical_Replicate excluded from primary and longitudinal analyses",
      "All 972 non-technical samples retained in a separate prediction table",
      paste0(
        sum(!pair_present), " of 390 pairs unavailable in the processed matrix; ",
        "imputed with each pair's modal -1/1 encoding in the original development matrix"
      ),
      "Higher IBI indicates SLE",
      "Two-sided P < 0.05; BH-adjusted values also reported for clinical correlations",
      "PDF only"
    )
  ),
  file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    dataset = "GSE65391", total_arrays = nrow(metadata),
    technical_replicate_arrays = sum(metadata$technical_replicate),
    nontechnical_arrays = nrow(nontechnical), unique_participants = nrow(primary),
    SLE_participants = sum(primary$clinical_group == "SLE"),
    healthy_participants = sum(primary$clinical_group == "Healthy"),
    mapped_model_genes = sum(model_genes %in% rownames(gene_expression)),
    total_model_genes = length(model_genes),
    gene_pairs_present = sum(pair_present), gene_pairs_imputed = sum(!pair_present),
    total_gene_pairs = length(features),
    missing_genes = paste(missing_genes, collapse = ";")
  ),
  file.path(output_dir, "sample_and_feature_coverage.csv"), row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("GSE65391 independent validation completed:", output_dir, "\n")
