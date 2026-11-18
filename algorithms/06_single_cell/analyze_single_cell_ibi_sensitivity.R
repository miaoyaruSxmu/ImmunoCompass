options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(e1071)
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

root <- "C:/temp/ibi_tool"
output_dir <- file.path(
  root, "Result", "revision_781_single_cell_sensitivity_top390_20261028"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

feature_file <- file.path(
  root, "Result", "revision_781_fixed_top390_final",
  "fixed_top390_all_samples_781_development.rda"
)
model_file <- file.path(
  root, "Result", "revision_781_final_svr_external_validation_20261028_corrected",
  "final_781_fixed_top390_svr_model_and_predictions.rda"
)
single_cell_file <- file.path(
  root, "Aging", "scRNA-seqProcessedLabelledObject.rds"
)

base_seed <- 20261025L
prediction_chunk_size <- 10000L
minimum_cells_per_donor_celltype <- 50L
downsample_cells_per_donor_celltype <- 10L
downsample_repeats <- 10L
downsample_fractions <- c(0.25, 0.50, 0.75)

development_cohorts <- c(
  "GTEx", "GSE232381", "GSE297198", "GSE264125", "GSE222408",
  "GSE175839", "GSE181228", "GSE199911", "GSE140713",
  "GSE59216_batch1", "GSE51405"
)

# Figure contract
# Core conclusion: quantify whether sparse single-cell detection and sequencing
# depth materially alter IBI, and whether donor-cell-type pseudobulk aggregation
# recovers concordant immune-state estimates.
# Evidence chain: feature coverage -> alternative missing-value encoding ->
# pseudobulk concordance -> controlled UMI downsampling.
# Archetype: quantitative grid. Backend: R only. Export: editable PDF only,
# 183-mm width, with plotted source data written as CSV.

message("Loading fixed Top390 feature matrix and final 781-sample reference model")
feature_environment <- new.env(parent = emptyenv())
load(feature_file, envir = feature_environment)
model_environment <- new.env(parent = emptyenv())
load(model_file, envir = model_environment)

features <- feature_environment$features_top390
all_feature_matrix <- as.matrix(feature_environment$top390_feature_matrix)
sample_manifest <- feature_environment$sample_manifest
reference_model <- model_environment$final_svr_model
final_prepared <- model_environment$final_prepared

stopifnot(
  length(features) == 390L,
  identical(colnames(all_feature_matrix), features),
  inherits(reference_model, "svm"),
  identical(final_prepared$kept_features, features),
  all(all_feature_matrix[is.finite(all_feature_matrix)] %in% c(0, 1))
)

development_index <- which(sample_manifest$cohort %in% development_cohorts)
stopifnot(
  length(development_index) == 781L,
  sum(sample_manifest$cohort[development_index] == "GSE59216_batch1") == 91L
)
feature_prevalence <- colMeans(all_feature_matrix[development_index, , drop = FALSE])
feature_modes <- ifelse(feature_prevalence >= 0.5, 1, 0)
names(feature_modes) <- features

apply_preprocessing <- function(x) {
  x <- as.matrix(x)
  for (j in seq_len(ncol(x))) {
    x[!is.finite(x[, j]), j] <- final_prepared$modes[j]
  }
  x <- x[, final_prepared$kept_features, drop = FALSE]
  sweep(
    sweep(x, 2L, final_prepared$centers, "-"),
    2L, final_prepared$scales, "/"
  )
}

predict_ibi <- function(x) {
  as.numeric(predict(reference_model, apply_preprocessing(x)))
}

pair_split <- strsplit(features, "_vs_", fixed = TRUE)
pair_table <- data.frame(
  rank = seq_along(features),
  gene_pair = features,
  gene_a = vapply(pair_split, `[`, character(1), 1L),
  gene_b = vapply(pair_split, `[`, character(1), 2L),
  development_prevalence = as.numeric(feature_prevalence[features]),
  development_modal_encoding = as.numeric(feature_modes[features]),
  stringsAsFactors = FALSE
)

message("Loading single-cell object")
single_cell <- readRDS(single_cell_file)
stopifnot(inherits(single_cell, "Seurat"), "RNA" %in% Assays(single_cell))

metadata <- single_cell[[]]
metadata$cell_id <- rownames(metadata)
metadata$sampleName <- as.character(metadata$sampleName)
metadata$primary_type <- as.character(metadata$primary_type)
metadata$primary_type[is.na(metadata$primary_type) | metadata$primary_type == ""] <-
  "Unclassified"
stopifnot(
  !anyNA(metadata$sampleName),
  all(c("nCount_RNA", "nFeature_RNA", "percent.mt", "age") %in% colnames(metadata))
)

all_feature_genes <- unique(c(pair_table$gene_a, pair_table$gene_b))
present_feature_genes <- intersect(all_feature_genes, rownames(single_cell))
counts <- LayerData(
  single_cell[["RNA"]], layer = "counts", features = present_feature_genes
)
stopifnot(identical(colnames(counts), metadata$cell_id))
rm(single_cell)
gc()

pair_table$gene_a_present <- pair_table$gene_a %in% rownames(counts)
pair_table$gene_b_present <- pair_table$gene_b %in% rownames(counts)
pair_table$pair_present <- pair_table$gene_a_present & pair_table$gene_b_present
pair_table$gene_a_index <- match(pair_table$gene_a, rownames(counts))
pair_table$gene_b_index <- match(pair_table$gene_b, rownames(counts))
pair_table$single_cell_default_encoding <- ifelse(
  pair_table$pair_present, NA_real_, 0
)
pair_table$missing_reason <- ifelse(
  pair_table$pair_present, "Both genes present in the single-cell matrix",
  "At least one gene absent from the single-cell matrix"
)
write.csv(
  pair_table,
  file.path(output_dir, "single_cell_feature_gene_coverage.csv"),
  row.names = FALSE
)

available_pair_n <- sum(pair_table$pair_present)
stopifnot(available_pair_n > 0L)

build_pair_matrices <- function(dense_expression, include_modal = TRUE) {
  n_profiles <- ncol(dense_expression)
  current <- matrix(
    0, nrow = n_profiles, ncol = nrow(pair_table),
    dimnames = list(colnames(dense_expression), pair_table$gene_pair)
  )
  modal <- if (include_modal) current else NULL
  both_detected_n <- integer(n_profiles)
  one_detected_n <- integer(n_profiles)
  both_zero_n <- integer(n_profiles)
  positive_tie_n <- integer(n_profiles)

  for (pair_index in seq_len(nrow(pair_table))) {
    if (!pair_table$pair_present[pair_index]) {
      if (include_modal) {
        modal[, pair_index] <- pair_table$development_modal_encoding[pair_index]
      }
      next
    }

    gene_a <- dense_expression[pair_table$gene_a_index[pair_index], ]
    gene_b <- dense_expression[pair_table$gene_b_index[pair_index], ]
    both_detected <- gene_a > 0 & gene_b > 0
    one_detected <- xor(gene_a > 0, gene_b > 0)
    both_zero <- gene_a == 0 & gene_b == 0
    positive_tie <- gene_a == gene_b & gene_a > 0
    encoded <- as.numeric(gene_a > gene_b)

    current[, pair_index] <- encoded
    if (include_modal) {
      modal_encoded <- encoded
      if (pair_table$development_modal_encoding[pair_index] == 1) {
        modal_encoded[both_zero] <- 1
      }
      modal[, pair_index] <- modal_encoded
    }

    both_detected_n <- both_detected_n + both_detected
    one_detected_n <- one_detected_n + one_detected
    both_zero_n <- both_zero_n + both_zero
    positive_tie_n <- positive_tie_n + positive_tie
  }

  list(
    current = current,
    modal = modal,
    both_detected_fraction = both_detected_n / available_pair_n,
    one_detected_fraction = one_detected_n / available_pair_n,
    both_zero_fraction = both_zero_n / available_pair_n,
    positive_tie_fraction = positive_tie_n / available_pair_n,
    at_least_one_detected_fraction =
      (both_detected_n + one_detected_n) / available_pair_n
  )
}

predict_sparse_counts <- function(
    sparse_counts, include_modal = TRUE, include_coverage = TRUE,
    progress_label = "profiles") {
  profile_n <- ncol(sparse_counts)
  current_prediction <- numeric(profile_n)
  modal_prediction <- if (include_modal) numeric(profile_n) else NULL
  coverage_names <- c(
    "both_detected_fraction", "one_detected_fraction",
    "both_zero_fraction", "positive_tie_fraction",
    "at_least_one_detected_fraction"
  )
  coverage <- if (include_coverage) {
    matrix(
      NA_real_, nrow = profile_n, ncol = length(coverage_names),
      dimnames = list(colnames(sparse_counts), coverage_names)
    )
  } else {
    NULL
  }

  chunk_starts <- seq.int(1L, profile_n, by = prediction_chunk_size)
  for (chunk_number in seq_along(chunk_starts)) {
    start <- chunk_starts[chunk_number]
    end <- min(profile_n, start + prediction_chunk_size - 1L)
    profile_index <- start:end
    dense_chunk <- as.matrix(sparse_counts[, profile_index, drop = FALSE])
    built <- build_pair_matrices(dense_chunk, include_modal = include_modal)
    current_prediction[profile_index] <- predict_ibi(
      built$current[, features, drop = FALSE]
    )
    if (include_modal) {
      modal_prediction[profile_index] <- predict_ibi(
        built$modal[, features, drop = FALSE]
      )
    }
    if (include_coverage) {
      for (coverage_name in coverage_names) {
        coverage[profile_index, coverage_name] <- built[[coverage_name]]
      }
    }
    if (chunk_number %% 10L == 0L || chunk_number == length(chunk_starts)) {
      message(
        "Completed ", progress_label, " chunk ", chunk_number,
        " of ", length(chunk_starts)
      )
    }
  }

  list(
    current = current_prediction,
    modal = modal_prediction,
    coverage = coverage
  )
}

metric_row <- function(reference, estimate, scope, n_override = NULL) {
  complete <- is.finite(reference) & is.finite(estimate)
  reference <- reference[complete]
  estimate <- estimate[complete]
  n_value <- if (is.null(n_override)) length(reference) else n_override
  data.frame(
    scope = scope,
    n = n_value,
    spearman_rho = suppressWarnings(cor(reference, estimate, method = "spearman")),
    pearson_r = suppressWarnings(cor(reference, estimate, method = "pearson")),
    mae = mean(abs(reference - estimate)),
    median_absolute_difference = median(abs(reference - estimate)),
    mean_difference = mean(estimate - reference),
    stringsAsFactors = FALSE
  )
}

cell_results_file <- file.path(
  output_dir, "cell_level_qc_and_ibi_sensitivity.rds"
)
if (file.exists(cell_results_file)) {
  message("Loading cached cell-level IBI sensitivity results")
  cell_results <- readRDS(cell_results_file)
  stopifnot(
    nrow(cell_results) == nrow(metadata),
    identical(cell_results$cell_id, metadata$cell_id)
  )
} else {
  message("Predicting cell-level IBI under default and modal missing-value rules")
  cell_prediction <- predict_sparse_counts(
    counts, include_modal = TRUE, include_coverage = TRUE,
    progress_label = "single-cell"
  )

  cell_results <- data.frame(
    cell_id = metadata$cell_id,
    sampleName = metadata$sampleName,
    primary_type = metadata$primary_type,
    secondary_type = as.character(metadata$secondary_type),
    age = metadata$age,
    nCount_RNA = metadata$nCount_RNA,
    nFeature_RNA = metadata$nFeature_RNA,
    percent.mt = metadata$percent.mt,
    IBI_default_zero = cell_prediction$current,
    IBI_training_mode = cell_prediction$modal,
    IBI_difference_mode_minus_zero =
      cell_prediction$modal - cell_prediction$current,
    cell_prediction$coverage,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cell_results$group_id <- paste(
    cell_results$sampleName, cell_results$primary_type, sep = "|||"
  )
  saveRDS(cell_results, cell_results_file, compress = "xz")
}

group_summary <- cell_results %>%
  group_by(group_id, sampleName, primary_type) %>%
  summarise(
    n_cells = n(),
    age = first(age),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    IBI_default_zero_median = median(IBI_default_zero, na.rm = TRUE),
    IBI_training_mode_median = median(IBI_training_mode, na.rm = TRUE),
    IBI_difference_median = median(
      IBI_difference_mode_minus_zero, na.rm = TRUE
    ),
    both_detected_fraction_median = median(
      both_detected_fraction, na.rm = TRUE
    ),
    one_detected_fraction_median = median(
      one_detected_fraction, na.rm = TRUE
    ),
    both_zero_fraction_median = median(both_zero_fraction, na.rm = TRUE),
    positive_tie_fraction_median = median(
      positive_tie_fraction, na.rm = TRUE
    ),
    at_least_one_detected_fraction_median = median(
      at_least_one_detected_fraction, na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(n_cells >= minimum_cells_per_donor_celltype)

write.csv(
  group_summary,
  file.path(output_dir, "donor_celltype_summary.csv"), row.names = FALSE
)

celltype_summary <- group_summary %>%
  group_by(primary_type) %>%
  summarise(
    donor_n = n_distinct(sampleName),
    donor_celltype_group_n = n(),
    median_cell_n = median(n_cells),
    median_at_least_one_detected_fraction = median(
      at_least_one_detected_fraction_median
    ),
    median_both_detected_fraction = median(both_detected_fraction_median),
    median_both_zero_fraction = median(both_zero_fraction_median),
    median_default_ibi = median(IBI_default_zero_median),
    median_modal_ibi = median(IBI_training_mode_median),
    median_encoding_difference = median(IBI_difference_median),
    .groups = "drop"
  )
write.csv(
  celltype_summary,
  file.path(output_dir, "celltype_sensitivity_summary.csv"), row.names = FALSE
)

encoding_metrics <- bind_rows(
  metric_row(
    cell_results$IBI_default_zero, cell_results$IBI_training_mode,
    "All cells", n_override = nrow(cell_results)
  ),
  metric_row(
    group_summary$IBI_default_zero_median,
    group_summary$IBI_training_mode_median,
    "Donor-by-primary-cell-type medians",
    n_override = nrow(group_summary)
  )
)

encoding_metrics_by_celltype <- bind_rows(lapply(
  split(group_summary, group_summary$primary_type),
  function(current_data) {
    output <- metric_row(
      current_data$IBI_default_zero_median,
      current_data$IBI_training_mode_median,
      unique(current_data$primary_type),
      n_override = n_distinct(current_data$sampleName)
    )
    output$primary_type <- unique(current_data$primary_type)
    output
  }
))
write.csv(
  encoding_metrics,
  file.path(output_dir, "missing_encoding_overall_metrics.csv"),
  row.names = FALSE
)
write.csv(
  encoding_metrics_by_celltype,
  file.path(output_dir, "missing_encoding_metrics_by_celltype.csv"),
  row.names = FALSE
)

z_score <- function(x) as.numeric(scale(x))
mixed_data <- group_summary %>%
  mutate(
    primary_type = factor(primary_type),
    sampleName = factor(sampleName),
    z_ibi_default = z_score(IBI_default_zero_median),
    z_abs_encoding_difference = z_score(abs(IBI_difference_median)),
    z_log_ncount = z_score(log10(median_nCount_RNA + 1)),
    z_log_nfeature = z_score(log10(median_nFeature_RNA + 1)),
    z_both_zero = z_score(both_zero_fraction_median)
  )

fit_specs <- list(
  IBI_by_library_size = c("z_ibi_default", "z_log_ncount"),
  IBI_by_detected_genes = c("z_ibi_default", "z_log_nfeature"),
  IBI_by_both_zero_fraction = c("z_ibi_default", "z_both_zero"),
  Encoding_difference_by_both_zero_fraction =
    c("z_abs_encoding_difference", "z_both_zero")
)

mixed_model_rows <- lapply(names(fit_specs), function(model_name) {
  outcome <- fit_specs[[model_name]][1]
  predictor <- fit_specs[[model_name]][2]
  formula <- as.formula(
    paste(outcome, "~", predictor, "+ primary_type + (1 | sampleName)")
  )
  fit <- lmer(
    formula, data = mixed_data, REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      check.conv.singular = "ignore"
    )
  )
  coefficient_table <- coef(summary(fit))
  estimate <- coefficient_table[predictor, "Estimate"]
  standard_error <- coefficient_table[predictor, "Std. Error"]
  z_value <- estimate / standard_error
  data.frame(
    model = model_name,
    outcome = outcome,
    predictor = predictor,
    standardized_beta = estimate,
    standard_error = standard_error,
    ci_lower = estimate - 1.96 * standard_error,
    ci_upper = estimate + 1.96 * standard_error,
    wald_p = 2 * pnorm(abs(z_value), lower.tail = FALSE),
    donor_n = n_distinct(mixed_data$sampleName),
    donor_celltype_group_n = nrow(mixed_data),
    singular_fit = isSingular(fit, tol = 1e-4),
    stringsAsFactors = FALSE
  )
})
mixed_model_results <- bind_rows(mixed_model_rows)
mixed_model_results$bh_fdr <- p.adjust(mixed_model_results$wald_p, method = "BH")
write.csv(
  mixed_model_results,
  file.path(output_dir, "technical_covariate_mixed_models.csv"),
  row.names = FALSE
)

message("Constructing donor-by-cell-type pseudobulk profiles")
all_group_factor <- factor(cell_results$group_id)
aggregation_matrix <- sparseMatrix(
  i = seq_along(all_group_factor),
  j = as.integer(all_group_factor),
  x = 1,
  dims = c(length(all_group_factor), nlevels(all_group_factor)),
  dimnames = list(cell_results$cell_id, levels(all_group_factor))
)
pseudobulk_all <- counts %*% aggregation_matrix
valid_group_ids <- group_summary$group_id
pseudobulk <- pseudobulk_all[, valid_group_ids, drop = FALSE]
rm(aggregation_matrix, pseudobulk_all)
gc()

pseudobulk_built <- build_pair_matrices(
  as.matrix(pseudobulk), include_modal = TRUE
)
pseudobulk_default <- predict_ibi(
  pseudobulk_built$current[, features, drop = FALSE]
)
pseudobulk_modal <- predict_ibi(
  pseudobulk_built$modal[, features, drop = FALSE]
)

pseudobulk_results <- group_summary
pseudobulk_results$IBI_pseudobulk_default_zero <- pseudobulk_default
pseudobulk_results$IBI_pseudobulk_training_mode <- pseudobulk_modal
pseudobulk_results$pseudobulk_both_detected_fraction <-
  pseudobulk_built$both_detected_fraction
pseudobulk_results$pseudobulk_at_least_one_detected_fraction <-
  pseudobulk_built$at_least_one_detected_fraction
pseudobulk_results$pseudobulk_both_zero_fraction <-
  pseudobulk_built$both_zero_fraction
write.csv(
  pseudobulk_results,
  file.path(output_dir, "donor_celltype_pseudobulk_predictions.csv"),
  row.names = FALSE
)

pseudobulk_metrics <- bind_rows(
  metric_row(
    pseudobulk_results$IBI_default_zero_median,
    pseudobulk_results$IBI_pseudobulk_default_zero,
    "Default-zero encoding"
  ),
  metric_row(
    pseudobulk_results$IBI_training_mode_median,
    pseudobulk_results$IBI_pseudobulk_training_mode,
    "Training-mode encoding"
  )
)

pseudobulk_metrics_by_celltype <- bind_rows(lapply(
  split(pseudobulk_results, pseudobulk_results$primary_type),
  function(current_data) {
    default_metrics <- metric_row(
      current_data$IBI_default_zero_median,
      current_data$IBI_pseudobulk_default_zero,
      "Default-zero encoding",
      n_override = n_distinct(current_data$sampleName)
    )
    modal_metrics <- metric_row(
      current_data$IBI_training_mode_median,
      current_data$IBI_pseudobulk_training_mode,
      "Training-mode encoding",
      n_override = n_distinct(current_data$sampleName)
    )
    output <- bind_rows(default_metrics, modal_metrics)
    output$primary_type <- unique(current_data$primary_type)
    output
  }
))
write.csv(
  pseudobulk_metrics,
  file.path(output_dir, "pseudobulk_concordance_overall.csv"),
  row.names = FALSE
)
write.csv(
  pseudobulk_metrics_by_celltype,
  file.path(output_dir, "pseudobulk_concordance_by_celltype.csv"),
  row.names = FALSE
)

message("Running stratified UMI downsampling sensitivity analysis")
set.seed(base_seed)
cell_indices_by_group <- split(seq_len(nrow(cell_results)), cell_results$group_id)
downsample_indices <- sort(unlist(lapply(
  cell_indices_by_group,
  function(index) sample(index, min(length(index), downsample_cells_per_donor_celltype))
), use.names = FALSE))
downsample_counts <- counts[, downsample_indices, drop = FALSE]
downsample_metadata <- cell_results[
  downsample_indices,
  c("cell_id", "sampleName", "primary_type", "group_id", "IBI_default_zero")
]
write.csv(
  downsample_metadata,
  file.path(output_dir, "downsampling_cell_manifest.csv"),
  row.names = FALSE
)

baseline_group_median <- tapply(
  downsample_metadata$IBI_default_zero,
  downsample_metadata$group_id,
  median
)

downsampling_rows <- list()
row_counter <- 1L
for (retained_fraction in downsample_fractions) {
  for (iteration in seq_len(downsample_repeats)) {
    set.seed(base_seed + as.integer(retained_fraction * 1000) + iteration)
    thinned_counts <- downsample_counts
    thinned_counts@x <- as.numeric(rbinom(
      length(thinned_counts@x),
      size = as.integer(round(thinned_counts@x)),
      prob = retained_fraction
    ))
    thinned_counts <- drop0(thinned_counts)
    thinned_prediction <- predict_sparse_counts(
      thinned_counts, include_modal = FALSE, include_coverage = FALSE,
      progress_label = paste0(
        "downsampling ", retained_fraction, " repeat ", iteration
      )
    )$current

    cell_metrics <- metric_row(
      downsample_metadata$IBI_default_zero,
      thinned_prediction,
      "Stratified cells"
    )
    cell_metrics$retained_fraction <- retained_fraction
    cell_metrics$iteration <- iteration
    cell_metrics$analysis_level <- "Cell"
    downsampling_rows[[row_counter]] <- cell_metrics
    row_counter <- row_counter + 1L

    thinned_group_median <- tapply(
      thinned_prediction, downsample_metadata$group_id, median
    )
    common_groups <- intersect(
      names(baseline_group_median), names(thinned_group_median)
    )
    group_metrics <- metric_row(
      baseline_group_median[common_groups],
      thinned_group_median[common_groups],
      "Donor-by-primary-cell-type medians"
    )
    group_metrics$retained_fraction <- retained_fraction
    group_metrics$iteration <- iteration
    group_metrics$analysis_level <- "Donor-cell-type"
    downsampling_rows[[row_counter]] <- group_metrics
    row_counter <- row_counter + 1L

    message(
      "Completed retained fraction ", retained_fraction,
      ", repeat ", iteration, " of ", downsample_repeats
    )
  }
}

downsampling_metrics <- bind_rows(downsampling_rows)
write.csv(
  downsampling_metrics,
  file.path(output_dir, "downsampling_repeat_metrics.csv"),
  row.names = FALSE
)

downsampling_summary <- downsampling_metrics %>%
  group_by(analysis_level, retained_fraction) %>%
  summarise(
    repeat_n = n(),
    spearman_rho_median = median(spearman_rho),
    spearman_rho_q025 = quantile(spearman_rho, 0.025),
    spearman_rho_q975 = quantile(spearman_rho, 0.975),
    mae_median = median(mae),
    mae_q025 = quantile(mae, 0.025),
    mae_q975 = quantile(mae, 0.975),
    median_absolute_difference_median = median(median_absolute_difference),
    .groups = "drop"
  )
write.csv(
  downsampling_summary,
  file.path(output_dir, "downsampling_metric_summary.csv"),
  row.names = FALSE
)

palette_contract <- c(
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  signal_blue = "#3F729B",
  signal_teal = "#3A8F82",
  accent_red = "#B44D4D",
  accent_orange = "#D08A35"
)

theme_nature <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 1.0, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank()
    )
}

save_pdf <- function(plot, filename, width_mm = 183, height_mm = 130) {
  cairo_pdf(
    file.path(output_dir, filename),
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    family = "Arial"
  )
  print(plot)
  dev.off()
}

save_qa_png <- function(plot, filename, width_mm, height_mm) {
  qa_dir <- Sys.getenv("SINGLE_CELL_QA_DIR", unset = "")
  if (!nzchar(qa_dir)) return(invisible(NULL))
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
  png(
    file.path(qa_dir, filename), width = width_mm / 25.4,
    height = height_mm / 25.4, units = "in", res = 180,
    type = "cairo", family = "Arial"
  )
  print(plot)
  dev.off()
  invisible(NULL)
}

celltype_order <- celltype_summary$primary_type[
  order(celltype_summary$median_at_least_one_detected_fraction)
]
group_summary$primary_type_plot <- factor(
  group_summary$primary_type, levels = celltype_order
)

coverage_long <- group_summary %>%
  select(
    sampleName, primary_type_plot,
    at_least_one_detected_fraction_median,
    both_detected_fraction_median
  ) %>%
  pivot_longer(
    cols = c(
      at_least_one_detected_fraction_median,
      both_detected_fraction_median
    ),
    names_to = "coverage_metric", values_to = "fraction"
  ) %>%
  mutate(
    coverage_metric = recode(
      coverage_metric,
      at_least_one_detected_fraction_median = "At least one gene detected",
      both_detected_fraction_median = "Both genes detected"
    )
  )
write.csv(
  coverage_long,
  file.path(output_dir, "source_data_coverage_panel.csv"), row.names = FALSE
)

p_coverage <- ggplot(
  coverage_long,
  aes(primary_type_plot, fraction, fill = coverage_metric)
) +
  geom_boxplot(
    width = 0.66, outlier.shape = NA, linewidth = 0.3,
    position = position_dodge(width = 0.72)
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    "At least one gene detected" = unname(palette_contract["signal_blue"]),
    "Both genes detected" = unname(palette_contract["signal_teal"])
  )) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Gene-pair detectability varies across cell types",
    x = NULL, y = "Fraction of 388 evaluable pairs", fill = NULL
  ) +
  theme_nature() +
  theme(legend.position = "top")

encoding_group_metrics <- metric_row(
  group_summary$IBI_default_zero_median,
  group_summary$IBI_training_mode_median,
  "Donor-by-primary-cell-type medians"
)

p_encoding_scatter <- ggplot(
  group_summary,
  aes(IBI_default_zero_median, IBI_training_mode_median)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 0.9, alpha = 0.55, colour = palette_contract["signal_teal"]
  ) +
  coord_equal() +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
    label = sprintf(
      "Spearman rho = %.3f\nMedian |difference| = %.3f\nn = %d donor-cell-type groups",
      encoding_group_metrics$spearman_rho,
      encoding_group_metrics$median_absolute_difference,
      nrow(group_summary)
    ),
    size = 2.2
  ) +
  labs(
    title = "Alternative missing-value encoding",
    x = "IBI: unavailable or double-zero pair = 0",
    y = "IBI: development-set modal encoding"
  ) +
  theme_nature()

p_encoding_delta <- ggplot(
  group_summary,
  aes(primary_type_plot, IBI_difference_median)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  geom_boxplot(
    width = 0.62, outlier.shape = NA, linewidth = 0.3,
    fill = "#E8D6C0", colour = palette_contract["accent_orange"]
  ) +
  coord_flip() +
  labs(
    title = "Encoding-related shift by cell type",
    x = NULL, y = "Median IBI difference (modal - zero)"
  ) +
  theme_nature()

p_zero_dependence <- ggplot(
  group_summary,
  aes(both_zero_fraction_median, abs(IBI_difference_median))
) +
  geom_point(
    size = 0.9, alpha = 0.5, colour = palette_contract["accent_red"]
  ) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.55,
    colour = palette_contract["neutral_dark"], fill = "#D9D9D9"
  ) +
  labs(
    title = "Encoding sensitivity increases with double-zero pairs",
    x = "Median fraction of double-zero pairs",
    y = "Absolute median IBI difference"
  ) +
  theme_nature()

coverage_figure <-
  (p_coverage | p_encoding_scatter) /
  (p_encoding_delta | p_zero_dependence) +
  plot_annotation(tag_levels = "a")
save_pdf(
  coverage_figure,
  "single_cell_coverage_and_missing_encoding_sensitivity.pdf",
  width_mm = 183, height_mm = 150
)
save_qa_png(
  coverage_figure,
  "single_cell_coverage_and_missing_encoding_sensitivity.png",
  width_mm = 183, height_mm = 150
)

p_pb_default <- ggplot(
  pseudobulk_results,
  aes(IBI_default_zero_median, IBI_pseudobulk_default_zero)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 1.0, alpha = 0.55, colour = palette_contract["signal_blue"]
  ) +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
    label = sprintf(
      "Spearman rho = %.3f\nMAE = %.3f\nn = %d donor-cell-type groups",
      pseudobulk_metrics$spearman_rho[1], pseudobulk_metrics$mae[1],
      nrow(pseudobulk_results)
    ),
    size = 2.2
  ) +
  labs(
    title = "Default-zero encoding",
    x = "Median cell-level IBI", y = "Pseudobulk IBI"
  ) +
  theme_nature()

p_pb_modal <- ggplot(
  pseudobulk_results,
  aes(IBI_training_mode_median, IBI_pseudobulk_training_mode)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 1.0, alpha = 0.55, colour = palette_contract["signal_teal"]
  ) +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
    label = sprintf(
      "Spearman rho = %.3f\nMAE = %.3f\nn = %d donor-cell-type groups",
      pseudobulk_metrics$spearman_rho[2], pseudobulk_metrics$mae[2],
      nrow(pseudobulk_results)
    ),
    size = 2.2
  ) +
  labs(
    title = "Development-mode encoding",
    x = "Median cell-level IBI", y = "Pseudobulk IBI"
  ) +
  theme_nature()

p_pb_celltype <- pseudobulk_metrics_by_celltype %>%
  mutate(
    primary_type = factor(primary_type, levels = rev(celltype_order)),
    scope = factor(
      scope,
      levels = c("Default-zero encoding", "Training-mode encoding")
    )
  ) %>%
  ggplot(aes(spearman_rho, primary_type, colour = scope, shape = scope)) +
  geom_vline(
    xintercept = 0, linewidth = 0.3,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(
    "Default-zero encoding" = unname(palette_contract["signal_blue"]),
    "Training-mode encoding" = unname(palette_contract["signal_teal"])
  )) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    title = "Pseudobulk concordance within cell types",
    x = "Spearman rho across donors", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_nature() +
  theme(legend.position = "top")

p_pb_coverage <- ggplot(
  pseudobulk_results,
  aes(
    at_least_one_detected_fraction_median,
    pseudobulk_at_least_one_detected_fraction
  )
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 1.0, alpha = 0.5, colour = palette_contract["accent_orange"]
  ) +
  labs(
    title = "Aggregation restores feature detectability",
    x = "Median cell-level detectable-pair fraction",
    y = "Pseudobulk detectable-pair fraction"
  ) +
  theme_nature()

pseudobulk_figure <-
  (p_pb_default | p_pb_modal) /
  (p_pb_celltype | p_pb_coverage) +
  plot_annotation(
    title = "Pseudobulk aggregation provides an orthogonal IBI sensitivity check",
    tag_levels = "a"
  )
save_pdf(
  pseudobulk_figure,
  "single_cell_pseudobulk_concordance.pdf",
  width_mm = 183, height_mm = 145
)
save_qa_png(
  pseudobulk_figure,
  "single_cell_pseudobulk_concordance.png",
  width_mm = 183, height_mm = 145
)

downsampling_plot_data <- downsampling_metrics %>%
  filter(analysis_level == "Donor-cell-type") %>%
  mutate(
    retained_label = factor(
      paste0(as.integer(retained_fraction * 100), "%"),
      levels = c("25%", "50%", "75%")
    )
  )
write.csv(
  downsampling_plot_data,
  file.path(output_dir, "source_data_downsampling_panels.csv"),
  row.names = FALSE
)

p_down_rho <- ggplot(
  downsampling_plot_data,
  aes(retained_label, spearman_rho)
) +
  geom_hline(
    yintercept = 1, linetype = "dotted", linewidth = 0.35,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_boxplot(
    width = 0.56, outlier.shape = NA, linewidth = 0.35,
    fill = "#C9D9E6", colour = palette_contract["signal_blue"]
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0), size = 1.5,
    colour = palette_contract["signal_blue"]
  ) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(
    title = "IBI ranking under UMI downsampling",
    subtitle = paste0(
      downsample_repeats,
      " repeats; donor-cell-type medians from a stratified cell subset"
    ),
    x = "Retained UMI fraction", y = "Spearman rho vs. unthinned IBI"
  ) +
  theme_nature()

p_down_mae <- ggplot(
  downsampling_plot_data,
  aes(retained_label, mae)
) +
  geom_hline(
    yintercept = 0, linetype = "dotted", linewidth = 0.35,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_boxplot(
    width = 0.56, outlier.shape = NA, linewidth = 0.35,
    fill = "#E8D6C0", colour = palette_contract["accent_orange"]
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0), size = 1.5,
    colour = palette_contract["accent_orange"]
  ) +
  labs(
    title = "Absolute IBI error under UMI downsampling",
    subtitle = paste0(
      downsample_repeats,
      " repeats; donor-cell-type medians from a stratified cell subset"
    ),
    x = "Retained UMI fraction", y = "MAE vs. unthinned IBI"
  ) +
  theme_nature()

downsampling_figure <- (p_down_rho | p_down_mae) +
  plot_annotation(tag_levels = "a")
save_pdf(
  downsampling_figure,
  "single_cell_umi_downsampling_sensitivity.pdf",
  width_mm = 183, height_mm = 88
)

analysis_provenance <- data.frame(
  item = c(
    "Single-cell object", "Cell count", "Donor count",
    "Primary cell-type count", "Feature source", "Reference model",
    "Feature encoding", "Globally available feature genes",
    "Globally available gene pairs", "Default missing rule",
    "Alternative missing rule", "Pseudobulk unit",
    "Minimum cells per pseudobulk unit", "Downsampling subset rule",
    "Downsampling fractions", "Downsampling repeats", "Random seed",
    "Inference boundary"
  ),
  value = c(
    "Healthy-donor PBMC aging Seurat object",
    nrow(cell_results),
    n_distinct(cell_results$sampleName),
    n_distinct(cell_results$primary_type),
    "Previously generated all-data exploratory Top390 gene-pair panel",
    "Final 781-sample radial epsilon-SVR reference fit",
    "1 when gene A expression is greater than gene B expression; otherwise 0",
    paste0(sum(all_feature_genes %in% rownames(counts)), "/", length(all_feature_genes)),
    paste0(sum(pair_table$pair_present), "/", nrow(pair_table)),
    "Gene pair encoded as 0 when globally unavailable or both genes have zero counts",
    "Globally unavailable and double-zero pairs replaced by the modal encoding in the 781-sample development set",
    "Donor by primary cell type",
    minimum_cells_per_donor_celltype,
    paste0(
      "Up to ", downsample_cells_per_donor_celltype,
      " cells sampled per donor-by-primary-cell-type group"
    ),
    paste(downsample_fractions, collapse = "; "),
    downsample_repeats,
    base_seed,
    "Cells are technical observations; donor or donor-by-cell-type summaries are used for inferential interpretation"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  analysis_provenance,
  file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE
)

legend_text <- c(
  "Figure: single_cell_coverage_and_missing_encoding_sensitivity.pdf",
  "Boxes show the median and interquartile range across donors; whiskers extend to 1.5 times the interquartile range. Each donor-by-primary-cell-type group contained at least 50 cells. Cell-level values are not treated as independent biological replicates.",
  "",
  "Figure: single_cell_pseudobulk_concordance.pdf",
  "Pseudobulk profiles were constructed by summing raw UMI counts within each donor-by-primary-cell-type group containing at least 50 cells. Because library-size normalization and log1p transformation are monotonic within a profile, they do not change within-profile gene ordering. Spearman correlations and MAE are descriptive concordance measures across donor-by-cell-type groups.",
  "",
  "Figure: single_cell_umi_downsampling_sensitivity.pdf",
  paste0(
    "Up to ", downsample_cells_per_donor_celltype,
    " cells were sampled per donor-by-primary-cell-type group. Raw UMI counts were independently thinned to 25%, 50%, or 75% of their original depth in ",
    downsample_repeats,
    " repeats. Boxes show the median and interquartile range across repeats; whiskers extend to 1.5 times the interquartile range. Metrics compare donor-by-cell-type median IBI values with the corresponding unthinned values."
  )
)
writeLines(
  legend_text,
  file.path(output_dir, "figure_legend_draft.txt"), useBytes = TRUE
)

save(
  pair_table, group_summary, celltype_summary, encoding_metrics,
  mixed_model_results, pseudobulk_results, pseudobulk_metrics,
  pseudobulk_metrics_by_celltype, downsampling_metrics,
  downsampling_summary,
  file = file.path(output_dir, "single_cell_sensitivity_analysis_objects.rda")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "R_sessionInfo.txt"), useBytes = TRUE
)

message("Single-cell sensitivity analysis complete: ", output_dir)
