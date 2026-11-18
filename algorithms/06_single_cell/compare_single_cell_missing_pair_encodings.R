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

arguments <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(arguments) >= 1L) {
  arguments[[1L]]
} else {
  Sys.getenv("IBI_TOOL_ROOT", unset = "C:/temp/ibi_tool")
}
output_root <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  Sys.getenv("IBI_TOOL_OUTPUT_ROOT", unset = data_root)
}

data_root <- gsub("\\\\", "/", path.expand(data_root))
output_root <- gsub("\\\\", "/", path.expand(output_root))
output_dir <- file.path(
  output_root, "Result", "revision_711_single_cell_modal_primary_20261113"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_file <- file.path(
  data_root, "Result", "revision_711_top390_feature_selection_LOCO_20261109",
  "full711_selection_and_model.rda"
)
single_cell_file <- file.path(
  data_root, "Aging", "scRNA-seqProcessedLabelledObject.rds"
)

stopifnot(file.exists(model_file), file.exists(single_cell_file))

base_seed <- 20261025L
prediction_chunk_size <- 10000L
minimum_cells_per_donor_celltype <- 50L
downsample_cells_per_donor_celltype <- 10L
downsample_repeats <- 10L
downsample_fractions <- c(0.25, 0.50, 0.75)

# Figure contract
# Core conclusion: IBI rankings from the fixed 711-sample model remain
# interpretable under sparse single-cell measurement when development-set
# modal encoding is the primary rule, while absolute values retain measurable
# sensitivity to dropout, sequencing depth, and the missing-pair rule.
# Evidence chain: encoding agreement -> technical covariates -> controlled UMI
# thinning -> donor-by-cell-type pseudobulk concordance.
# Archetype: quantitative grid. Backend: R only. Final size: 183 mm.

message("Loading the fixed 711-sample Top390 model")
model_environment <- new.env(parent = emptyenv())
load(model_file, envir = model_environment)

selection <- model_environment$selection
reference_model <- model_environment$full_svr_model
development_feature_matrix <- as.matrix(
  model_environment$full_feature_matrix
)
features <- as.character(selection$top390$gene_pair)
feature_modes <- as.numeric(selection$training_modes)
names(feature_modes) <- features

stopifnot(
  length(features) == 390L,
  nrow(development_feature_matrix) == 711L,
  identical(colnames(development_feature_matrix), features),
  inherits(reference_model, "svm"),
  all(development_feature_matrix %in% c(0, 1)),
  all(feature_modes %in% c(0, 1)),
  identical(
    unname(feature_modes),
    unname(ifelse(colMeans(development_feature_matrix) >= 0.5, 1, 0))
  )
)

predict_ibi <- function(feature_matrix) {
  feature_matrix <- as.matrix(feature_matrix)
  stopifnot(identical(colnames(feature_matrix), features))
  as.numeric(predict(reference_model, feature_matrix))
}

pair_table <- data.frame(
  rank = seq_along(features),
  gene_pair = features,
  gene_a = as.character(selection$top390$geneA),
  gene_b = as.character(selection$top390$geneB),
  development_prevalence = colMeans(development_feature_matrix),
  development_modal_encoding = feature_modes,
  stringsAsFactors = FALSE
)

message("Loading the single-cell object")
single_cell <- readRDS(single_cell_file)
stopifnot(inherits(single_cell, "Seurat"), "RNA" %in% Assays(single_cell))

metadata <- single_cell[[]]
metadata$cell_id <- rownames(metadata)
metadata$sampleName <- as.character(metadata$sampleName)
metadata$primary_type <- as.character(metadata$primary_type)
metadata$primary_type[
  is.na(metadata$primary_type) | metadata$primary_type == ""
] <- "Unclassified"
metadata$secondary_type <- if ("secondary_type" %in% colnames(metadata)) {
  as.character(metadata$secondary_type)
} else {
  NA_character_
}

stopifnot(
  !anyNA(metadata$sampleName),
  all(c("nCount_RNA", "nFeature_RNA", "percent.mt") %in% colnames(metadata))
)

all_feature_genes <- unique(c(pair_table$gene_a, pair_table$gene_b))
present_feature_genes <- intersect(all_feature_genes, rownames(single_cell))
counts <- LayerData(
  single_cell[["RNA"]], layer = "counts", features = present_feature_genes
)
stopifnot(identical(colnames(counts), metadata$cell_id))
rm(single_cell)
invisible(gc())

pair_table$gene_a_present <- pair_table$gene_a %in% rownames(counts)
pair_table$gene_b_present <- pair_table$gene_b %in% rownames(counts)
pair_table$pair_present <-
  pair_table$gene_a_present & pair_table$gene_b_present
pair_table$gene_a_index <- match(pair_table$gene_a, rownames(counts))
pair_table$gene_b_index <- match(pair_table$gene_b, rownames(counts))
pair_table$primary_encoding_if_unavailable <-
  pair_table$development_modal_encoding
pair_table$sensitivity_encoding_if_unavailable <- 0
pair_table$missing_reason <- ifelse(
  pair_table$pair_present,
  "Both genes present in the single-cell matrix",
  "At least one gene absent from the single-cell matrix"
)
write.csv(
  pair_table,
  file.path(output_dir, "single_cell_feature_gene_coverage.csv"),
  row.names = FALSE
)

available_pair_n <- sum(pair_table$pair_present)
stopifnot(available_pair_n > 0L)

build_pair_matrices <- function(dense_expression) {
  profile_n <- ncol(dense_expression)
  zero_matrix <- matrix(
    0,
    nrow = profile_n,
    ncol = nrow(pair_table),
    dimnames = list(colnames(dense_expression), pair_table$gene_pair)
  )
  modal_matrix <- matrix(
    rep(feature_modes, each = profile_n),
    nrow = profile_n,
    ncol = nrow(pair_table),
    dimnames = list(colnames(dense_expression), pair_table$gene_pair)
  )
  both_detected_n <- integer(profile_n)
  one_detected_n <- integer(profile_n)
  both_zero_n <- integer(profile_n)
  positive_tie_n <- integer(profile_n)

  for (pair_index in which(pair_table$pair_present)) {
    gene_a <- dense_expression[pair_table$gene_a_index[pair_index], ]
    gene_b <- dense_expression[pair_table$gene_b_index[pair_index], ]
    both_detected <- gene_a > 0 & gene_b > 0
    one_detected <- xor(gene_a > 0, gene_b > 0)
    both_zero <- gene_a == 0 & gene_b == 0
    positive_tie <- gene_a == gene_b & gene_a > 0
    encoded <- as.numeric(gene_a > gene_b)

    zero_matrix[, pair_index] <- encoded
    modal_matrix[, pair_index] <- encoded
    modal_matrix[both_zero, pair_index] <-
      feature_modes[pair_index]

    both_detected_n <- both_detected_n + both_detected
    one_detected_n <- one_detected_n + one_detected
    both_zero_n <- both_zero_n + both_zero
    positive_tie_n <- positive_tie_n + positive_tie
  }

  list(
    modal = modal_matrix,
    zero = zero_matrix,
    both_detected_fraction = both_detected_n / available_pair_n,
    one_detected_fraction = one_detected_n / available_pair_n,
    both_zero_fraction = both_zero_n / available_pair_n,
    positive_tie_fraction = positive_tie_n / available_pair_n,
    at_least_one_detected_fraction =
      (both_detected_n + one_detected_n) / available_pair_n
  )
}

predict_sparse_counts <- function(
    sparse_counts, include_coverage = TRUE, progress_label = "profiles") {
  profile_n <- ncol(sparse_counts)
  modal_prediction <- numeric(profile_n)
  zero_prediction <- numeric(profile_n)
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
    built <- build_pair_matrices(dense_chunk)
    modal_prediction[profile_index] <- predict_ibi(
      built$modal[, features, drop = FALSE]
    )
    zero_prediction[profile_index] <- predict_ibi(
      built$zero[, features, drop = FALSE]
    )
    if (include_coverage) {
      for (coverage_name in coverage_names) {
        coverage[profile_index, coverage_name] <- built[[coverage_name]]
      }
    }
    if (chunk_number %% 10L == 0L ||
        chunk_number == length(chunk_starts)) {
      message(
        "Completed ", progress_label, " chunk ", chunk_number,
        " of ", length(chunk_starts)
      )
    }
  }

  list(
    modal = modal_prediction,
    zero = zero_prediction,
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
    spearman_rho = suppressWarnings(
      cor(reference, estimate, method = "spearman")
    ),
    pearson_r = suppressWarnings(
      cor(reference, estimate, method = "pearson")
    ),
    mae = mean(abs(reference - estimate)),
    median_absolute_difference = median(abs(reference - estimate)),
    mean_difference = mean(estimate - reference),
    stringsAsFactors = FALSE
  )
}

cell_results_file <- file.path(
  output_dir, "cell_level_qc_and_ibi_modal_primary.rds"
)
if (file.exists(cell_results_file)) {
  message("Loading cached cell-level predictions")
  cell_results <- readRDS(cell_results_file)
  stopifnot(
    nrow(cell_results) == nrow(metadata),
    identical(cell_results$cell_id, metadata$cell_id),
    all(c("IBI_modal_primary", "IBI_zero_sensitivity") %in%
      colnames(cell_results))
  )
} else {
  message("Predicting cell-level IBI under modal-primary and zero rules")
  cell_prediction <- predict_sparse_counts(
    counts, include_coverage = TRUE, progress_label = "single-cell"
  )
  cell_results <- data.frame(
    cell_id = metadata$cell_id,
    sampleName = metadata$sampleName,
    primary_type = metadata$primary_type,
    secondary_type = metadata$secondary_type,
    nCount_RNA = metadata$nCount_RNA,
    nFeature_RNA = metadata$nFeature_RNA,
    percent.mt = metadata$percent.mt,
    IBI_modal_primary = cell_prediction$modal,
    IBI_zero_sensitivity = cell_prediction$zero,
    IBI_difference_zero_minus_modal =
      cell_prediction$zero - cell_prediction$modal,
    cell_prediction$coverage,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cell_results$group_id <- paste(
    cell_results$sampleName, cell_results$primary_type, sep = "|||"
  )
  saveRDS(cell_results, cell_results_file, compress = "gzip")
}

group_summary <- cell_results %>%
  group_by(group_id, sampleName, primary_type) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    IBI_modal_primary_median = median(IBI_modal_primary, na.rm = TRUE),
    IBI_zero_sensitivity_median = median(
      IBI_zero_sensitivity, na.rm = TRUE
    ),
    IBI_difference_median = median(
      IBI_difference_zero_minus_modal, na.rm = TRUE
    ),
    both_detected_fraction_median = median(
      both_detected_fraction, na.rm = TRUE
    ),
    one_detected_fraction_median = median(
      one_detected_fraction, na.rm = TRUE
    ),
    both_zero_fraction_median = median(
      both_zero_fraction, na.rm = TRUE
    ),
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
  file.path(output_dir, "donor_celltype_summary.csv"),
  row.names = FALSE
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
    median_both_detected_fraction = median(
      both_detected_fraction_median
    ),
    median_both_zero_fraction = median(both_zero_fraction_median),
    median_modal_primary_ibi = median(IBI_modal_primary_median),
    median_zero_sensitivity_ibi = median(IBI_zero_sensitivity_median),
    median_encoding_difference = median(IBI_difference_median),
    .groups = "drop"
  )
write.csv(
  celltype_summary,
  file.path(output_dir, "celltype_sensitivity_summary.csv"),
  row.names = FALSE
)

encoding_metrics <- bind_rows(
  metric_row(
    cell_results$IBI_modal_primary,
    cell_results$IBI_zero_sensitivity,
    "All cells: modal primary versus zero sensitivity",
    n_override = nrow(cell_results)
  ),
  metric_row(
    group_summary$IBI_modal_primary_median,
    group_summary$IBI_zero_sensitivity_median,
    "Donor-by-primary-cell-type medians: modal primary versus zero sensitivity",
    n_override = nrow(group_summary)
  )
)

encoding_metrics_by_celltype <- bind_rows(lapply(
  split(group_summary, group_summary$primary_type),
  function(current_data) {
    output <- metric_row(
      current_data$IBI_modal_primary_median,
      current_data$IBI_zero_sensitivity_median,
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
    z_ibi_modal_primary = z_score(IBI_modal_primary_median),
    z_abs_encoding_difference = z_score(abs(IBI_difference_median)),
    z_log_ncount = z_score(log10(median_nCount_RNA + 1)),
    z_log_nfeature = z_score(log10(median_nFeature_RNA + 1)),
    z_both_zero = z_score(both_zero_fraction_median)
  )

fit_specs <- list(
  IBI_by_library_size = list(
    outcome = "z_ibi_modal_primary",
    predictor = "z_log_ncount",
    label = "Library size"
  ),
  IBI_by_detected_genes = list(
    outcome = "z_ibi_modal_primary",
    predictor = "z_log_nfeature",
    label = "Detected-gene count"
  ),
  IBI_by_double_zero_fraction = list(
    outcome = "z_ibi_modal_primary",
    predictor = "z_both_zero",
    label = "Double-zero fraction"
  ),
  Encoding_difference_by_double_zero_fraction = list(
    outcome = "z_abs_encoding_difference",
    predictor = "z_both_zero",
    label = "Double-zero fraction vs. encoding difference"
  )
)

mixed_model_rows <- lapply(names(fit_specs), function(model_name) {
  specification <- fit_specs[[model_name]]
  formula <- as.formula(paste(
    specification$outcome, "~", specification$predictor,
    "+ primary_type + (1 | sampleName)"
  ))
  fit <- lmer(
    formula,
    data = mixed_data,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      check.conv.singular = "ignore"
    )
  )
  coefficient_table <- coef(summary(fit))
  estimate <- coefficient_table[specification$predictor, "Estimate"]
  standard_error <- coefficient_table[
    specification$predictor, "Std. Error"
  ]
  z_value <- estimate / standard_error
  data.frame(
    model = model_name,
    label = specification$label,
    outcome = specification$outcome,
    predictor = specification$predictor,
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
mixed_model_results$bh_fdr <- p.adjust(
  mixed_model_results$wald_p, method = "BH"
)
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
invisible(gc())

pseudobulk_built <- build_pair_matrices(as.matrix(pseudobulk))
pseudobulk_modal <- predict_ibi(
  pseudobulk_built$modal[, features, drop = FALSE]
)
pseudobulk_zero <- predict_ibi(
  pseudobulk_built$zero[, features, drop = FALSE]
)

pseudobulk_results <- group_summary
pseudobulk_results$IBI_pseudobulk_modal_primary <- pseudobulk_modal
pseudobulk_results$IBI_pseudobulk_zero_sensitivity <- pseudobulk_zero
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
    pseudobulk_results$IBI_modal_primary_median,
    pseudobulk_results$IBI_pseudobulk_modal_primary,
    "711-development-set modal encoding"
  ),
  metric_row(
    pseudobulk_results$IBI_zero_sensitivity_median,
    pseudobulk_results$IBI_pseudobulk_zero_sensitivity,
    "Zero encoding sensitivity"
  )
)

pseudobulk_metrics_by_celltype <- bind_rows(lapply(
  split(pseudobulk_results, pseudobulk_results$primary_type),
  function(current_data) {
    modal_metrics <- metric_row(
      current_data$IBI_modal_primary_median,
      current_data$IBI_pseudobulk_modal_primary,
      "711-development-set modal encoding",
      n_override = n_distinct(current_data$sampleName)
    )
    zero_metrics <- metric_row(
      current_data$IBI_zero_sensitivity_median,
      current_data$IBI_pseudobulk_zero_sensitivity,
      "Zero encoding sensitivity",
      n_override = n_distinct(current_data$sampleName)
    )
    output <- bind_rows(modal_metrics, zero_metrics)
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

message("Running stratified UMI-thinning sensitivity analysis")
set.seed(base_seed)
valid_cell_index <- which(cell_results$group_id %in% valid_group_ids)
cell_indices_by_group <- split(
  valid_cell_index,
  cell_results$group_id[valid_cell_index]
)
downsample_indices <- sort(unlist(lapply(
  cell_indices_by_group,
  function(index) {
    sample(index, min(length(index), downsample_cells_per_donor_celltype))
  }
), use.names = FALSE))
downsample_counts <- counts[, downsample_indices, drop = FALSE]
downsample_metadata <- cell_results[
  downsample_indices,
  c(
    "cell_id", "sampleName", "primary_type", "group_id",
    "IBI_modal_primary", "IBI_zero_sensitivity"
  )
]
write.csv(
  downsample_metadata,
  file.path(output_dir, "downsampling_cell_manifest.csv"),
  row.names = FALSE
)

baseline_group_medians <- list(
  modal = tapply(
    downsample_metadata$IBI_modal_primary,
    downsample_metadata$group_id,
    median
  ),
  zero = tapply(
    downsample_metadata$IBI_zero_sensitivity,
    downsample_metadata$group_id,
    median
  )
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
      thinned_counts,
      include_coverage = FALSE,
      progress_label = paste0(
        "UMI thinning ", retained_fraction, " repeat ", iteration
      )
    )

    for (encoding_rule in c("modal", "zero")) {
      baseline_cell <- if (encoding_rule == "modal") {
        downsample_metadata$IBI_modal_primary
      } else {
        downsample_metadata$IBI_zero_sensitivity
      }
      thinned_cell <- thinned_prediction[[encoding_rule]]
      rule_label <- if (encoding_rule == "modal") {
        "711-development-set modal encoding"
      } else {
        "Zero encoding sensitivity"
      }

      cell_metrics <- metric_row(
        baseline_cell, thinned_cell, "Stratified cells"
      )
      cell_metrics$retained_fraction <- retained_fraction
      cell_metrics$iteration <- iteration
      cell_metrics$analysis_level <- "Cell"
      cell_metrics$encoding_rule <- rule_label
      downsampling_rows[[row_counter]] <- cell_metrics
      row_counter <- row_counter + 1L

      thinned_group_median <- tapply(
        thinned_cell, downsample_metadata$group_id, median
      )
      common_groups <- intersect(
        names(baseline_group_medians[[encoding_rule]]),
        names(thinned_group_median)
      )
      group_metrics <- metric_row(
        baseline_group_medians[[encoding_rule]][common_groups],
        thinned_group_median[common_groups],
        "Donor-by-primary-cell-type medians"
      )
      group_metrics$retained_fraction <- retained_fraction
      group_metrics$iteration <- iteration
      group_metrics$analysis_level <- "Donor-cell-type"
      group_metrics$encoding_rule <- rule_label
      downsampling_rows[[row_counter]] <- group_metrics
      row_counter <- row_counter + 1L
    }

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
  group_by(analysis_level, encoding_rule, retained_fraction) %>%
  summarise(
    repeat_n = n(),
    spearman_rho_median = median(spearman_rho),
    spearman_rho_q025 = quantile(spearman_rho, 0.025),
    spearman_rho_q975 = quantile(spearman_rho, 0.975),
    mae_median = median(mae),
    mae_q025 = quantile(mae, 0.025),
    mae_q975 = quantile(mae, 0.975),
    median_absolute_difference_median = median(
      median_absolute_difference
    ),
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
  modal = "#2F7F73",
  zero = "#477DA8",
  accent_red = "#B44D4D",
  accent_orange = "#C98535"
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
      plot.title = element_text(size = base_size + 0.6, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

save_pub_r <- function(
    plot, filename, width_mm = 183, height_mm = 145, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(
    file.path(output_dir, paste0(filename, ".svg")),
    width = width_in, height = height_in
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    file.path(output_dir, paste0(filename, ".pdf")),
    width = width_in, height = height_in,
    family = "Arial"
  )
  print(plot)
  dev.off()

  temporary_tiff <- tempfile(fileext = ".tiff")
  ragg::agg_tiff(
    temporary_tiff,
    width = width_in, height = height_in,
    units = "in", res = dpi, background = "white"
  )
  print(plot)
  dev.off()
  stopifnot(file.copy(
    temporary_tiff,
    file.path(output_dir, paste0(filename, ".tiff")),
    overwrite = TRUE
  ))
  unlink(temporary_tiff)

  temporary_png <- tempfile(fileext = ".png")
  ragg::agg_png(
    temporary_png,
    width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  dev.off()
  stopifnot(file.copy(
    temporary_png,
    file.path(output_dir, paste0(filename, ".png")),
    overwrite = TRUE
  ))
  unlink(temporary_png)
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
    names_to = "coverage_metric",
    values_to = "fraction"
  ) %>%
  mutate(
    coverage_metric = recode(
      coverage_metric,
      at_least_one_detected_fraction_median =
        "At least one gene detected",
      both_detected_fraction_median = "Both genes detected"
    )
  )
write.csv(
  coverage_long,
  file.path(output_dir, "source_data_coverage_panel.csv"),
  row.names = FALSE
)

encoding_group_metrics <- encoding_metrics[2L, ]
encoding_difference_model <- mixed_model_results[
  mixed_model_results$model ==
    "Encoding_difference_by_double_zero_fraction",
]

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
    "At least one gene detected" = unname(palette_contract["zero"]),
    "Both genes detected" = unname(palette_contract["modal"])
  )) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Gene-pair detectability across cell types",
    x = NULL,
    y = paste0("Fraction of ", available_pair_n, " evaluable pairs"),
    fill = NULL
  ) +
  theme_nature() +
  theme(legend.position = "top")

p_encoding_scatter <- ggplot(
  group_summary,
  aes(IBI_modal_primary_median, IBI_zero_sensitivity_median)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 0.9, alpha = 0.55, colour = palette_contract["modal"]
  ) +
  coord_equal() +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf(
      "Spearman rho = %.3f\nMedian |difference| = %.3f\nn = %d groups",
      encoding_group_metrics$spearman_rho,
      encoding_group_metrics$median_absolute_difference,
      nrow(group_summary)
    ),
    size = 2.2
  ) +
  labs(
    title = "Sensitivity to the encoding rule",
    x = "IBI: 711-development-set modal encoding",
    y = "IBI: zero encoding sensitivity"
  ) +
  theme_nature()

p_encoding_delta <- ggplot(
  group_summary,
  aes(primary_type_plot, IBI_difference_median)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  geom_boxplot(
    width = 0.62, outlier.shape = NA, linewidth = 0.3,
    fill = "#D9E5EE", colour = palette_contract["zero"]
  ) +
  coord_flip() +
  labs(
    title = "Encoding-related shift by cell type",
    x = NULL,
    y = "Median IBI difference (zero - modal)"
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
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf(
      "Adjusted beta = %.3f (95%% CI %.3f to %.3f)\nP = %.2g",
      encoding_difference_model$standardized_beta,
      encoding_difference_model$ci_lower,
      encoding_difference_model$ci_upper,
      encoding_difference_model$wald_p
    ),
    size = 2.2
  ) +
  labs(
    title = "Encoding sensitivity and double-zero pairs",
    x = "Median fraction of double-zero pairs",
    y = "Absolute median IBI difference"
  ) +
  theme_nature()

coverage_figure <-
  (p_coverage | p_encoding_scatter) /
  (p_encoding_delta | p_zero_dependence) +
  plot_annotation(tag_levels = "A")
save_pub_r(
  coverage_figure,
  "single_cell_coverage_and_encoding_sensitivity_711_modal_primary",
  width_mm = 183, height_mm = 150
)

model_plot_data <- mixed_model_results %>%
  filter(outcome == "z_ibi_modal_primary") %>%
  mutate(
    label = factor(
      label,
      levels = rev(c(
        "Library size", "Detected-gene count", "Double-zero fraction"
      ))
    )
  )
write.csv(
  model_plot_data,
  file.path(output_dir, "source_data_technical_covariate_forest.csv"),
  row.names = FALSE
)

p_technical <- ggplot(
  model_plot_data,
  aes(standardized_beta, label)
) +
  geom_vline(
    xintercept = 0, linetype = "dashed", linewidth = 0.35,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    orientation = "y", width = 0.18, linewidth = 0.55,
    colour = palette_contract["modal"]
  ) +
  geom_point(size = 2.2, colour = palette_contract["modal"]) +
  labs(
    title = "Technical covariates of modal-primary IBI",
    x = "Standardized beta (95% CI)", y = NULL
  ) +
  theme_nature()
save_pub_r(
  p_technical,
  "single_cell_modal_primary_technical_covariates",
  width_mm = 120, height_mm = 72
)

pseudobulk_long <- bind_rows(
  pseudobulk_results %>%
    transmute(
      group_id, sampleName, primary_type,
      encoding_rule = "711-development-set modal encoding",
      cell_level_ibi = IBI_modal_primary_median,
      pseudobulk_ibi = IBI_pseudobulk_modal_primary
    ),
  pseudobulk_results %>%
    transmute(
      group_id, sampleName, primary_type,
      encoding_rule = "Zero encoding sensitivity",
      cell_level_ibi = IBI_zero_sensitivity_median,
      pseudobulk_ibi = IBI_pseudobulk_zero_sensitivity
    )
)
pseudobulk_long$encoding_rule <- factor(
  pseudobulk_long$encoding_rule,
  levels = c(
    "711-development-set modal encoding",
    "Zero encoding sensitivity"
  )
)
write.csv(
  pseudobulk_long,
  file.path(output_dir, "source_data_pseudobulk_scatter.csv"),
  row.names = FALSE
)

pseudobulk_annotations <- pseudobulk_metrics %>%
  transmute(
    encoding_rule = factor(
      scope,
      levels = c(
        "711-development-set modal encoding",
        "Zero encoding sensitivity"
      )
    ),
    label = sprintf(
      "Spearman rho = %.3f\nMAE = %.3f\nn = %d groups",
      spearman_rho, mae, n
    )
  )

p_pb_scatter <- ggplot(
  pseudobulk_long,
  aes(cell_level_ibi, pseudobulk_ibi)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 0.9, alpha = 0.5, colour = palette_contract["modal"]
  ) +
  geom_text(
    data = pseudobulk_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.15, size = 2.2
  ) +
  facet_wrap(
    ~encoding_rule,
    scales = "free",
    labeller = as_labeller(c(
      "711-development-set modal encoding" = "Modal encoding (primary)",
      "Zero encoding sensitivity" = "Zero encoding (sensitivity)"
    ))
  ) +
  labs(
    title = "Cell-level medians versus pseudobulk IBI",
    x = "Median cell-level IBI", y = "Pseudobulk IBI"
  ) +
  theme_nature()

p_pb_celltype <- pseudobulk_metrics_by_celltype %>%
  mutate(
    primary_type = factor(primary_type, levels = rev(celltype_order)),
    scope = factor(
      recode(
        scope,
        "711-development-set modal encoding" = "Modal",
        "Zero encoding sensitivity" = "Zero"
      ),
      levels = c("Modal", "Zero")
    )
  ) %>%
  ggplot(aes(spearman_rho, primary_type, colour = scope, shape = scope)) +
  geom_vline(
    xintercept = 0, linewidth = 0.3,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(
    "Modal" = unname(palette_contract["modal"]),
    "Zero" = unname(palette_contract["zero"])
  )) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(
    title = "Within-cell-type concordance",
    x = "Spearman rho across donors", y = NULL,
    colour = "Encoding", shape = "Encoding"
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
    title = "Pair detectability",
    x = "Median cell-level detectable-pair fraction",
    y = "Pseudobulk detectable-pair fraction"
  ) +
  theme_nature()

p_pb_zero <- ggplot(
  pseudobulk_results,
  aes(both_zero_fraction_median, pseudobulk_both_zero_fraction)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    linewidth = 0.35, colour = palette_contract["neutral_mid"]
  ) +
  geom_point(
    size = 1.0, alpha = 0.5, colour = palette_contract["accent_red"]
  ) +
  labs(
    title = "Double-zero fraction",
    x = "Median cell-level double-zero fraction",
    y = "Pseudobulk double-zero fraction"
  ) +
  theme_nature()

pseudobulk_figure <-
  p_pb_scatter /
  (p_pb_celltype | p_pb_coverage | p_pb_zero) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A")
save_pub_r(
  pseudobulk_figure,
  "single_cell_pseudobulk_concordance_711_modal_primary",
  width_mm = 183, height_mm = 145
)

downsampling_plot_data <- downsampling_metrics %>%
  filter(analysis_level == "Donor-cell-type") %>%
  mutate(
    retained_label = factor(
      paste0(as.integer(retained_fraction * 100), "%"),
      levels = c("25%", "50%", "75%")
    ),
    encoding_rule = factor(
      encoding_rule,
      levels = c(
        "711-development-set modal encoding",
        "Zero encoding sensitivity"
      )
    )
  )
write.csv(
  downsampling_plot_data,
  file.path(output_dir, "source_data_downsampling_panels.csv"),
  row.names = FALSE
)

encoding_colours <- c(
  "711-development-set modal encoding" =
    unname(palette_contract["modal"]),
  "Zero encoding sensitivity" = unname(palette_contract["zero"])
)

p_down_rho <- ggplot(
  downsampling_plot_data,
  aes(retained_label, spearman_rho, colour = encoding_rule)
) +
  geom_hline(
    yintercept = 1, linetype = "dotted", linewidth = 0.35,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_boxplot(
    aes(group = interaction(retained_label, encoding_rule)),
    width = 0.54, outlier.shape = NA, linewidth = 0.35,
    position = position_dodge(width = 0.62)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.06, jitter.height = 0, dodge.width = 0.62
    ),
    size = 1.2, alpha = 0.75
  ) +
  scale_colour_manual(values = encoding_colours) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(
    title = "IBI ranking under UMI thinning",
    subtitle = paste0(
      downsample_repeats,
      " repeats; donor-cell-type medians from a stratified cell subset"
    ),
    x = "Retained UMI fraction",
    y = "Spearman rho versus unthinned IBI",
    colour = NULL
  ) +
  theme_nature() +
  theme(legend.position = "top")

p_down_mae <- ggplot(
  downsampling_plot_data,
  aes(retained_label, mae, colour = encoding_rule)
) +
  geom_hline(
    yintercept = 0, linetype = "dotted", linewidth = 0.35,
    colour = palette_contract["neutral_mid"]
  ) +
  geom_boxplot(
    aes(group = interaction(retained_label, encoding_rule)),
    width = 0.54, outlier.shape = NA, linewidth = 0.35,
    position = position_dodge(width = 0.62)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.06, jitter.height = 0, dodge.width = 0.62
    ),
    size = 1.2, alpha = 0.75
  ) +
  scale_colour_manual(values = encoding_colours) +
  labs(
    title = "Absolute IBI error under UMI thinning",
    subtitle = paste0(
      downsample_repeats,
      " repeats; donor-cell-type medians from a stratified cell subset"
    ),
    x = "Retained UMI fraction",
    y = "MAE versus unthinned IBI",
    colour = NULL
  ) +
  theme_nature() +
  theme(legend.position = "top")

downsampling_figure <- (p_down_rho | p_down_mae) +
  plot_annotation(tag_levels = "A")
save_pub_r(
  downsampling_figure,
  "single_cell_umi_thinning_711_modal_primary",
  width_mm = 183, height_mm = 88
)

# Reviewer-response figure: one panel for each requested validation question.
p_response_a <- p_encoding_scatter +
  labs(title = "Encoding-rule sensitivity")
p_response_b <- p_technical +
  labs(title = "Sequencing depth and gene detection")
p_response_c <- p_down_rho +
  labs(title = "Controlled UMI thinning", subtitle = NULL) +
  theme(legend.position = "top")
p_response_d <- p_pb_scatter +
  labs(title = "Pseudobulk comparison")

response_figure <-
  (p_response_a | p_response_b) /
  (p_response_c | p_response_d) +
  plot_annotation(tag_levels = "A")
save_pub_r(
  response_figure,
  "Response_Figure_R6_modal_primary_single_cell_validation",
  width_mm = 183, height_mm = 150
)

analysis_provenance <- data.frame(
  item = c(
    "Single-cell object", "Cell count", "Donor count",
    "Primary cell-type count", "Feature source", "Reference model",
    "Feature encoding", "Globally available feature genes",
    "Globally evaluable gene pairs", "Primary missing rule",
    "Sensitivity rule", "Pseudobulk unit",
    "Minimum cells per pseudobulk unit", "UMI-thinning subset rule",
    "UMI-thinning fractions", "UMI-thinning repeats", "Random seed",
    "Inference boundary"
  ),
  value = c(
    "Healthy-donor PBMC aging Seurat object",
    nrow(cell_results),
    n_distinct(cell_results$sampleName),
    n_distinct(cell_results$primary_type),
    "Fixed Top390 gene-pair panel selected in the 711-sample development set",
    "Final 711-sample radial epsilon-SVR reference fit",
    "1 when gene A count is greater than gene B count; otherwise 0",
    paste0(
      sum(all_feature_genes %in% rownames(counts)),
      "/", length(all_feature_genes)
    ),
    paste0(sum(pair_table$pair_present), "/", nrow(pair_table)),
    paste0(
      "Globally unavailable or double-zero pairs replaced by the ",
      "pair-specific modal state from the 711-sample development set"
    ),
    "Globally unavailable or double-zero pairs encoded as literal zero",
    "Donor by primary cell type",
    minimum_cells_per_donor_celltype,
    paste0(
      "Up to ", downsample_cells_per_donor_celltype,
      " cells sampled per eligible donor-by-primary-cell-type group"
    ),
    paste(downsample_fractions, collapse = "; "),
    downsample_repeats,
    base_seed,
    paste0(
      "Cells are technical observations; donor or donor-by-cell-type ",
      "summaries define the inferential interpretation"
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(
  analysis_provenance,
  file.path(output_dir, "analysis_provenance.csv"),
  row.names = FALSE
)

legend_text <- c(
  "Response Figure R6 | Validation of the fixed bulk-trained IBI model in single-cell RNA-seq data.",
  paste0(
    "(A) Comparison of donor-by-primary-cell-type median IBI values under ",
    "the primary pair-specific modal encoding derived from the 711-sample ",
    "development set and the literal-zero sensitivity rule."
  ),
  paste0(
    "(B) Standardized coefficients and 95% Wald confidence intervals from ",
    "separate linear mixed-effects models relating modal-primary IBI to ",
    "library size, detected-gene count or double-zero fraction. Models were ",
    "adjusted for primary cell type and included donor as a random intercept."
  ),
  paste0(
    "(C) Spearman concordance between unthinned and UMI-thinned ",
    "donor-by-primary-cell-type median IBI estimates at 25%, 50% and 75% ",
    "retained depth under both encoding rules. Ten independent thinning ",
    "repeats were performed at each depth."
  ),
  paste0(
    "(D) Comparison of median cell-level and pseudobulk IBI estimates under ",
    "both encoding rules. Pseudobulk profiles were constructed by summing ",
    "raw UMI counts within donor-by-primary-cell-type groups containing at ",
    "least 50 cells. Correlations and errors are descriptive across groups."
  ),
  paste0(
    "Across all panels, cells were treated as technical observations; the ",
    "61 donors and 475 eligible donor-by-cell-type groups define the ",
    "biological and grouped analysis structure, respectively."
  )
)
writeLines(
  legend_text,
  file.path(output_dir, "response_figure_legend_draft.txt"),
  useBytes = TRUE
)

qa_notes <- c(
  "Core conclusion: relative IBI patterns are evaluated under modal-primary encoding, controlled UMI thinning, and pseudobulk aggregation.",
  "Archetype: quantitative grid.",
  "Backend: R only.",
  "Final width: 183 mm for multi-panel figures.",
  paste0("Independent donor count: ", n_distinct(group_summary$sampleName), "."),
  paste0("Eligible donor-by-cell-type groups: ", nrow(group_summary), "."),
  "Cell-level estimates are descriptive technical observations, not independent biological replicates.",
  "Mixed models adjust for primary cell type and include donor as a random intercept.",
  "All quantitative panels have CSV source data in the same output directory.",
  "SVG and PDF retain editable text; TIFF is exported at 600 dpi; PNG is a QA preview."
)
writeLines(
  qa_notes,
  file.path(output_dir, "figure_qa_notes.txt"),
  useBytes = TRUE
)

save(
  pair_table, group_summary, celltype_summary, encoding_metrics,
  mixed_model_results, pseudobulk_results, pseudobulk_metrics,
  pseudobulk_metrics_by_celltype, downsampling_metrics,
  downsampling_summary,
  file = file.path(
    output_dir, "single_cell_modal_primary_analysis_objects.rda"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "R_sessionInfo.txt"),
  useBytes = TRUE
)

message("Modal-primary single-cell sensitivity analysis complete: ", output_dir)
