options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(e1071)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

# Reviewer-facing BAS/CAS/IAS robustness analysis.
# The feature set is the exploratory top390 selected from all available cohorts.
# IBI is predicted with the final fixed-feature SVR trained in 781 development
# samples. The irAE cohort is prediction-only.

root <- Sys.getenv("IBI_TOOL_ROOT", unset = "C:/temp/ibi_tool")
feature_file <- file.path(
  root, "Result", "revision_781_fixed_top390_final",
  "fixed_top390_all_samples_781_development.rda"
)
model_file <- file.path(
  root, "Result", "revision_781_final_svr_external_validation_20261028_corrected",
  "final_781_fixed_top390_svr_model_and_predictions.rda"
)
expression_file <- file.path(root, "Result", "validation", "GSE186143_bulk_TPM.txt")
metadata_file <- file.path(root, "Result", "validation", "GSE186143_sam_info.txt")
output_dir <- file.path(root, "Result", "revision_781_BAS_CAS_IAS_sensitivity_top390_20261028")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(feature_file, model_file, expression_file, metadata_file)
if (!all(file.exists(required_files))) {
  stop("Missing input file(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

feature_environment <- new.env(parent = emptyenv())
load(feature_file, envir = feature_environment)
feature_reference <- list(
  features_top390 = feature_environment$features_top390,
  top390 = feature_environment$top390,
  sample_manifest = feature_environment$sample_manifest
)
model_environment <- new.env(parent = emptyenv())
load(model_file, envir = model_environment)
features_top390 <- model_environment$features_top390
final_svr_model <- model_environment$final_svr_model
final_prepared <- model_environment$final_prepared
development_predictions <- model_environment$development_predictions

stopifnot(
  length(feature_reference$features_top390) == 390L,
  identical(feature_reference$features_top390, features_top390),
  inherits(final_svr_model, "svm"),
  nrow(development_predictions) == 781L,
  length(final_prepared$modes) == 390L,
  identical(final_prepared$kept_features, features_top390)
)

expression_all <- read.delim(
  expression_file, row.names = 1, check.names = FALSE, quote = "",
  comment.char = ""
)
expression_all <- as.matrix(expression_all)
storage.mode(expression_all) <- "double"
metadata_all <- read.delim(metadata_file, row.names = 1, check.names = FALSE)

grade_text <- as.character(metadata_all[colnames(expression_all), "Grade"])
grade_number <- suppressWarnings(as.numeric(sub("^Grade", "", grade_text)))
analysis_samples <- colnames(expression_all)[is.finite(grade_number)]
excluded_samples <- colnames(expression_all)[!is.finite(grade_number)]
expression <- expression_all[, analysis_samples, drop = FALSE]
grade_text <- grade_text[is.finite(grade_number)]
grade_number <- grade_number[is.finite(grade_number)]

stopifnot(
  ncol(expression) == 60L,
  all(grade_number %in% 0:4),
  identical(colnames(expression), analysis_samples)
)

make_pair_matrix <- function(expression_matrix, pair_table, pair_names) {
  result <- matrix(
    NA_real_, nrow = ncol(expression_matrix), ncol = length(pair_names),
    dimnames = list(colnames(expression_matrix), pair_names)
  )
  available <- logical(length(pair_names))
  tie_count <- integer(length(pair_names))
  for (i in seq_along(pair_names)) {
    gene_a <- pair_table$geneA[i]
    gene_b <- pair_table$geneB[i]
    if (gene_a %in% rownames(expression_matrix) && gene_b %in% rownames(expression_matrix)) {
      a <- expression_matrix[gene_a, ]
      b <- expression_matrix[gene_b, ]
      observed <- is.finite(a) & is.finite(b)
      result[observed, i] <- as.numeric(a[observed] > b[observed])
      available[i] <- all(observed)
      tie_count[i] <- sum(observed & a == b)
    }
  }
  list(matrix = result, available = available, tie_count = tie_count)
}

pair_result <- make_pair_matrix(
  expression, feature_reference$top390,
  feature_reference$features_top390
)
pair_matrix <- pair_result$matrix
apply_preprocessing <- function(x, preprocessing) {
  x <- as.matrix(x)
  for (j in seq_len(ncol(x))) {
    x[!is.finite(x[, j]), j] <- preprocessing$modes[j]
  }
  x <- x[, preprocessing$kept_features, drop = FALSE]
  sweep(
    sweep(x, 2L, preprocessing$centers, "-"),
    2L, preprocessing$scales, "/"
  )
}
pair_matrix_scaled <- apply_preprocessing(pair_matrix, final_prepared)
ibi <- as.numeric(predict(final_svr_model, pair_matrix_scaled))
names(ibi) <- rownames(pair_matrix)

sample_data <- data.frame(
  sample_id = analysis_samples,
  grade = grade_text,
  grade_number = grade_number,
  IBI = ibi,
  stringsAsFactors = FALSE
)
write.csv(sample_data, file.path(output_dir, "GSE186143_top390_IBI_and_grade.csv"), row.names = FALSE)

pair_coverage <- data.frame(
  gene_pair = feature_reference$features_top390,
  geneA = feature_reference$top390$geneA,
  geneB = feature_reference$top390$geneB,
  both_genes_available = pair_result$available,
  tie_count = pair_result$tie_count,
  tie_fraction = pair_result$tie_count / nrow(pair_matrix),
  stringsAsFactors = FALSE
)
write.csv(pair_coverage, file.path(output_dir, "GSE186143_top390_pair_coverage.csv"), row.names = FALSE)

# Match the legacy expression filter, then standardize every modeled quantity.
expression_filter <- rowMeans(expression > 0.5, na.rm = TRUE) > 0.20
expression_log <- log2(expression[expression_filter, , drop = FALSE] + 1)
gene_sd <- apply(expression_log, 1L, sd, na.rm = TRUE)
expression_log <- expression_log[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
gene_mean <- rowMeans(expression_log)
gene_sd <- apply(expression_log, 1L, sd)
expression_z <- sweep(expression_log, 1L, gene_mean, "-")
expression_z <- sweep(expression_z, 1L, gene_sd, "/")

z_score <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s <= 0) stop("Cannot standardize a constant or invalid vector")
  (x - mean(x)) / s
}

ibi_z <- z_score(ibi)
grade_z <- z_score(grade_number)
interaction_z <- z_score(ibi_z * grade_z)

design_simple <- cbind(Intercept = 1, Grade = grade_z)
design_adjusted <- cbind(Intercept = 1, IBI = ibi_z, Grade = grade_z)
design_full <- cbind(
  Intercept = 1, IBI = ibi_z, Grade = grade_z,
  Interaction = interaction_z
)

fit_limma <- function(y, design) {
  fit <- limma::lmFit(y, design)
  limma::eBayes(fit)
}

fit_simple <- fit_limma(expression_z, design_simple)
fit_adjusted <- fit_limma(expression_z, design_adjusted)
fit_full <- fit_limma(expression_z, design_full)

genes <- rownames(expression_z)
component_results <- data.frame(
  gene = genes,
  mean_log2_TPM = gene_mean[genes],
  detection_fraction_TPM_gt_0_5 = rowMeans(expression[genes, , drop = FALSE] > 0.5),
  beta_BAS = fit_full$coefficients[, "IBI"],
  beta_CAS = fit_full$coefficients[, "Grade"],
  beta_IAS = fit_full$coefficients[, "Interaction"],
  BAS = abs(fit_full$coefficients[, "IBI"]),
  CAS = abs(fit_full$coefficients[, "Grade"]),
  IAS = abs(fit_full$coefficients[, "Interaction"]),
  t_BAS = fit_full$t[, "IBI"],
  t_CAS = fit_full$t[, "Grade"],
  t_IAS = fit_full$t[, "Interaction"],
  p_BAS = fit_full$p.value[, "IBI"],
  p_CAS = fit_full$p.value[, "Grade"],
  p_IAS = fit_full$p.value[, "Interaction"],
  stringsAsFactors = FALSE
)
component_results$FDR_BAS <- p.adjust(component_results$p_BAS, method = "BH")
component_results$FDR_CAS <- p.adjust(component_results$p_CAS, method = "BH")
component_results$FDR_IAS <- p.adjust(component_results$p_IAS, method = "BH")

classify_dominance <- function(result_table, threshold) {
  score <- as.matrix(result_table[, c("BAS", "CAS", "IAS")])
  output <- rep("None", nrow(score))
  for (component in colnames(score)) {
    others <- setdiff(colnames(score), component)
    dominant <- score[, component] > threshold * score[, others[1L]] &
      score[, component] > threshold * score[, others[2L]]
    output[dominant] <- component
  }
  output
}

thresholds <- c(1.2, 1.5, 2.0)
threshold_labels <- paste0("dominance_", gsub("\\.", "_", thresholds))
for (i in seq_along(thresholds)) {
  component_results[[threshold_labels[i]]] <- classify_dominance(
    component_results, thresholds[i]
  )
}

write.csv(
  component_results,
  file.path(output_dir, "standardized_BAS_CAS_IAS_gene_level.csv"),
  row.names = FALSE
)

dominance_summary <- do.call(rbind, lapply(seq_along(thresholds), function(i) {
  classification <- component_results[[threshold_labels[i]]]
  do.call(rbind, lapply(c("BAS", "CAS", "IAS", "None"), function(component) {
    selected <- classification == component
    significant <- if (component == "None") {
      rep(FALSE, nrow(component_results))
    } else {
      component_results[[paste0("FDR_", component)]] < 0.05
    }
    data.frame(
      threshold = thresholds[i], component = component,
      n_genes = sum(selected),
      n_genes_component_FDR_lt_0_05 = sum(selected & significant),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(dominance_summary, file.path(output_dir, "dominance_threshold_counts.csv"), row.names = FALSE)

jaccard <- function(a, b) {
  union_size <- length(union(a, b))
  if (!union_size) return(NA_real_)
  length(intersect(a, b)) / union_size
}

threshold_jaccard <- do.call(rbind, lapply(c("BAS", "CAS", "IAS"), function(component) {
  sets <- lapply(seq_along(thresholds), function(i) {
    component_results$gene[component_results[[threshold_labels[i]]] == component]
  })
  do.call(rbind, lapply(seq_along(thresholds), function(i) {
    do.call(rbind, lapply(seq_along(thresholds), function(j) {
      data.frame(
        component = component,
        threshold_1 = thresholds[i], threshold_2 = thresholds[j],
        n_set_1 = length(sets[[i]]), n_set_2 = length(sets[[j]]),
        intersection = length(intersect(sets[[i]], sets[[j]])),
        union = length(union(sets[[i]], sets[[j]])),
        jaccard = jaccard(sets[[i]], sets[[j]]),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
write.csv(threshold_jaccard, file.path(output_dir, "dominance_threshold_jaccard.csv"), row.names = FALSE)

# Conventional model comparison: condition-only, IBI-adjusted, and interaction.
coefficient_table <- data.frame(
  gene = genes,
  beta_grade_condition_only = fit_simple$coefficients[, "Grade"],
  p_grade_condition_only = fit_simple$p.value[, "Grade"],
  beta_ibi_adjusted = fit_adjusted$coefficients[, "IBI"],
  p_ibi_adjusted = fit_adjusted$p.value[, "IBI"],
  beta_grade_ibi_adjusted = fit_adjusted$coefficients[, "Grade"],
  p_grade_ibi_adjusted = fit_adjusted$p.value[, "Grade"],
  beta_ibi_full = fit_full$coefficients[, "IBI"],
  p_ibi_full = fit_full$p.value[, "IBI"],
  beta_grade_full = fit_full$coefficients[, "Grade"],
  p_grade_full = fit_full$p.value[, "Grade"],
  beta_interaction_full = fit_full$coefficients[, "Interaction"],
  p_interaction_full = fit_full$p.value[, "Interaction"],
  stringsAsFactors = FALSE
)
for (p_column in grep("^p_", names(coefficient_table), value = TRUE)) {
  coefficient_table[[sub("^p_", "FDR_", p_column)]] <- p.adjust(
    coefficient_table[[p_column]], method = "BH"
  )
}

calculate_fit_metrics <- function(y, design, coefficients) {
  fitted <- coefficients %*% t(design)
  rss <- rowSums((y - fitted)^2)
  tss <- rowSums((y - rowMeans(y))^2)
  r_squared <- 1 - rss / tss
  n <- ncol(y)
  predictors <- ncol(design) - 1L
  adjusted_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - predictors - 1)
  list(rss = rss, r_squared = r_squared, adjusted_r_squared = adjusted_r_squared)
}

metrics_simple <- calculate_fit_metrics(expression_z, design_simple, fit_simple$coefficients)
metrics_adjusted <- calculate_fit_metrics(expression_z, design_adjusted, fit_adjusted$coefficients)
metrics_full <- calculate_fit_metrics(expression_z, design_full, fit_full$coefficients)

nested_p_value <- function(rss_reduced, rss_full, parameters_reduced, parameters_full, n) {
  df_num <- parameters_full - parameters_reduced
  df_den <- n - parameters_full
  f_stat <- ((rss_reduced - rss_full) / df_num) / (rss_full / df_den)
  f_stat[!is.finite(f_stat) | f_stat < 0] <- 0
  stats::pf(f_stat, df_num, df_den, lower.tail = FALSE)
}

p_add_ibi <- nested_p_value(
  metrics_simple$rss, metrics_adjusted$rss,
  ncol(design_simple), ncol(design_adjusted), ncol(expression_z)
)
p_add_interaction <- nested_p_value(
  metrics_adjusted$rss, metrics_full$rss,
  ncol(design_adjusted), ncol(design_full), ncol(expression_z)
)

model_fit_table <- data.frame(
  gene = genes,
  R2_condition_only = metrics_simple$r_squared,
  adjusted_R2_condition_only = metrics_simple$adjusted_r_squared,
  R2_IBI_adjusted = metrics_adjusted$r_squared,
  adjusted_R2_IBI_adjusted = metrics_adjusted$adjusted_r_squared,
  R2_full_interaction = metrics_full$r_squared,
  adjusted_R2_full_interaction = metrics_full$adjusted_r_squared,
  delta_R2_add_IBI = metrics_adjusted$r_squared - metrics_simple$r_squared,
  p_add_IBI = p_add_ibi,
  FDR_add_IBI = p.adjust(p_add_ibi, method = "BH"),
  delta_R2_add_interaction = metrics_full$r_squared - metrics_adjusted$r_squared,
  p_add_interaction = p_add_interaction,
  FDR_add_interaction = p.adjust(p_add_interaction, method = "BH"),
  stringsAsFactors = FALSE
)
write.csv(coefficient_table, file.path(output_dir, "conventional_model_coefficients.csv"), row.names = FALSE)
write.csv(model_fit_table, file.path(output_dir, "conventional_model_fit_metrics.csv"), row.names = FALSE)

top_genes <- function(score, n = 100L) {
  names(sort(setNames(score, genes), decreasing = TRUE))[seq_len(min(n, length(score)))]
}

top100_simple_grade <- top_genes(abs(fit_simple$coefficients[, "Grade"]))
top100_adjusted_grade <- top_genes(abs(fit_adjusted$coefficients[, "Grade"]))
top100_full_grade <- top_genes(component_results$CAS)

conventional_summary <- data.frame(
  comparison = c(
    "Condition-only Grade vs IBI-adjusted Grade",
    "Condition-only Grade vs full-model CAS",
    "IBI-adjusted Grade vs full-model CAS"
  ),
  spearman_signed_beta = c(
    cor(fit_simple$coefficients[, "Grade"], fit_adjusted$coefficients[, "Grade"], method = "spearman"),
    cor(fit_simple$coefficients[, "Grade"], fit_full$coefficients[, "Grade"], method = "spearman"),
    cor(fit_adjusted$coefficients[, "Grade"], fit_full$coefficients[, "Grade"], method = "spearman")
  ),
  spearman_absolute_beta = c(
    cor(abs(fit_simple$coefficients[, "Grade"]), abs(fit_adjusted$coefficients[, "Grade"]), method = "spearman"),
    cor(abs(fit_simple$coefficients[, "Grade"]), component_results$CAS, method = "spearman"),
    cor(abs(fit_adjusted$coefficients[, "Grade"]), component_results$CAS, method = "spearman")
  ),
  top100_jaccard = c(
    jaccard(top100_simple_grade, top100_adjusted_grade),
    jaccard(top100_simple_grade, top100_full_grade),
    jaccard(top100_adjusted_grade, top100_full_grade)
  ),
  stringsAsFactors = FALSE
)
write.csv(conventional_summary, file.path(output_dir, "conventional_model_comparison_summary.csv"), row.names = FALSE)

nested_summary <- data.frame(
  added_term = c("IBI", "IBI-by-Grade interaction"),
  reduced_model = c("Expression ~ Grade", "Expression ~ IBI + Grade"),
  full_model = c("Expression ~ IBI + Grade", "Expression ~ IBI + Grade + IBI:Grade"),
  n_genes_FDR_lt_0_05 = c(sum(model_fit_table$FDR_add_IBI < 0.05), sum(model_fit_table$FDR_add_interaction < 0.05)),
  median_delta_R2 = c(median(model_fit_table$delta_R2_add_IBI), median(model_fit_table$delta_R2_add_interaction)),
  q1_delta_R2 = c(quantile(model_fit_table$delta_R2_add_IBI, 0.25), quantile(model_fit_table$delta_R2_add_interaction, 0.25)),
  q3_delta_R2 = c(quantile(model_fit_table$delta_R2_add_IBI, 0.75), quantile(model_fit_table$delta_R2_add_interaction, 0.75)),
  stringsAsFactors = FALSE
)
write.csv(nested_summary, file.path(output_dir, "nested_model_comparison_summary.csv"), row.names = FALSE)

# Exclude every gene directly represented by the top390 gene-pair features.
ibi_input_genes <- unique(c(feature_reference$top390$geneA, feature_reference$top390$geneB))
ibi_input_genes_in_analysis <- intersect(ibi_input_genes, component_results$gene)
excluded_component_results <- component_results[
  !component_results$gene %in% ibi_input_genes, , drop = FALSE
]
write.csv(
  data.frame(gene = sort(ibi_input_genes), in_decomposition_analysis = sort(ibi_input_genes) %in% genes),
  file.path(output_dir, "top390_unique_input_genes.csv"), row.names = FALSE
)
write.csv(
  excluded_component_results,
  file.path(output_dir, "standardized_BAS_CAS_IAS_excluding_IBI_input_genes.csv"),
  row.names = FALSE
)

exclusion_top100_summary <- do.call(rbind, lapply(c("BAS", "CAS", "IAS"), function(component) {
  full_rank <- component_results$gene[order(component_results[[component]], decreasing = TRUE)]
  excluded_rank <- excluded_component_results$gene[
    order(excluded_component_results[[component]], decreasing = TRUE)
  ]
  full_top <- head(full_rank, 100L)
  excluded_top <- head(excluded_rank, 100L)
  shared_genes <- intersect(component_results$gene, excluded_component_results$gene)
  data.frame(
    component = component,
    n_unique_IBI_input_genes = length(ibi_input_genes),
    n_IBI_input_genes_in_analysis = length(ibi_input_genes_in_analysis),
    n_removed_from_full_top100 = sum(full_top %in% ibi_input_genes),
    top100_intersection = length(intersect(full_top, excluded_top)),
    top100_jaccard = jaccard(full_top, excluded_top),
    shared_gene_score_spearman = cor(
      component_results[[component]][match(shared_genes, component_results$gene)],
      excluded_component_results[[component]][match(shared_genes, excluded_component_results$gene)],
      method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}))
write.csv(exclusion_top100_summary, file.path(output_dir, "IBI_gene_exclusion_top100_summary.csv"), row.names = FALSE)

exclusion_dominance_summary <- do.call(rbind, lapply(seq_along(thresholds), function(i) {
  do.call(rbind, lapply(c("BAS", "CAS", "IAS"), function(component) {
    full_set <- component_results$gene[
      component_results[[threshold_labels[i]]] == component
    ]
    excluded_set <- excluded_component_results$gene[
      excluded_component_results[[threshold_labels[i]]] == component
    ]
    data.frame(
      threshold = thresholds[i], component = component,
      full_n = length(full_set), excluded_n = length(excluded_set),
      removed_IBI_input_genes = sum(full_set %in% ibi_input_genes),
      jaccard = jaccard(full_set, excluded_set),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(exclusion_dominance_summary, file.path(output_dir, "IBI_gene_exclusion_dominance_summary.csv"), row.names = FALSE)

# Genuine null-hypothesis permutation analysis.
# BAS null: permute IBI labels while retaining Grade and expression.
# CAS/IAS null: permute Grade labels while retaining IBI and expression.
# The statistic is the mean absolute standardized coefficient among the top 100 genes.
coefficient_for_term <- function(y, ibi_vector, grade_vector, term) {
  interaction_vector <- z_score(ibi_vector * grade_vector)
  design <- cbind(
    Intercept = 1, IBI = ibi_vector, Grade = grade_vector,
    Interaction = interaction_vector
  )
  target <- design[, term]
  other <- design[, setdiff(colnames(design), term), drop = FALSE]
  residual_target <- qr.resid(qr(other), target)
  denominator <- sum(residual_target^2)
  if (!is.finite(denominator) || denominator <= 1e-12) {
    return(rep(NA_real_, nrow(y)))
  }
  as.numeric(y %*% residual_target / denominator)
}

n_permutations <- as.integer(Sys.getenv("BAS_N_PERMUTATIONS", unset = "1000"))
permutation_seed <- 20261024L
set.seed(permutation_seed)
ibi_permutation_index <- replicate(n_permutations, sample.int(length(ibi_z)))
grade_permutation_index <- replicate(n_permutations, sample.int(length(grade_z)))

observed_top_indices <- list(
  BAS = order(component_results$BAS, decreasing = TRUE)[1:100],
  CAS = order(component_results$CAS, decreasing = TRUE)[1:100],
  IAS = order(component_results$IAS, decreasing = TRUE)[1:100]
)
observed_statistic <- c(
  BAS = mean(component_results$BAS[observed_top_indices$BAS]),
  CAS = mean(component_results$CAS[observed_top_indices$CAS]),
  IAS = mean(component_results$IAS[observed_top_indices$IAS])
)

one_permutation <- function(iteration, y, ibi_vector, grade_vector,
                            ibi_indices, grade_indices, observed_indices) {
  permuted_ibi <- ibi_vector[ibi_indices[, iteration]]
  permuted_grade <- grade_vector[grade_indices[, iteration]]
  beta_bas <- abs(coefficient_for_term(y, permuted_ibi, grade_vector, "IBI"))
  beta_cas <- abs(coefficient_for_term(y, ibi_vector, permuted_grade, "Grade"))
  beta_ias <- abs(coefficient_for_term(y, ibi_vector, permuted_grade, "Interaction"))
  top_bas <- order(beta_bas, decreasing = TRUE)[1:100]
  top_cas <- order(beta_cas, decreasing = TRUE)[1:100]
  top_ias <- order(beta_ias, decreasing = TRUE)[1:100]
  c(
    BAS_statistic = mean(beta_bas[top_bas]),
    CAS_statistic = mean(beta_cas[top_cas]),
    IAS_statistic = mean(beta_ias[top_ias]),
    BAS_jaccard = jaccard(observed_indices$BAS, top_bas),
    CAS_jaccard = jaccard(observed_indices$CAS, top_cas),
    IAS_jaccard = jaccard(observed_indices$IAS, top_ias)
  )
}

permutation_checkpoint <- file.path(output_dir, "permutation_null_1000_checkpoint.rds")
if (file.exists(permutation_checkpoint)) {
  permutation_results <- readRDS(permutation_checkpoint)
  if (!is.matrix(permutation_results) || ncol(permutation_results) != 6L ||
      nrow(permutation_results) > n_permutations) {
    stop("Existing permutation checkpoint is incompatible with this analysis")
  }
} else {
  permutation_results <- matrix(numeric(), nrow = 0L, ncol = 6L)
  colnames(permutation_results) <- c(
    "BAS_statistic", "CAS_statistic", "IAS_statistic",
    "BAS_jaccard", "CAS_jaccard", "IAS_jaccard"
  )
}

n_workers_requested <- as.integer(Sys.getenv("BAS_N_WORKERS", unset = "4"))
n_workers <- max(1L, min(n_workers_requested, parallel::detectCores(logical = TRUE)))
cluster <- parallel::makeCluster(n_workers)
parallel::clusterExport(
  cluster,
  c("expression_z", "ibi_z", "grade_z", "ibi_permutation_index",
    "grade_permutation_index", "observed_top_indices", "z_score",
    "coefficient_for_term", "jaccard", "one_permutation"),
  envir = environment()
)

remaining <- if (nrow(permutation_results) < n_permutations) {
  seq.int(nrow(permutation_results) + 1L, n_permutations)
} else {
  integer()
}
if (length(remaining)) {
  batches <- split(remaining, ceiling(seq_along(remaining) / 100L))
  for (batch in batches) {
    batch_result <- parallel::parLapply(
      cluster, batch, one_permutation,
      y = expression_z, ibi_vector = ibi_z, grade_vector = grade_z,
      ibi_indices = ibi_permutation_index,
      grade_indices = grade_permutation_index,
      observed_indices = observed_top_indices
    )
    batch_matrix <- do.call(rbind, batch_result)
    permutation_results <- rbind(permutation_results, batch_matrix)
    saveRDS(permutation_results, permutation_checkpoint)
    cat("Completed", nrow(permutation_results), "of", n_permutations, "permutations\n")
  }
}
parallel::stopCluster(cluster)

permutation_long <- do.call(rbind, lapply(c("BAS", "CAS", "IAS"), function(component) {
  data.frame(
    permutation = seq_len(n_permutations), component = component,
    null_statistic = permutation_results[, paste0(component, "_statistic")],
    observed_top100_jaccard = permutation_results[, paste0(component, "_jaccard")],
    stringsAsFactors = FALSE
  )
}))
write.csv(permutation_long, file.path(output_dir, "permutation_null_1000.csv"), row.names = FALSE)

permutation_summary <- do.call(rbind, lapply(c("BAS", "CAS", "IAS"), function(component) {
  null_value <- permutation_results[, paste0(component, "_statistic")]
  null_jaccard <- permutation_results[, paste0(component, "_jaccard")]
  data.frame(
    component = component,
    hypothesis = if (component == "BAS") {
      "IBI labels are exchangeable with respect to expression conditional on Grade"
    } else {
      "irAE Grade labels are exchangeable with respect to expression conditional on IBI"
    },
    observed_top100_mean_absolute_beta = observed_statistic[component],
    null_median = median(null_value),
    null_q025 = quantile(null_value, 0.025),
    null_q975 = quantile(null_value, 0.975),
    empirical_p = (1 + sum(null_value >= observed_statistic[component])) /
      (n_permutations + 1),
    median_observed_vs_null_top100_jaccard = median(null_jaccard),
    q025_observed_vs_null_top100_jaccard = quantile(null_jaccard, 0.025),
    q975_observed_vs_null_top100_jaccard = quantile(null_jaccard, 0.975),
    n_permutations = n_permutations,
    stringsAsFactors = FALSE
  )
}))
write.csv(permutation_summary, file.path(output_dir, "permutation_null_summary.csv"), row.names = FALSE)

# Optional offline GO-BP stability analysis for the top 100 dominant genes.
go_results <- data.frame()
go_stability <- data.frame()
run_go <- !identical(Sys.getenv("BAS_SKIP_GO"), "1")
if (!run_go) {
  existing_go <- file.path(output_dir, "GO_BP_threshold_enrichment.csv")
  existing_go_stability <- file.path(output_dir, "GO_BP_threshold_stability.csv")
  if (file.exists(existing_go)) go_results <- read.csv(existing_go, check.names = FALSE)
  if (file.exists(existing_go_stability)) {
    go_stability <- read.csv(existing_go_stability, check.names = FALSE)
  }
}
if (run_go &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  gene_map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, keys = unique(genes),
    columns = "ENTREZID", keytype = "SYMBOL"
  ))
  gene_map <- gene_map[!is.na(gene_map$ENTREZID), , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map$SYMBOL), , drop = FALSE]
  universe_entrez <- unique(gene_map$ENTREZID)
  enrichment_sets <- list()
  for (component in c("BAS", "CAS", "IAS")) {
    for (i in seq_along(thresholds)) {
      dominant_genes <- component_results$gene[
        component_results[[threshold_labels[i]]] == component
      ]
      dominant_genes <- dominant_genes[
        order(component_results[[component]][match(dominant_genes, component_results$gene)],
              decreasing = TRUE)
      ]
      selected_genes <- head(dominant_genes, 100L)
      selected_entrez <- unique(gene_map$ENTREZID[match(selected_genes, gene_map$SYMBOL)])
      selected_entrez <- selected_entrez[!is.na(selected_entrez)]
      key <- paste(component, thresholds[i], sep = "_")
      if (length(selected_entrez) >= 10L) {
        enrichment <- suppressMessages(clusterProfiler::enrichGO(
          gene = selected_entrez, universe = universe_entrez,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1,
          qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
          readable = TRUE
        ))
        enrichment_table <- as.data.frame(enrichment)
      } else {
        enrichment_table <- data.frame()
      }
      if (nrow(enrichment_table)) {
        enrichment_table$component <- component
        enrichment_table$threshold <- thresholds[i]
        enrichment_table$n_dominant_genes <- length(dominant_genes)
        enrichment_table$n_selected_genes <- length(selected_genes)
        go_results <- rbind(go_results, enrichment_table)
        enrichment_sets[[key]] <- enrichment_table$ID[enrichment_table$p.adjust < 0.05]
      } else {
        enrichment_sets[[key]] <- character()
      }
    }
  }
  if (nrow(go_results)) {
    write.csv(go_results, file.path(output_dir, "GO_BP_threshold_enrichment.csv"), row.names = FALSE)
  }
  go_stability <- do.call(rbind, lapply(c("BAS", "CAS", "IAS"), function(component) {
    do.call(rbind, lapply(seq_along(thresholds), function(i) {
      do.call(rbind, lapply(seq_along(thresholds), function(j) {
        set_i <- enrichment_sets[[paste(component, thresholds[i], sep = "_")]]
        set_j <- enrichment_sets[[paste(component, thresholds[j], sep = "_")]]
        data.frame(
          component = component,
          threshold_1 = thresholds[i], threshold_2 = thresholds[j],
          n_terms_1 = length(set_i), n_terms_2 = length(set_j),
          jaccard_FDR_lt_0_05 = jaccard(set_i, set_j),
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  write.csv(go_stability, file.path(output_dir, "GO_BP_threshold_stability.csv"), row.names = FALSE)
}

# Figure contract:
# Core conclusion: standardized decomposition is robust to reasonable dominance
# thresholds and exclusion of direct IBI-input genes, with observed effect
# strength evaluated against component-specific permutation nulls.
# Archetype: quantitative robustness grid; target: Genome Medicine revision.
# Backend/export: R-only, double-column PDF, sans-serif, editable vector text.
component_palette <- c(BAS = "#3B6FB6", CAS = "#D9822B", IAS = "#2A9D8F")
theme_revision <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      plot.tag = element_text(face = "bold", size = 9),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank()
    )
}

dominance_plot_data <- dominance_summary[dominance_summary$component != "None", ]
dominance_plot_data$threshold_label <- factor(
  paste0(dominance_plot_data$threshold, "-fold"),
  levels = paste0(thresholds, "-fold")
)
p_a <- ggplot(dominance_plot_data, aes(threshold_label, n_genes, fill = component)) +
  geom_col(
    aes(group = component), position = position_dodge(width = 0.75),
    width = 0.68, alpha = 0.28
  ) +
  geom_col(
    aes(y = n_genes_component_FDR_lt_0_05, group = component),
    position = position_dodge(width = 0.75), width = 0.68
  ) +
  scale_fill_manual(values = component_palette) +
  labs(
    x = "Dominance threshold", y = "Dominant genes", fill = NULL,
    caption = "Pale bars: all classified genes; solid bars: component FDR < 0.05"
  ) +
  theme_revision() +
  theme(
    legend.position = "top",
    plot.caption = element_text(size = 5.5, hjust = 0, colour = "#4A4A4A")
  )

threshold_reference <- threshold_jaccard[
  threshold_jaccard$threshold_2 == 1.5, , drop = FALSE
]
threshold_reference$threshold_label <- factor(
  paste0(threshold_reference$threshold_1, "-fold"),
  levels = paste0(thresholds, "-fold")
)
p_b <- ggplot(
  threshold_reference,
  aes(threshold_label, jaccard, colour = component, group = component)
) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "#B5B5B5") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = component_palette) +
  scale_y_continuous(limits = c(0, 1.03), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Compared with the 1.5-fold set", y = "Jaccard index", colour = NULL) +
  theme_revision() +
  theme(legend.position = "none")

permutation_plot_data <- permutation_long
permutation_plot_data$component <- factor(
  permutation_plot_data$component, levels = c("BAS", "CAS", "IAS")
)
observed_plot_data <- permutation_summary
observed_plot_data$component <- factor(observed_plot_data$component, levels = c("BAS", "CAS", "IAS"))
observed_plot_data$label <- paste0(
  "Empirical P = ", formatC(observed_plot_data$empirical_p, format = "g", digits = 2)
)
p_c <- ggplot(permutation_plot_data, aes(null_statistic, fill = component)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, alpha = 0.75, colour = "white", linewidth = 0.15) +
  geom_vline(
    data = observed_plot_data,
    aes(xintercept = observed_top100_mean_absolute_beta),
    linewidth = 0.65, linetype = "dashed", colour = "black"
  ) +
  geom_text(
    data = observed_plot_data,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.08, vjust = 1.35, size = 2.2
  ) +
  facet_wrap(~ component, scales = "free", nrow = 1) +
  scale_fill_manual(values = component_palette, guide = "none") +
  labs(x = "Top-100 mean |standardized beta| under the null", y = "Density") +
  theme_revision()

exclusion_caption <- paste0(
  "Removed IBI-input genes:\n",
  paste0(
    exclusion_top100_summary$component, " ",
    exclusion_top100_summary$n_removed_from_full_top100,
    collapse = "; "
  )
)
p_d <- ggplot(
  exclusion_top100_summary,
  aes(factor(component, levels = c("BAS", "CAS", "IAS")), top100_jaccard, fill = component)
) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = sprintf("%.2f", top100_jaccard)),
    vjust = -0.45, size = 2.2, colour = "black"
  ) +
  scale_fill_manual(values = component_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "Top-100 Jaccard after exclusion", caption = exclusion_caption) +
  theme_revision() +
  theme(plot.caption = element_text(size = 5.2, hjust = 0, colour = "#4A4A4A"))

robustness_design <- "
AABB
CCCD
"
robustness_figure <- p_a + p_b + p_c + p_d +
  plot_layout(design = robustness_design, heights = c(1, 1.08)) +
  plot_annotation(tag_levels = "a")

grDevices::cairo_pdf(
  file.path(output_dir, "BAS_CAS_IAS_robustness_top390.pdf"),
  width = 183 / 25.4, height = 145 / 25.4, family = "sans"
)
print(robustness_figure)
dev.off()

fit_plot_data <- rbind(
  data.frame(gene = genes, model = "Condition only", adjusted_R2 = metrics_simple$adjusted_r_squared),
  data.frame(gene = genes, model = "IBI adjusted", adjusted_R2 = metrics_adjusted$adjusted_r_squared),
  data.frame(gene = genes, model = "Full interaction", adjusted_R2 = metrics_full$adjusted_r_squared)
)
fit_plot_data$model <- factor(
  fit_plot_data$model,
  levels = c("Condition only", "IBI adjusted", "Full interaction")
)
p_e <- ggplot(fit_plot_data, aes(model, adjusted_R2, fill = model)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Condition only" = "#9A9A9A", "IBI adjusted" = "#3B6FB6",
    "Full interaction" = "#2A9D8F"
  ), guide = "none") +
  labs(x = NULL, y = "Adjusted R-squared") +
  theme_revision() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

model_pair_labels <- c(
  "Condition-only Grade vs IBI-adjusted Grade" = "Condition only vs\nIBI adjusted",
  "Condition-only Grade vs full-model CAS" = "Condition only vs\nfull-model CAS",
  "IBI-adjusted Grade vs full-model CAS" = "IBI adjusted vs\nfull-model CAS"
)
conventional_summary$comparison_label <- factor(
  model_pair_labels[conventional_summary$comparison],
  levels = unname(model_pair_labels)
)
p_f <- ggplot(conventional_summary, aes(comparison_label, top100_jaccard)) +
  geom_col(width = 0.62, fill = "#D9822B") +
  geom_text(aes(label = sprintf("J = %.3f", top100_jaccard)), vjust = -0.5, size = 2.2) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "Top-100 Jaccard") +
  theme_revision() +
  theme(axis.text.x = element_text(size = 5.6))

nested_plot_data <- nested_summary
nested_plot_data$label <- c("Add IBI", "Add interaction")
p_g <- ggplot(nested_plot_data, aes(label, n_genes_FDR_lt_0_05, fill = label)) +
  geom_col(width = 0.58) +
  geom_text(
    aes(label = n_genes_FDR_lt_0_05),
    vjust = ifelse(nested_plot_data$n_genes_FDR_lt_0_05 > 0, -0.45, -0.8),
    size = 2.3
  ) +
  scale_fill_manual(values = c("Add IBI" = "#3B6FB6", "Add interaction" = "#2A9D8F"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Genes with nested-model FDR < 0.05") +
  theme_revision()

conventional_figure <- (p_e | p_f | p_g) +
  plot_layout(widths = c(1.15, 1.2, 0.9)) +
  plot_annotation(tag_levels = "a")

grDevices::cairo_pdf(
  file.path(output_dir, "BAS_CAS_IAS_conventional_model_comparison.pdf"),
  width = 183 / 25.4, height = 72 / 25.4, family = "sans"
)
print(conventional_figure)
dev.off()

# Optional R-native raster previews for visual QA. This switch is off by
# default, and the temporary PNG files are not part of the formal deliverables.
if (identical(Sys.getenv("BAS_QA_PNG"), "1")) {
  grDevices::png(
    file.path(output_dir, "BAS_CAS_IAS_robustness_top390_QA.png"),
    width = 183 / 25.4, height = 145 / 25.4, units = "in", res = 300,
    family = "sans", type = "cairo"
  )
  print(robustness_figure)
  dev.off()
  grDevices::png(
    file.path(output_dir, "BAS_CAS_IAS_conventional_model_comparison_QA.png"),
    width = 183 / 25.4, height = 72 / 25.4, units = "in", res = 300,
    family = "sans", type = "cairo"
  )
  print(conventional_figure)
  dev.off()
}

provenance <- data.frame(
  item = c(
    "Feature source", "Feature count", "Unique IBI input genes",
    "Gene-pair encoding", "IBI training samples", "IBI training composition",
    "irAE source samples", "irAE analyzed patients", "Excluded samples",
    "Expression transformation", "Expression filter", "Standardization",
    "Clinical variable", "Interaction construction", "Dominance thresholds",
    "Permutation count", "Permutation seed", "Permutation statistic",
    "Formal figure format"
  ),
  value = c(
    basename(feature_file), length(features_top390), length(ibi_input_genes),
    "1 if gene A > gene B; otherwise 0 (including ties); unavailable pairs imputed with the 781-sample development mode",
    nrow(development_predictions), "179 immunodeficiency; 279 healthy; 323 autoimmune",
    ncol(expression_all), ncol(expression), paste(excluded_samples, collapse = "; "),
    "log2(TPM + 1)", "TPM > 0.5 in >20% of analyzed patients; nonzero gene SD",
    "Gene outcome, IBI, Grade, and IBI-by-Grade product each z-scored",
    "irAE Grade 0-4 modeled as a continuous ordinal variable",
    "z-score(z-score(IBI) multiplied by z-score(Grade))",
    paste(thresholds, collapse = "; "), n_permutations, permutation_seed,
    "Mean absolute standardized coefficient among the top 100 genes",
    "PDF only"
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE)

analysis_summary <- data.frame(
  measure = c(
    "Analyzed patients", "Grade 0", "Grade 1", "Grade 2", "Grade 3", "Grade 4",
    "Modeled genes", "Complete top390 pairs", "Unavailable top390 pairs",
    "Spearman IBI vs Grade", "Unique IBI-input genes in modeled genes"
  ),
  value = c(
    nrow(sample_data), as.integer(table(factor(grade_number, levels = 0:4))),
    nrow(component_results), sum(pair_result$available), sum(!pair_result$available),
    cor(ibi, grade_number, method = "spearman"), length(ibi_input_genes_in_analysis)
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_summary, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

save(
  sample_data, component_results, excluded_component_results,
  dominance_summary, threshold_jaccard, conventional_summary,
  nested_summary, exclusion_top100_summary, exclusion_dominance_summary,
  permutation_results, permutation_summary, go_results, go_stability,
  file = file.path(output_dir, "BAS_CAS_IAS_sensitivity_analysis_objects.rda")
)

capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("Completed BAS/CAS/IAS sensitivity analysis in:", output_dir, "\n")
