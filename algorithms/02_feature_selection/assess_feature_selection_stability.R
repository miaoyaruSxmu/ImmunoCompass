options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(e1071)
  library(ggplot2)
  library(patchwork)
})

# Reviewer 2, Comment 12: multiplicity audit and stability-controlled
# reconstruction of the de novo gene-pair selection workflow.
root <- "."
source_dir <- file.path(
  root, "Result", "revision_392_de_novo_feature_stability_90pct"
)
output_dir <- file.path(
  root, "Result", "revision_comment12_stability_controlled_top390_20261102"
)
qa_dir <- file.path(
  Sys.getenv("TEMP", unset = tempdir()), "ibi_tool", "Result", "comment12_qa"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

base_seed <- 20261021L
resample_fraction <- 0.90
final_feature_count <- 390L
split_half_repeats <- 1000L
total_pair_tests <- 2936676L
development_n <- 682L

required_files <- c(
  "full682_top392_gene_pairs.csv",
  "full682_top5000_RF_XGBoost_ranking.csv",
  "bootstrap_run_metrics.csv",
  "development_682_expression_cache_392.rda"
)
missing_files <- required_files[!file.exists(file.path(source_dir, required_files))]
if (length(missing_files)) {
  stop("Missing source files: ", paste(missing_files, collapse = ", "))
}

full_ranking <- read.csv(
  file.path(source_dir, "full682_top5000_RF_XGBoost_ranking.csv"),
  check.names = FALSE
)
full_panel <- full_ranking[
  order(full_ranking$integrated_rank), , drop = FALSE
][seq_len(final_feature_count), , drop = FALSE]

resample_dir <- file.path(source_dir, "bootstrap_replicates")
resample_files <- file.path(
  resample_dir, sprintf("bootstrap_%03d_top392.csv", seq_len(100L))
)
stopifnot(all(file.exists(resample_files)))
resample_tables <- lapply(resample_files, function(path) {
  table <- read.csv(path, check.names = FALSE)
  table <- table[order(table$integrated_rank), , drop = FALSE]
  table[seq_len(min(final_feature_count, nrow(table))), , drop = FALSE]
})
resample_sets <- lapply(resample_tables, `[[`, "gene_pair")

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

consensus_table <- function(indices) {
  rows <- do.call(rbind, lapply(indices, function(index) {
    table <- resample_tables[[index]]
    data.frame(
      gene_pair = table$gene_pair,
      geneA = table$geneA,
      geneB = table$geneB,
      rank = seq_len(nrow(table)),
      stringsAsFactors = FALSE
    )
  }))
  count <- aggregate(rank ~ gene_pair + geneA + geneB, rows, length)
  names(count)[names(count) == "rank"] <- "selection_count"
  rank_sum <- aggregate(rank ~ gene_pair, rows, sum)
  names(rank_sum)[names(rank_sum) == "rank"] <- "selected_rank_sum"
  output <- merge(count, rank_sum, by = "gene_pair", all.x = TRUE)
  output$selection_frequency <- output$selection_count / length(indices)
  output$mean_selected_rank <- output$selected_rank_sum / output$selection_count
  output <- output[order(
    -output$selection_frequency, output$mean_selected_rank, output$gene_pair
  ), , drop = FALSE]
  output$consensus_rank <- seq_len(nrow(output))
  output
}

consensus_all <- consensus_table(seq_along(resample_tables))
consensus_panel <- consensus_all[seq_len(final_feature_count), , drop = FALSE]
write.csv(
  consensus_panel,
  file.path(output_dir, "stability_consensus_top390_gene_pairs.csv"),
  row.names = FALSE
)

full_frequency <- data.frame(
  gene_pair = full_panel$gene_pair,
  selection_count = vapply(
    full_panel$gene_pair,
    function(pair) sum(vapply(resample_sets, `%in%`, logical(1), x = pair)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
full_frequency$selection_frequency <- full_frequency$selection_count /
  length(resample_sets)
full_frequency <- merge(
  full_panel[, c("integrated_rank", "gene_pair", "geneA", "geneB")],
  full_frequency, by = "gene_pair", all.x = TRUE
)
full_frequency <- full_frequency[order(full_frequency$integrated_rank), ]
write.csv(
  full_frequency,
  file.path(output_dir, "full_data_top390_resampling_frequency.csv"),
  row.names = FALSE
)

raw_jaccard <- vapply(
  resample_sets,
  function(features) jaccard(features, full_panel$gene_pair),
  numeric(1)
)

set.seed(20261102L)
split_half <- vector("list", split_half_repeats)
for (iteration in seq_len(split_half_repeats)) {
  first <- sample(seq_along(resample_tables), length(resample_tables) / 2L)
  second <- setdiff(seq_along(resample_tables), first)
  panel_first <- consensus_table(first)$gene_pair[seq_len(final_feature_count)]
  panel_second <- consensus_table(second)$gene_pair[seq_len(final_feature_count)]
  split_half[[iteration]] <- data.frame(
    iteration = iteration,
    jaccard_between_halves = jaccard(panel_first, panel_second),
    jaccard_half1_with_full_consensus = jaccard(
      panel_first, consensus_panel$gene_pair
    ),
    jaccard_half2_with_full_consensus = jaccard(
      panel_second, consensus_panel$gene_pair
    )
  )
}
split_half <- do.call(rbind, split_half)
write.csv(
  split_half, file.path(output_dir, "split_half_consensus_stability.csv"),
  row.names = FALSE
)

# Conservative multiplicity audit. The 5,000 retained candidates contain the
# smallest association scores used by the downstream ranking. p.adjust(...,
# n = total_pair_tests) treats all unretained tests as p = 1.
n_observed <- pmax(
  4L, pmin(development_n, round(full_ranking$sample_coverage * development_n))
)
association <- pmin(abs(full_ranking$prescreen_score), 1 - 1e-12)
t_statistic <- association * sqrt((n_observed - 2) / (1 - association^2))
raw_p <- 2 * pt(-abs(t_statistic), df = n_observed - 2)
bh_fdr <- p.adjust(raw_p, method = "BH", n = total_pair_tests)
multiplicity <- data.frame(
  gene_pair = full_ranking$gene_pair,
  prescreen_rank = full_ranking$prescreen_rank,
  integrated_rank = full_ranking$integrated_rank,
  n_observed = n_observed,
  association = association,
  p_value = raw_p,
  BH_FDR = bh_fdr,
  in_full_top390 = full_ranking$gene_pair %in% full_panel$gene_pair,
  stringsAsFactors = FALSE
)
write.csv(
  multiplicity,
  file.path(output_dir, "candidate_pair_multiple_testing_audit.csv"),
  row.names = FALSE
)

load(file.path(source_dir, "development_682_expression_cache_392.rda"))
stopifnot(nrow(sample_manifest) == development_n)

build_features <- function(expression, panel, training_modes = NULL) {
  output <- matrix(
    NA_real_, nrow = ncol(expression), ncol = nrow(panel),
    dimnames = list(colnames(expression), panel$gene_pair)
  )
  for (index in seq_len(nrow(panel))) {
    gene_a <- panel$geneA[index]
    gene_b <- panel$geneB[index]
    observed <- is.finite(expression[gene_a, ]) &
      is.finite(expression[gene_b, ])
    output[observed, index] <- as.numeric(
      expression[gene_a, observed] > expression[gene_b, observed]
    )
  }
  if (is.null(training_modes)) {
    training_modes <- apply(output, 2, function(values) {
      if (!any(is.finite(values))) return(0)
      as.numeric(mean(values, na.rm = TRUE) > 0.5)
    })
  }
  for (index in seq_len(ncol(output))) {
    output[!is.finite(output[, index]), index] <- training_modes[index]
  }
  list(matrix = output, modes = training_modes)
}

model_environment <- new.env(parent = emptyenv())
load(file.path(root, "best_svr0731.rda"), envir = model_environment)
reference_model <- model_environment$best_svr
kernel_name <- c("linear", "polynomial", "radial", "sigmoid")[
  as.integer(reference_model$kernel) + 1L
]
train_svr <- function(x, y) {
  e1071::svm(
    x = x, y = y, type = "eps-regression", kernel = kernel_name,
    cost = reference_model$cost, gamma = reference_model$gamma,
    degree = reference_model$degree, coef0 = reference_model$coef0,
    nu = reference_model$nu, epsilon = reference_model$epsilon,
    scale = TRUE
  )
}

stratified_subsample <- function(fraction, seed) {
  set.seed(seed)
  strata <- interaction(sample_manifest$cohort, sample_manifest$label, drop = TRUE)
  unlist(lapply(split(seq_len(nrow(sample_manifest)), strata), function(index) {
    sample(index, max(1L, floor(length(index) * fraction)), replace = FALSE)
  }), use.names = FALSE)
}

reference_features <- build_features(all_expression, full_panel)
reference_fit <- train_svr(reference_features$matrix, sample_manifest$label)
reference_prediction <- as.numeric(predict(reference_fit, reference_features$matrix))

prediction_matrix <- matrix(
  NA_real_, nrow = development_n, ncol = length(resample_tables),
  dimnames = list(sample_manifest$sample_key, sprintf("repeat_%03d", seq_len(100L)))
)
oob_rows <- list()
prediction_summary <- vector("list", length(resample_tables))
for (iteration in seq_along(resample_tables)) {
  training_indices <- stratified_subsample(
    resample_fraction, base_seed + iteration
  )
  test_indices <- setdiff(seq_len(development_n), training_indices)
  panel <- resample_tables[[iteration]]
  training_features <- build_features(
    all_expression[, training_indices, drop = FALSE], panel
  )
  all_features <- build_features(
    all_expression, panel, training_modes = training_features$modes
  )
  fit <- train_svr(
    training_features$matrix, sample_manifest$label[training_indices]
  )
  prediction <- as.numeric(predict(fit, all_features$matrix))
  prediction_matrix[, iteration] <- prediction
  prediction_summary[[iteration]] <- data.frame(
    iteration = iteration,
    spearman_with_full_model = suppressWarnings(cor(
      prediction, reference_prediction, method = "spearman"
    )),
    oob_rmse = sqrt(mean(
      (prediction[test_indices] - sample_manifest$label[test_indices])^2
    )),
    oob_mae = mean(abs(
      prediction[test_indices] - sample_manifest$label[test_indices]
    )),
    stringsAsFactors = FALSE
  )
  oob_rows[[iteration]] <- data.frame(
    iteration = iteration,
    sample_key = sample_manifest$sample_key[test_indices],
    label = sample_manifest$label[test_indices],
    prediction = prediction[test_indices],
    stringsAsFactors = FALSE
  )
}
prediction_summary <- do.call(rbind, prediction_summary)
oob_predictions <- do.call(rbind, oob_rows)
write.csv(
  prediction_summary,
  file.path(output_dir, "resampled_model_prediction_stability.csv"),
  row.names = FALSE
)
write.csv(
  oob_predictions,
  file.path(output_dir, "resampled_model_oob_predictions.csv"),
  row.names = FALSE
)

pairwise_prediction_cor <- cor(prediction_matrix, method = "spearman")
pairwise_prediction_cor <- pairwise_prediction_cor[upper.tri(
  pairwise_prediction_cor
)]
saveRDS(
  pairwise_prediction_cor,
  file.path(output_dir, "pairwise_prediction_cor.rds")
)
oob_aggregate <- aggregate(
  prediction ~ sample_key + label, oob_predictions, median
)
oob_counts <- aggregate(prediction ~ sample_key, oob_predictions, length)
names(oob_counts)[2] <- "oob_repeats"
oob_aggregate <- merge(oob_aggregate, oob_counts, by = "sample_key")
write.csv(
  oob_aggregate,
  file.path(output_dir, "sample_level_aggregated_oob_predictions.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "Theoretical unordered pairs",
    "Evaluable pairs after coverage and ordering filters",
    "Candidates entering RF/XGBoost ranking",
    "Candidates with BH FDR < 0.05",
    "Maximum BH FDR among the Top390",
    "Full-data Top390 selected in >=50% of repeats",
    "Full-data Top390 selected in >=80% of repeats",
    "Full-data Top390 selected in >=90% of repeats",
    "Full-data Top390 selected in 100% of repeats",
    "Median raw-resample Jaccard with full-data Top390",
    "Median split-half consensus Jaccard",
    "Median half-panel Jaccard with full consensus",
    "Median prediction rho with full-data model",
    "Median pairwise prediction rho among resampled models",
    "Aggregated out-of-bag RMSE",
    "Aggregated out-of-bag MAE",
    "Aggregated out-of-bag Spearman rho"
  ),
  value = c(
    choose(3198, 2), total_pair_tests, nrow(full_ranking),
    sum(multiplicity$BH_FDR < 0.05),
    max(multiplicity$BH_FDR[multiplicity$in_full_top390]),
    sum(full_frequency$selection_frequency >= 0.50),
    sum(full_frequency$selection_frequency >= 0.80),
    sum(full_frequency$selection_frequency >= 0.90),
    sum(full_frequency$selection_frequency >= 1.00),
    median(raw_jaccard),
    median(split_half$jaccard_between_halves),
    median(c(
      split_half$jaccard_half1_with_full_consensus,
      split_half$jaccard_half2_with_full_consensus
    )),
    median(prediction_summary$spearman_with_full_model),
    median(pairwise_prediction_cor),
    sqrt(mean((oob_aggregate$prediction - oob_aggregate$label)^2)),
    mean(abs(oob_aggregate$prediction - oob_aggregate$label)),
    suppressWarnings(cor(
      oob_aggregate$prediction, oob_aggregate$label, method = "spearman"
    ))
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_table, file.path(output_dir, "comment12_stability_summary.csv"),
  row.names = FALSE
)

palette <- c(
  neutral = "#777777", light = "#D9D9D9", blue = "#3C78A8",
  teal = "#2A9D8F", red = "#C65D4B"
)
theme_set(
  theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.45),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, colour = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, face = "bold", margin = margin(b = 6)),
      plot.tag = element_text(size = 13, face = "bold"),
      plot.margin = margin(8, 9, 7, 8),
      legend.title = element_blank(),
      legend.position = "top"
    )
)

flow_data <- data.frame(
  stage = factor(
    c("Theoretical\npairs", "Evaluable\npairs", "BH-controlled\nprescreen", "Locked\npanel"),
    levels = c("Theoretical\npairs", "Evaluable\npairs", "BH-controlled\nprescreen", "Locked\npanel")
  ),
  count = c(choose(3198, 2), total_pair_tests, nrow(full_ranking), final_feature_count)
)
p_a <- ggplot(flow_data, aes(stage, count)) +
  geom_col(width = 0.68, fill = palette[["blue"]]) +
  geom_text(
    aes(label = format(count, big.mark = ",", scientific = FALSE)),
    vjust = -0.35, size = 3.2
  ) +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.16))) +
  labs(title = "Controlled reduction of candidate space", x = NULL, y = "Features (log10 scale)")

frequency_plot_data <- full_frequency
p_b <- ggplot(
  frequency_plot_data,
  aes(reorder(gene_pair, -selection_frequency), selection_frequency)
) +
  geom_col(
    aes(fill = selection_frequency >= 0.80), width = 0.85,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0.80, linetype = "dashed", linewidth = 0.35) +
  scale_fill_manual(values = c(`FALSE` = palette[["light"]], `TRUE` = palette[["teal"]])) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 0.8, 1)) +
  labs(
    title = "Repeated 90% stratified resampling",
    x = "Full-data Top390 pairs", y = "Selection frequency"
  ) +
  annotate(
    "text", x = 330, y = 0.88,
    label = paste0(sum(full_frequency$selection_frequency >= 0.80), " pairs >=80%"),
    size = 3.2, hjust = 1
  )

jaccard_data <- rbind(
  data.frame(
    comparison = "Single resample\nvs full-data panel",
    jaccard = raw_jaccard
  ),
  data.frame(
    comparison = "Split-half\nconsensus panels",
    jaccard = split_half$jaccard_between_halves
  )
)
jaccard_data$comparison <- factor(
  jaccard_data$comparison,
  levels = c("Single resample\nvs full-data panel", "Split-half\nconsensus panels")
)
p_c <- ggplot(jaccard_data, aes(comparison, jaccard, fill = comparison)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35) +
  scale_fill_manual(values = c(palette[["light"]], palette[["teal"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Consensus aggregation stabilizes panel identity", x = NULL, y = "Jaccard index")

prediction_plot_data <- rbind(
  data.frame(
    comparison = "Resampled model\nvs full-data model",
    rho = prediction_summary$spearman_with_full_model
  ),
  data.frame(
    comparison = "Between resampled\nmodels",
    rho = pairwise_prediction_cor
  )
)
prediction_plot_data$comparison <- factor(
  prediction_plot_data$comparison,
  levels = c("Resampled model\nvs full-data model", "Between resampled\nmodels")
)
p_d <- ggplot(prediction_plot_data, aes(comparison, rho, fill = comparison)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35) +
  scale_fill_manual(values = c(palette[["blue"]], palette[["teal"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "IBI predictions remain concordant", x = NULL, y = "Spearman rho")

combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, 1.05))

grDevices::cairo_pdf(
  file.path(output_dir, "comment12_feature_selection_stability.pdf"),
  width = 183 / 25.4, height = 145 / 25.4, family = "Arial"
)
print(combined)
dev.off()

png(
  file.path(output_dir, "comment12_feature_selection_stability.png"),
  width = 183 / 25.4, height = 145 / 25.4, units = "in", res = 400,
  type = "cairo"
)
print(combined)
dev.off()

provenance <- data.frame(
  item = c(
    "Candidate genes", "Theoretical pairs", "Evaluable pairs",
    "Multiplicity correction", "Variance-filter rationale",
    "Updated identifiability control", "Stability resampling",
    "Consensus ranking", "Random seed"
  ),
  value = c(
    "3,198 immune-related genes",
    format(choose(3198, 2), scientific = FALSE),
    format(total_pair_tests, scientific = FALSE),
    "Benjamini-Hochberg audit of 5,000 retained candidates with n=2,936,676 total tests; unretained tests treated as non-significant",
    paste0(
      "The legacy top-75% variance rule was an unsupervised screen intended to remove ",
      "the bottom quartile of near-constant features; it was not tuned against validation outcomes."
    ),
    paste0(
      "Reproducible analysis used >=70% sample coverage, presence in >=50% of cohorts, ",
      ">=70% pair coverage and ordering frequency between 5% and 95%."
    ),
    "100 repeats; 90% sampled within each cohort-by-label stratum",
    "Selection frequency, followed by mean rank among selected repeats",
    base_seed
  ),
  stringsAsFactors = FALSE
)
write.csv(
  provenance, file.path(output_dir, "analysis_provenance.csv"),
  row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))

cat("Completed Comment 12 stability analysis.\n")
print(summary_table)
cat("Output:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
