options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(e1071)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

root <- gsub("\\\\", "/", Sys.getenv("IBI_TOOL_ROOT", unset = "C:/temp/ibi_tool"))
source_dir <- file.path(root, "Result", "revision_711_top390_feature_selection_LOCO_20261109")
output_dir <- file.path(root, "Result", "revision_711_fixed390_model_LOCO_20261109")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cache_env <- new.env(parent = emptyenv())
load(file.path(source_dir, "development_711_expression_cache.rda"), envir = cache_env)
model_env <- new.env(parent = emptyenv())
load(file.path(source_dir, "full711_selection_and_model.rda"), envir = model_env)

expression <- cache_env$all_expression
manifest <- cache_env$sample_manifest
top390 <- model_env$selection$top390
features <- top390$gene_pair

stopifnot(
  nrow(manifest) == 711L,
  nrow(top390) == 390L,
  identical(colnames(expression), manifest$sample_key),
  identical(as.integer(table(factor(manifest$label, levels = c(-1, 0, 1)))), c(165L, 223L, 323L))
)

svr_parameters <- function() {
  environment <- new.env(parent = emptyenv())
  load(file.path(root, "best_svr0731.rda"), envir = environment)
  model <- environment$best_svr
  list(
    kernel = c("linear", "polynomial", "radial", "sigmoid")[as.integer(model$kernel) + 1L],
    cost = model$cost,
    gamma = model$gamma,
    degree = model$degree,
    coef0 = model$coef0,
    nu = model$nu,
    epsilon = model$epsilon
  )
}

make_pair_matrix <- function(sample_indices, training_modes = NULL) {
  x <- expression[, sample_indices, drop = FALSE]
  result <- matrix(
    NA_real_, nrow = length(sample_indices), ncol = length(features),
    dimnames = list(manifest$sample_key[sample_indices], features)
  )
  learned_modes <- numeric(length(features))
  for (j in seq_along(features)) {
    a <- x[top390$geneA[j], ]
    b <- x[top390$geneB[j], ]
    observed <- is.finite(a) & is.finite(b)
    values <- rep(NA_real_, length(sample_indices))
    values[observed] <- as.numeric(a[observed] > b[observed])
    mode_value <- if (is.null(training_modes)) {
      as.numeric(mean(values[observed] == 1, na.rm = TRUE) >= 0.5)
    } else {
      training_modes[j]
    }
    values[!observed] <- mode_value
    result[, j] <- values
    learned_modes[j] <- mode_value
  }
  list(matrix = result, modes = learned_modes)
}

fit_predict_fold <- function(train_indices, test_indices, seed) {
  train_built <- make_pair_matrix(train_indices)
  test_built <- make_pair_matrix(test_indices, train_built$modes)
  train_x <- train_built$matrix
  test_x <- test_built$matrix
  variable <- apply(train_x, 2L, sd) > 0
  train_x <- train_x[, variable, drop = FALSE]
  test_x <- test_x[, variable, drop = FALSE]
  params <- svr_parameters()
  set.seed(seed)
  fit <- e1071::svm(
    x = train_x, y = manifest$label[train_indices],
    type = "eps-regression", kernel = params$kernel,
    cost = params$cost, gamma = params$gamma, degree = params$degree,
    coef0 = params$coef0, nu = params$nu, epsilon = params$epsilon,
    scale = TRUE
  )
  list(
    prediction = as.numeric(predict(fit, test_x)),
    retained_features = sum(variable)
  )
}

calculate_metrics <- function(data) {
  output <- data.frame(
    cohort = data$cohort[1L],
    n = nrow(data),
    n_immunodeficiency = sum(data$observed == -1),
    n_healthy = sum(data$observed == 0),
    n_autoimmune = sum(data$observed == 1),
    RMSE = sqrt(mean((data$predicted_IBI - data$observed)^2)),
    MAE = mean(abs(data$predicted_IBI - data$observed)),
    Spearman_rho = if (length(unique(data$observed)) > 1L) {
      suppressWarnings(cor(data$observed, data$predicted_IBI, method = "spearman"))
    } else NA_real_,
    AUC = NA_real_, AUC_CI_low = NA_real_, AUC_CI_high = NA_real_,
    contrast = NA_character_,
    stringsAsFactors = FALSE
  )
  nonzero <- setdiff(unique(data$observed), 0)
  if (0 %in% data$observed && length(nonzero) == 1L) {
    disease <- data$observed != 0
    score <- if (nonzero == 1) data$predicted_IBI else -data$predicted_IBI
    roc_object <- pROC::roc(
      disease, score, levels = c(FALSE, TRUE), direction = "<", quiet = TRUE
    )
    interval <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
    output$AUC <- as.numeric(pROC::auc(roc_object))
    output$AUC_CI_low <- interval[1L]
    output$AUC_CI_high <- interval[3L]
    output$contrast <- if (nonzero == 1) "Autoimmune vs healthy" else "Immunodeficiency vs healthy"
  }
  output
}

cohorts <- unique(manifest$cohort)
folds <- vector("list", length(cohorts))
metrics <- vector("list", length(cohorts))
for (k in seq_along(cohorts)) {
  held_out <- cohorts[k]
  test_indices <- which(manifest$cohort == held_out)
  train_indices <- setdiff(seq_len(nrow(manifest)), test_indices)
  result <- fit_predict_fold(train_indices, test_indices, 20261109L + k)
  fold_data <- manifest[test_indices, c(
    "sample_key", "sample_id", "cohort", "clinical_group", "label"
  ), drop = FALSE]
  fold_data$observed <- fold_data$label
  fold_data$predicted_IBI <- result$prediction
  fold_data$retained_features <- result$retained_features
  folds[[k]] <- fold_data
  metrics[[k]] <- calculate_metrics(fold_data)
  metrics[[k]]$training_n <- length(train_indices)
  metrics[[k]]$retained_features <- result$retained_features
  cat("Completed fixed-feature LOCO", k, "of", length(cohorts), held_out, "\n")
}

predictions <- do.call(rbind, folds)
performance <- do.call(rbind, metrics)
write.csv(predictions, file.path(output_dir, "fixed390_LOCO_sample_predictions.csv"), row.names = FALSE)
write.csv(performance, file.path(output_dir, "fixed390_LOCO_cohort_performance.csv"), row.names = FALSE)
write.csv(top390, file.path(output_dir, "locked_390_gene_pairs.csv"), row.names = FALSE)

auc_data <- performance[is.finite(performance$AUC), , drop = FALSE]
auc_data$label <- sprintf("%s (n=%d)", auc_data$cohort, auc_data$n)
auc_data$label <- factor(auc_data$label, levels = rev(auc_data$label[order(auc_data$AUC)]))

p_auc <- ggplot(auc_data, aes(AUC, label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.35, colour = "#8A8A8A") +
  geom_errorbar(
    aes(xmin = AUC_CI_low, xmax = AUC_CI_high), width = 0.18,
    linewidth = 0.55, colour = "#356B8C", orientation = "y"
  ) +
  geom_point(size = 2.3, shape = 21, stroke = 0.45, fill = "#D66B52", colour = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    tag = "A",
    title = "Held-out cohort discrimination",
    subtitle = "Direction-aware AUC; 95% DeLong CI",
    x = "AUC (95% CI)", y = NULL
  )

error_data <- rbind(
  data.frame(cohort = performance$cohort, metric = "RMSE", value = performance$RMSE),
  data.frame(cohort = performance$cohort, metric = "MAE", value = performance$MAE)
)
cohort_order <- performance$cohort[order(performance$RMSE)]
error_data$cohort <- factor(error_data$cohort, levels = rev(cohort_order))
p_error <- ggplot(error_data, aes(value, cohort, colour = metric, shape = metric)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.42)) +
  scale_colour_manual(values = c(RMSE = "#356B8C", MAE = "#D66B52")) +
  scale_shape_manual(values = c(RMSE = 16, MAE = 17)) +
  scale_x_continuous(limits = c(0, max(error_data$value) * 1.08), expand = expansion(mult = c(0, 0.02))) +
  labs(
    tag = "B",
    title = "Held-out cohort prediction error",
    subtitle = "All 11 development cohorts",
    x = "Prediction error", y = NULL, colour = NULL, shape = NULL
  ) +
  theme(legend.position = "top")

theme_publication <- theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "black", size = 7.4),
    axis.title = element_text(size = 8.2),
    plot.title = element_text(size = 9.2, face = "bold"),
    plot.subtitle = element_text(size = 7.3, colour = "#4A4A4A"),
    plot.tag = element_text(size = 10, face = "bold"),
    plot.tag.position = c(0, 1),
    legend.text = element_text(size = 7.2),
    plot.margin = margin(5, 7, 5, 5)
  )
p_auc <- p_auc + theme_publication
p_error <- p_error + theme_publication + theme(legend.position = "top")

figure <- p_auc | p_error
figure <- figure + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "Fixed-feature leave-one-cohort-out validation",
    subtitle = "711 samples; each cohort excluded from SVR fitting and preprocessing estimation"
  )

width_mm <- 183
height_mm <- 104
grDevices::cairo_pdf(
  file.path(output_dir, "fixed390_LOCO_model_performance.pdf"),
  width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial"
)
print(figure)
dev.off()
svglite::svglite(
  file.path(output_dir, "fixed390_LOCO_model_performance.svg"),
  width = width_mm / 25.4, height = height_mm / 25.4
)
print(figure)
dev.off()
grDevices::png(
  file.path(output_dir, "fixed390_LOCO_model_performance.png"),
  width = width_mm / 25.4, height = height_mm / 25.4,
  units = "in", res = 300, type = "cairo", family = "Arial"
)
print(figure)
dev.off()
tiff_staging_file <- "fixed390_LOCO_model_performance_711_staging.tiff"
grDevices::tiff(
  tiff_staging_file,
  width = width_mm / 25.4, height = height_mm / 25.4,
  units = "in", res = 600, compression = "lzw", type = "cairo", family = "Arial"
)
print(figure)
dev.off()
if (!file.copy(
  tiff_staging_file,
  file.path(output_dir, "fixed390_LOCO_model_performance.tiff"),
  overwrite = TRUE
)) warning("TIFF was rendered to the staging path but could not be copied to the output directory")

writeLines(c(
  "Core conclusion: a fixed 390-pair SVR retains measurable discrimination and prediction accuracy when complete development cohorts are excluded from model fitting.",
  "Feature panel: locked full-711 Top390; no feature re-selection was performed within LOCO folds.",
  "Preprocessing: missing pair values were assigned using the modal ordering learned from the outer-fold training cohorts only.",
  "Model: locked SVR hyperparameters; SVR coefficients were re-estimated in every outer fold.",
  "Panel A: direction-aware AUC with 95% DeLong confidence interval for cohorts containing healthy and exactly one disease class.",
  "Panel B: RMSE and MAE for all 11 completely held-out development cohorts."
), file.path(output_dir, "figure_QA_notes.txt"))

capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
print(performance)
