options(stringsAsFactors = FALSE)

root <- "C:/temp/ibi_tool/Result"
output_dir <- file.path(root, "revision_encoding_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE51404", "GSE44969")
summary_rows <- list()
sample_rows <- list()

for (dataset in datasets) {
  zero_one_dir <- file.path(root, paste0("revision_", dataset, "_validation"))
  minus_one_dir <- file.path(root, paste0("revision_", dataset, "_validation_pm1"))

  zero_one <- read.csv(file.path(zero_one_dir, paste0(dataset, "_sample_level_IBI.csv")))
  minus_one <- read.csv(file.path(minus_one_dir, paste0(dataset, "_sample_level_IBI.csv")))
  predictions <- merge(
    zero_one[, c("sample_id", "IBI")],
    minus_one[, c("sample_id", "IBI")],
    by = "sample_id", suffixes = c("_01", "_pm1"), sort = FALSE
  )
  predictions$dataset <- dataset
  predictions$difference_pm1_minus_01 <- predictions$IBI_pm1 - predictions$IBI_01
  sample_rows[[dataset]] <- predictions[, c(
    "dataset", "sample_id", "IBI_01", "IBI_pm1", "difference_pm1_minus_01"
  )]

  group_zero_one <- read.csv(file.path(zero_one_dir, paste0(dataset, "_IBI_group_tests.csv")))
  group_minus_one <- read.csv(file.path(minus_one_dir, paste0(dataset, "_IBI_group_tests.csv")))
  signature_zero_one <- read.csv(
    file.path(zero_one_dir, paste0(dataset, "_IBI_signature_correlations.csv"))
  )
  signature_minus_one <- read.csv(
    file.path(minus_one_dir, paste0(dataset, "_IBI_signature_correlations.csv"))
  )

  summary_rows[[dataset]] <- data.frame(
    dataset = dataset,
    n = nrow(predictions),
    pearson_prediction_correlation = cor(predictions$IBI_01, predictions$IBI_pm1),
    spearman_prediction_correlation = cor(
      predictions$IBI_01, predictions$IBI_pm1, method = "spearman"
    ),
    mean_absolute_prediction_difference = mean(abs(predictions$difference_pm1_minus_01)),
    maximum_absolute_prediction_difference = max(abs(predictions$difference_pm1_minus_01)),
    maximum_group_test_p_difference = max(abs(
      group_zero_one$p_value - group_minus_one$p_value
    )),
    maximum_signature_rho_difference = max(abs(
      signature_zero_one$rho - signature_minus_one$rho
    )),
    maximum_signature_p_difference = max(abs(
      signature_zero_one$p_value - signature_minus_one$p_value
    ))
  )
}

write.csv(
  do.call(rbind, summary_rows),
  file.path(output_dir, "IBI_encoding_sensitivity_summary.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, sample_rows),
  file.path(output_dir, "IBI_encoding_sensitivity_sample_level.csv"),
  row.names = FALSE
)

cat("Encoding sensitivity comparison completed:", output_dir, "\n")
