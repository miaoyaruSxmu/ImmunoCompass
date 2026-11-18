options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(e1071)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(
    paste(
      "Usage: Rscript train_final_svr_from_cohort_table.R",
      "<fixed_781_rda> <GSE122459_xlsx> <TableS2_xlsx>",
      "<output_dir> <package_data_dir> [package_extdata_dir]"
    ),
    call. = FALSE
  )
}

feature_file <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
gse122459_file <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
table_s2_file <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[4L]]
package_data_dir <- args[[5L]]
package_extdata_dir <- if (length(args) >= 6L) args[[6L]] else NA_character_

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(package_data_dir, recursive = TRUE, showWarnings = FALSE)
if (!is.na(package_extdata_dir)) {
  dir.create(package_extdata_dir, recursive = TRUE, showWarnings = FALSE)
}

development_cohorts <- c(
  "GTEx", "GSE264125", "GSE222408", "GSE175839", "GSE181228",
  "GSE199911", "GSE140713", "GSE297198", "GSE232381",
  "GSE59216_batch1", "GSE51405", "GSE122459"
)

expected_table_counts <- data.frame(
  cohort = development_cohorts,
  n = c(70L, 21L, 61L, 10L, 277L, 75L, 57L, 11L, 25L, 91L, 83L, 26L),
  stringsAsFactors = FALSE
)

normalize_cohort <- function(x) {
  x <- trimws(as.character(x))
  x[x == "GSE59216 (batch 1)"] <- "GSE59216_batch1"
  x
}

parse_analysis_n <- function(x) {
  value <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", trimws(x))))
  if (anyNA(value)) stop("Could not parse one or more Table S2 analysis counts")
  value
}

table_s2 <- as.data.frame(read_excel(table_s2_file), check.names = FALSE)
required_columns <- c(
  "Dataset/accession", "Analysis role", "Analysis n", "Training label"
)
if (!all(required_columns %in% names(table_s2))) {
  stop(
    "Table S2 is missing columns: ",
    paste(setdiff(required_columns, names(table_s2)), collapse = "; ")
  )
}

table_training <- table_s2[
  trimws(table_s2[["Analysis role"]]) == "Model development", , drop = FALSE
]
table_training$cohort <- normalize_cohort(table_training[["Dataset/accession"]])
table_training$n <- parse_analysis_n(table_training[["Analysis n"]])
table_training$label <- suppressWarnings(
  as.integer(trimws(as.character(table_training[["Training label"]])))
)
if (anyNA(table_training$label)) stop("Non-numeric training label in Table S2")

table_cohort_counts <- aggregate(n ~ cohort, table_training, sum)
table_cohort_counts <- table_cohort_counts[
  match(development_cohorts, table_cohort_counts$cohort), , drop = FALSE
]
rownames(table_cohort_counts) <- NULL
if (!identical(table_cohort_counts, expected_table_counts)) {
  stop("Table S2 model-development cohorts or sample counts do not match the expected 807-sample allocation")
}
if (sum(table_training$n) != 807L) stop("Table S2 development total is not 807")

table_label_counts <- aggregate(n ~ label, table_training, sum)
table_label_counts <- table_label_counts[order(table_label_counts$label), ]
expected_label_counts <- data.frame(
  label = c(-1L, 0L, 1L), n = c(179L, 285L, 343L)
)
rownames(table_label_counts) <- NULL
if (!identical(table_label_counts, expected_label_counts)) {
  stop("Table S2 development-label counts are not 179/285/343")
}

feature_environment <- new.env(parent = emptyenv())
load(feature_file, envir = feature_environment)
required_objects <- c(
  "features_top390", "top390", "top390_feature_matrix",
  "sample_manifest", "development_cohorts"
)
if (!all(required_objects %in% ls(feature_environment))) {
  stop(
    "The fixed-feature object is missing: ",
    paste(setdiff(required_objects, ls(feature_environment)), collapse = "; ")
  )
}

gene_pairs <- as.character(feature_environment$features_top390)
top390 <- feature_environment$top390
base_matrix <- as.matrix(feature_environment$top390_feature_matrix)
base_manifest <- as.data.frame(feature_environment$sample_manifest)
base_development_cohorts <- as.character(feature_environment$development_cohorts)

base_rows <- base_manifest$cohort %in% base_development_cohorts
base_manifest <- base_manifest[base_rows, , drop = FALSE]
base_matrix <- base_matrix[base_manifest$sample_key, gene_pairs, drop = FALSE]
storage.mode(base_matrix) <- "double"
if (any(!is.finite(base_matrix))) {
  stop("The 781-sample fixed-feature matrix contains missing values")
}
if (!all(base_matrix %in% c(0, 1))) {
  stop("The 781-sample fixed-feature matrix is not encoded as 0/1")
}
stopifnot(
  nrow(base_matrix) == 781L,
  ncol(base_matrix) == 390L,
  identical(colnames(base_matrix), gene_pairs)
)

gse_table <- as.data.frame(
  read_excel(gse122459_file, sheet = "merged normalized"),
  check.names = FALSE
)
if (!"name" %in% names(gse_table)) stop("GSE122459 workbook lacks the name column")
gse_samples <- grep("^(SLE|Healthy)[0-9]+$", names(gse_table), value = TRUE)
if (
  length(gse_samples) != 26L ||
  sum(grepl("^SLE", gse_samples)) != 20L ||
  sum(grepl("^Healthy", gse_samples)) != 6L
) {
  stop("GSE122459 must contain 20 SLE and 6 healthy samples")
}

gse_genes <- toupper(trimws(as.character(gse_table[["name"]])))
gse_expression <- as.matrix(gse_table[, gse_samples, drop = FALSE])
storage.mode(gse_expression) <- "double"
valid_gene <- !is.na(gse_genes) & nzchar(gse_genes)
gse_genes <- gse_genes[valid_gene]
gse_expression <- gse_expression[valid_gene, , drop = FALSE]

gene_sums <- rowsum(gse_expression, group = gse_genes, reorder = FALSE)
gene_counts <- table(gse_genes)
gse_expression <- sweep(
  gene_sums, 1L, as.numeric(gene_counts[rownames(gene_sums)]), "/"
)

base_modes <- as.integer(colMeans(base_matrix == 1L) >= 0.5)
names(base_modes) <- gene_pairs
gse_matrix <- matrix(
  rep(base_modes, each = length(gse_samples)),
  nrow = length(gse_samples), ncol = length(gene_pairs),
  dimnames = list(paste("GSE122459", gse_samples, sep = "::"), gene_pairs)
)
pair_available <- logical(length(gene_pairs))
for (j in seq_along(gene_pairs)) {
  genes <- strsplit(gene_pairs[[j]], "_vs_", fixed = TRUE)[[1L]]
  pair_available[[j]] <- all(genes %in% rownames(gse_expression))
  if (pair_available[[j]]) {
    a <- gse_expression[genes[[1L]], ]
    b <- gse_expression[genes[[2L]], ]
    observed <- is.finite(a) & is.finite(b)
    gse_matrix[observed, j] <- as.numeric(a[observed] > b[observed])
  }
}

gse_manifest <- data.frame(
  sample_key = rownames(gse_matrix),
  sample_id = gse_samples,
  cohort = "GSE122459",
  clinical_group = ifelse(grepl("^Healthy", gse_samples), "Healthy", "SLE"),
  label = ifelse(grepl("^Healthy", gse_samples), 0, 1),
  subject_id = gse_samples,
  age_group = NA_character_,
  sex = NA_character_,
  stringsAsFactors = FALSE
)

missing_manifest_columns <- setdiff(names(base_manifest), names(gse_manifest))
for (column in missing_manifest_columns) gse_manifest[[column]] <- NA
gse_manifest <- gse_manifest[, names(base_manifest), drop = FALSE]

training_manifest <- rbind(base_manifest, gse_manifest)
training_matrix <- rbind(base_matrix, gse_matrix)
training_matrix <- training_matrix[
  training_manifest$sample_key, gene_pairs, drop = FALSE
]
training_labels <- as.numeric(training_manifest$label)
gene_pair_modes <- as.integer(colMeans(training_matrix == 1L) >= 0.5)
names(gene_pair_modes) <- gene_pairs

cohort_counts <- as.data.frame(table(training_manifest$cohort))
names(cohort_counts) <- c("cohort", "n")
cohort_counts <- cohort_counts[match(development_cohorts, cohort_counts$cohort), ]
rownames(cohort_counts) <- NULL
if (!identical(as.character(cohort_counts$cohort), development_cohorts)) {
  stop("Final training cohorts are not in Table S2 order")
}
if (!identical(as.integer(cohort_counts$n), expected_table_counts$n)) {
  stop("Final training matrix does not match Table S2 cohort counts")
}
stopifnot(
  nrow(training_matrix) == 807L,
  ncol(training_matrix) == 390L,
  all(is.finite(training_matrix)),
  all(training_matrix %in% c(0, 1)),
  identical(
    as.integer(table(factor(training_labels, levels = c(-1, 0, 1)))),
    c(179L, 285L, 343L)
  )
)

set.seed(20261117L)
svr_model <- svm(
  x = training_matrix,
  y = training_labels,
  type = "eps-regression",
  kernel = "radial",
  cost = 10,
  gamma = 0.0001,
  epsilon = 0.1,
  scale = TRUE
)

training_predictions <- as.numeric(predict(svr_model, training_matrix))
if (any(!is.finite(training_predictions))) stop("Non-finite training predictions")

model_metadata <- list(
  model_version = "0.3.0",
  model_status = paste(
    "807-sample refit using the latest 781-sample Top390 feature panel"
  ),
  feature_panel = "latest 781-sample Top390",
  training_samples = 807L,
  training_cohorts = development_cohorts,
  cohort_counts = expected_table_counts,
  label_counts = data.frame(
    immune_category = c("Immunodeficiency", "Healthy", "Autoimmune"),
    label = c(-1L, 0L, 1L),
    n = c(179L, 285L, 343L),
    stringsAsFactors = FALSE
  ),
  feature_count = 390L,
  feature_panel_source = paste(
    "Locked Top390 panel used by Table S6 and the 781-sample model"
  ),
  feature_encoding = "1 when GeneA > GeneB; 0 otherwise, including ties",
  training_unavailable_pair_rule = paste(
    "Pairs unavailable in GSE122459 used the modal binary ordering",
    "estimated from the preceding 781 development samples"
  ),
  unavailable_pair_rule = paste(
    "Pair-specific modal binary ordering estimated from the 807",
    "model-development samples"
  ),
  gse122459_fully_evaluable_pairs = sum(pair_available),
  gse122459_total_pairs = length(pair_available),
  kernel = "radial",
  cost = 10,
  gamma = 0.0001,
  epsilon = 0.1,
  internal_scaling = TRUE,
  gene_count_in_panel = length(unique(unlist(strsplit(gene_pairs, "_vs_", fixed = TRUE)))),
  source = paste(
    "Latest 781-sample Top390 panel refitted on the Table-S2-aligned",
    "807-sample development set"
  ),
  source_table = basename(table_s2_file),
  build_date = as.character(Sys.Date())
)

attr(svr_model, "ImmunoCompass_model_version") <- model_metadata$model_version
attr(svr_model, "ImmunoCompass_training_samples") <- model_metadata$training_samples
attr(svr_model, "ImmunoCompass_feature_encoding") <- model_metadata$feature_encoding

save(svr_model, file = file.path(package_data_dir, "svr_model.rda"), compress = "xz")
save(gene_pairs, file = file.path(package_data_dir, "gene_pairs.rda"), compress = "xz")
save(
  gene_pair_modes,
  file = file.path(package_data_dir, "gene_pair_modes.rda"),
  compress = "xz"
)
save(model_metadata, file = file.path(package_data_dir, "model_metadata.rda"), compress = "xz")

set.seed(20261117L)
example_genes <- unique(unlist(strsplit(gene_pairs, "_vs_", fixed = TRUE)))
example_expression <- matrix(
  round(rnorm(length(example_genes) * 6L, mean = 8, sd = 2), 6),
  nrow = length(example_genes), ncol = 6L,
  dimnames = list(example_genes, paste0("Example", seq_len(6L)))
)
example_metadata <- data.frame(
  sample_id = colnames(example_expression),
  condition = c(0, 0, 0, 1, 1, 1),
  stringsAsFactors = FALSE
)
save(
  example_expression,
  file = file.path(package_data_dir, "example_expression.rda"),
  compress = "xz"
)
save(
  example_metadata,
  file = file.path(package_data_dir, "example_metadata.rda"),
  compress = "xz"
)

if (!is.na(package_extdata_dir)) {
  write.table(
    cbind(gene = rownames(example_expression), example_expression),
    file = file.path(package_extdata_dir, "example_expression.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.csv(
    example_metadata,
    file = file.path(package_extdata_dir, "example_metadata.csv"),
    row.names = FALSE
  )
}

write.csv(training_manifest, file.path(output_dir, "training_manifest_807.csv"), row.names = FALSE)
write.csv(expected_table_counts, file.path(output_dir, "training_cohort_counts.csv"), row.names = FALSE)
write.csv(
  model_metadata$label_counts,
  file.path(output_dir, "training_immune_category_counts.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    gene_pair = gene_pairs,
    GSE122459_evaluable = pair_available,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "GSE122459_pair_coverage.csv"),
  row.names = FALSE
)
write.csv(
  transform(training_manifest, IBI = training_predictions),
  file.path(output_dir, "training_predictions_807.csv"),
  row.names = FALSE
)
save(
  svr_model, gene_pairs, top390, model_metadata, training_manifest,
  file = file.path(output_dir, "final_807_fixed390_svr_model.rda"),
  compress = "xz"
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

cat("Built Table-S2-aligned ImmunoCompass model\n")
cat("  training samples: 807\n")
cat("  cohorts: 12\n")
cat("  labels (-1/0/1): 179/285/343\n")
cat("  features: 390 binary gene pairs\n")
cat("  GSE122459 coverage:", sum(pair_available), "/", length(pair_available), "\n")
cat("  output:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
