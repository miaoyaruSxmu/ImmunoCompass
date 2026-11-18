options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(e1071)
  library(pROC)
  library(ggplot2)
})

# This analysis compares three representations while keeping the samples,
# endpoint, SVR kernel and hyperparameters identical. The 390 gene pairs are
# loaded from the previously generated feature set; no feature reselection is
# performed here.
root <- "C:/temp/ibi_tool"
use_final_711 <- identical(Sys.getenv("ENCODING_COMPARISON_711", unset = "0"), "1")
manifest_file <- file.path(
  root, "Result", "revision_all_samples_de_novo_top390",
  "de_novo_top390_features_and_matrix.rda"
)
dataset_dir <- file.path(root, "dataset")
input_dir <- file.path(root, "revision_supp_data")
output_dir <- file.path(
  root, "Result",
  if (use_final_711) "revision_711_encoding_comparison_fixed390_20261110" else
    "revision_encoding_comparison_original390"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261018L
set.seed(seed)

load(manifest_file)
load(file.path(root, "X0723.rda"))
load(file.path(root, "y0723.rda"))
load(file.path(root, "best_svr0731.rda"))
pair_names <- colnames(best_svr[["SV"]])
stopifnot(length(pair_names) == 390L, all(pair_names %in% colnames(X)),
          nrow(X) == 599L, length(y) == 599L)
original_binary_input <- as.matrix(X[, pair_names, drop = FALSE])
pair_genes <- unique(unlist(strsplit(pair_names, "_vs_", fixed = TRUE)))

old_model <- best_svr
kernel_name <- function(x) c("linear", "polynomial", "radial", "sigmoid")[as.integer(x) + 1L]
svr_kernel <- kernel_name(old_model$kernel)
svr_params <- list(
  kernel = svr_kernel, degree = old_model$degree, gamma = old_model$gamma,
  coef0 = old_model$coef0, cost = old_model$cost, nu = old_model$nu,
  epsilon = old_model$epsilon
)

if (use_final_711) {
  final_711_dir <- file.path(
    root, "Result", "revision_711_top390_feature_selection_LOCO_20261109"
  )
  cache_env <- new.env(parent = emptyenv())
  load(file.path(final_711_dir, "development_711_expression_cache.rda"), envir = cache_env)
  model_env <- new.env(parent = emptyenv())
  load(file.path(final_711_dir, "full711_selection_and_model.rda"), envir = model_env)
  development_711_expression <- cache_env$all_expression
  development_711_manifest <- cache_env$sample_manifest
  original_binary_input <- as.matrix(model_env$full_feature_matrix)
  old_model <- model_env$full_svr_model
  pair_names <- colnames(original_binary_input)
  pair_genes <- unique(unlist(strsplit(pair_names, "_vs_", fixed = TRUE)))
  svr_kernel <- kernel_name(old_model$kernel)
  svr_params <- list(
    kernel = svr_kernel, degree = old_model$degree, gamma = old_model$gamma,
    coef0 = old_model$coef0, cost = old_model$cost, nu = old_model$nu,
    epsilon = old_model$epsilon
  )
  stopifnot(
    nrow(original_binary_input) == 711L,
    nrow(development_711_manifest) == 711L,
    length(pair_names) == 390L,
    identical(rownames(original_binary_input), development_711_manifest$sample_key)
  )
}

prepare_expression <- function(expression, genes = pair_genes) {
  expression <- as.matrix(expression)
  storage.mode(expression) <- "double"
  symbols <- toupper(trimws(rownames(expression)))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols %in% genes
  if (!any(keep)) {
    return(matrix(NA_real_, nrow = length(genes), ncol = ncol(expression),
                  dimnames = list(genes, colnames(expression))))
  }
  expression <- expression[keep, , drop = FALSE]
  symbols <- symbols[keep]
  mapped <- if (anyDuplicated(symbols)) avereps(expression, ID = symbols) else {
    rownames(expression) <- symbols
    expression
  }
  out <- matrix(NA_real_, nrow = length(genes), ncol = ncol(expression),
                dimnames = list(genes, colnames(expression)))
  present <- intersect(genes, rownames(mapped))
  out[present, ] <- mapped[present, , drop = FALSE]
  out
}

read_series <- function(path, characteristic_keys = character()) {
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
  header <- lines[seq_len(begin)]
  extract_field <- function(prefix) {
    line <- header[grep(prefix, header)[1L]]
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    sub('^"(.*)"$', "\\1", values)
  }
  characteristic_lines <- header[grep("^!Sample_characteristics_ch1", header)]
  extract_characteristic <- function(key) {
    line <- characteristic_lines[grep(
      paste0('^!Sample_characteristics_ch1\\t"', key, ':'),
      characteristic_lines, ignore.case = TRUE
    )][1L]
    if (is.na(line) || !length(line)) return(rep(NA_character_, ncol(expression)))
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    trimws(sub(
      paste0("^", key, ":"), "", sub('^"(.*)"$', "\\1", values),
      ignore.case = TRUE
    ))
  }
  metadata <- data.frame(
    sample_id = extract_field("^!Sample_geo_accession"),
    sample_title = extract_field("^!Sample_title"), stringsAsFactors = FALSE
  )
  for (key in characteristic_keys) {
    metadata[[gsub(" ", "_", key, fixed = TRUE)]] <- extract_characteristic(key)
  }
  stopifnot(identical(colnames(expression), metadata$sample_id))
  list(expression = expression, metadata = metadata)
}

read_annotation <- function(path, probe_column, skip = 0L) {
  annotation <- read.delim(
    path, skip = skip, comment.char = if (skip == 0L) "#" else "",
    quote = "", check.names = FALSE, fill = TRUE
  )
  data.frame(
    probe = trimws(annotation[[probe_column]]),
    gene = toupper(trimws(annotation[["Symbol"]])),
    stringsAsFactors = FALSE
  )
}

map_probes <- function(expression, annotation) {
  mapped_gene <- annotation$gene[match(trimws(rownames(expression)), annotation$probe)]
  keep <- !is.na(mapped_gene) & nzchar(mapped_gene) & mapped_gene %in% pair_genes
  if (!any(keep)) {
    return(matrix(NA_real_, nrow = length(pair_genes), ncol = ncol(expression),
                  dimnames = list(pair_genes, colnames(expression))))
  }
  mapped <- avereps(expression[keep, , drop = FALSE], ID = mapped_gene[keep])
  prepare_expression(mapped)
}

read_gse59216_batch2 <- function(annotation) {
  path <- file.path(dataset_dir,
                    "GSE59216_normalizedData_log2transformed_HT12v3_21samples.csv.gz")
  tab <- fread(path, data.table = FALSE, showProgress = FALSE)
  probe <- trimws(tab[[1L]])
  sample_columns <- !grepl("^detectionPval\\.", names(tab))
  sample_columns[1L] <- FALSE
  expression <- as.matrix(tab[, sample_columns, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe
  map_probes(expression, annotation)
}

# GSE65391 is large; read only probes mapping to genes used in the 390 pairs.
read_gse65391_targeted <- function(annotation) {
  path <- file.path(input_dir, "GSE65391_series_matrix.txt")
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  header <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("GSE65391 table marker not found")
    header <- c(header, line)
    if (grepl("^!series_matrix_table_begin", line)) break
  }
  close(con)
  on.exit(NULL, add = FALSE)
  extract_field <- function(prefix) {
    line <- header[grep(prefix, header)[1L]]
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    sub('^"(.*)"$', "\\1", values)
  }
  extract_characteristic <- function(key) {
    lines <- header[grep("^!Sample_characteristics_ch1", header)]
    line <- lines[grep(
      paste0('^!Sample_characteristics_ch1\\t"', key, ':'),
      lines, ignore.case = TRUE
    )][1L]
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    trimws(sub(paste0("^", key, ":"), "",
                sub('^"(.*)"$', "\\1", values), ignore.case = TRUE))
  }
  metadata <- data.frame(
    sample_id = extract_field("^!Sample_geo_accession"),
    sample_title = extract_field("^!Sample_title"),
    subject = extract_characteristic("subject"),
    visit = suppressWarnings(as.numeric(extract_characteristic("visit"))),
    set = extract_characteristic("set"),
    disease_state = extract_characteristic("disease state"),
    stringsAsFactors = FALSE
  )
  technical <- metadata$set == "Technical_Replicate"
  keep_meta <- !technical & metadata$disease_state %in% c("Healthy", "SLE") &
    !is.na(metadata$subject)
  eligible <- which(keep_meta)
  eligible <- eligible[order(
    metadata$subject[eligible],
    ifelse(is.na(metadata$visit[eligible]), Inf, metadata$visit[eligible]),
    metadata$sample_id[eligible]
  )]
  primary <- eligible[!duplicated(metadata$subject[eligible])]
  sample_ids <- metadata$sample_id[primary]
  # The saved de-novo feature manifest is authoritative for the sample keys.
  # The original GSE65391 import used a column-selection path that can shift
  # the retained GEO accessions by one position; using the saved keys keeps
  # the continuous-input comparison aligned to the exact binary-feature set.
  saved_index <- which(sample_manifest$cohort == "GSE65391")
  if (length(saved_index)) {
    sample_ids <- sub("^GSE65391::", "", sample_manifest$sample_key[saved_index])
    primary <- seq_along(sample_ids)
    metadata <- data.frame(
      sample_id = sample_ids,
      sample_title = NA_character_,
      subject = NA_character_, visit = NA_real_, set = NA_character_,
      disease_state = sample_manifest$clinical_group[saved_index],
      stringsAsFactors = FALSE
    )
  }
  target <- annotation$probe[annotation$gene %in% pair_genes & nzchar(annotation$probe)]
  table_begin <- grep("^!series_matrix_table_begin", header)
  raw <- fread(
    path, skip = table_begin, sep = "\t", header = TRUE, quote = "\"",
    data.table = FALSE, select = c(1L, match(sample_ids, extract_field("^!Sample_geo_accession"))),
    showProgress = TRUE
  )
  names(raw)[1L] <- "probe_id"
  raw <- raw[trimws(raw$probe_id) %in% target, , drop = FALSE]
  probe <- trimws(raw$probe_id)
  expression <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe
  stopifnot(ncol(expression) == length(sample_ids))
  colnames(expression) <- sample_ids
  mapped <- map_probes(expression, annotation)
  metadata <- metadata[primary, , drop = FALSE]
  rownames(metadata) <- NULL
  list(expression = mapped, metadata = metadata)
}

cohort_expr <- list()
cohort_meta <- list()
add_cohort <- function(name, expression, label, clinical_group = NULL) {
  expression <- prepare_expression(expression)
  if (is.null(clinical_group)) clinical_group <- as.character(label)
  stopifnot(ncol(expression) == length(label), length(clinical_group) == length(label))
  sample_id <- colnames(expression)
  sample_key <- paste(name, sample_id, sep = "::")
  colnames(expression) <- sample_key
  cohort_expr[[name]] <<- expression
  cohort_meta[[name]] <<- data.frame(
    sample_key = sample_key, sample_id = sample_id, cohort = name,
    clinical_group = as.character(clinical_group), label = as.numeric(label),
    stringsAsFactors = FALSE
  )
}

# Reconstruct the same 599-sample development compendium used by the saved
# feature matrix, preserving the original sample order and labels.
base <- new.env()
load(file.path(root, "gene_activator_inhibitor.rda"), envir = base)
raw_original <- list(
  GTEx = base$whole_blood,
  GSE232381 = base$XLA_exp,
  GSE297198 = base$AIDS_exp
)
original_files <- c(
  GSE264125 = "GSE264125_tpm_final.rda",
  GSE222408 = "GSE222408_tpm_final.rda",
  GSE175839 = "GSE175839_tpm_final.rda",
  GSE181228 = "GSE181228_exp.rda",
  GSE199911 = "GSE199911_exp_final.rda",
  GSE140713 = "GSE140713_exp.rda"
)
for (name in names(original_files)) {
  env <- new.env()
  load(file.path(dataset_dir, original_files[[name]]), envir = env)
  raw_original[[name]] <- env[[ls(env)[1L]]]
}
annotation_v3 <- read_annotation(
  file.path(input_dir, "GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx"),
  "Probe_Id", skip = 8L
)
annotation_v4 <- read_annotation(file.path(input_dir, "GPL10558-50081.txt"), "ID")
raw_original$GSE59216_batch2 <- read_gse59216_batch2(annotation_v3)
load(file.path(root, "y0723.rda"))
original_counts <- vapply(raw_original, ncol, integer(1))
stopifnot(sum(original_counts) == length(y), length(y) == 599L)
label_start <- 1L
for (name in names(raw_original)) {
  count <- original_counts[[name]]
  indices <- label_start:(label_start + count - 1L)
  # The historical data objects already have gene symbols for the original
  # development matrix; process them exactly as in feature screening.
  clinical <- as.character(y[indices])
  add_cohort(name, raw_original[[name]], as.numeric(y[indices]), clinical)
  label_start <- max(indices) + 1L
}
rm(raw_original, base)
gc()

# Cohorts added during revision and external cohorts.
gse51405 <- read_series(file.path(input_dir, "GSE51405_series_matrix.txt"),
                        c("illness", "gender", "age"))
gse51404 <- read_series(file.path(input_dir, "GSE51404_series_matrix.txt"),
                        c("illness", "gender", "age"))
gse44969 <- read_series(file.path(input_dir, "GSE44969_series_matrix.txt"),
                        c("group", "age", "sex"))
expr_51405 <- map_probes(gse51405$expression, annotation_v3)
expr_51404 <- map_probes(gse51404$expression, annotation_v4)
expr_44969 <- map_probes(gse44969$expression, annotation_v3)
group_51405 <- ifelse(tolower(gse51405$metadata$illness) == "healthy",
                      "Healthy", "Immunodeficiency")
group_51404 <- ifelse(
  grepl("_XLA$", gse51404$metadata$sample_title), "XLA",
  ifelse(tolower(gse51404$metadata$illness) == "healthy", "Healthy",
         gse51404$metadata$illness)
)
group_44969 <- ifelse(gse44969$metadata$group == "Healthy control", "Healthy",
                      ifelse(gse44969$metadata$group == "Chronic granulomatous patient",
                             "CGD", "Carrier"))
add_cohort("GSE51405", expr_51405,
           ifelse(group_51405 == "Healthy", 0, -1), group_51405)
add_cohort("GSE51404", expr_51404,
           ifelse(group_51404 == "Healthy", 0, -1), group_51404)
keep_44969 <- group_44969 != "Carrier"
add_cohort("GSE44969", expr_44969[, keep_44969, drop = FALSE],
           ifelse(group_44969[keep_44969] == "Healthy", 0, -1),
           group_44969[keep_44969])

gse65391 <- read_gse65391_targeted(annotation_v4)
add_cohort("GSE65391", gse65391$expression,
           ifelse(gse65391$metadata$disease_state == "SLE", 1, 0),
           gse65391$metadata$disease_state)

env <- new.env()
load(file.path(dataset_dir, "GSE53419_tpm_final.rda"), envir = env)
expr_53419 <- env$GSE53419_tpm_final
group_53419 <- ifelse(grepl("^Control_", colnames(expr_53419)), "Healthy", "Autoimmune")
add_cohort("GSE53419", expr_53419, ifelse(group_53419 == "Healthy", 0, 1), group_53419)

env <- new.env()
load(file.path(dataset_dir, "GSE59216_91_exp_final.rda"), envir = env)
expr_59216 <- env$GSE59216_91_exp_final
load(file.path(dataset_dir, "GSE59216_91_group_final.rda"))
group_59216 <- ifelse(as.character(group) == "Control", "Healthy", "22q11DS")
add_cohort("GSE59216_batch1", expr_59216,
           ifelse(group_59216 == "Healthy", 0, -1), group_59216)

env <- new.env()
load(file.path(dataset_dir, "GSE221786_exp_tpm.txt"), envir = env)
expr_221786 <- env$GSE221786_exp_tpm
group_221786 <- ifelse(grepl("^Ctrl_", colnames(expr_221786)), "Healthy", "Autoimmune")
add_cohort("GSE221786", expr_221786,
           ifelse(group_221786 == "Healthy", 0, 1), group_221786)

sample_manifest_rebuilt <- do.call(rbind, cohort_meta)
rownames(sample_manifest_rebuilt) <- NULL
stopifnot(nrow(sample_manifest_rebuilt) == nrow(sample_manifest))
manifest_order <- match(sample_manifest$sample_key, sample_manifest_rebuilt$sample_key)
if (anyNA(manifest_order)) {
  writeLines(sample_manifest$sample_key[is.na(manifest_order)],
             file.path(output_dir, "manifest_missing_keys.txt"))
  writeLines(setdiff(sample_manifest_rebuilt$sample_key, sample_manifest$sample_key),
             file.path(output_dir, "manifest_extra_keys.txt"))
  stop("Some saved sample keys were not reconstructed")
}
sample_manifest_rebuilt <- sample_manifest_rebuilt[manifest_order, , drop = FALSE]
rownames(sample_manifest_rebuilt) <- NULL

all_expression <- matrix(
  NA_real_, nrow = length(pair_genes), ncol = nrow(sample_manifest_rebuilt),
  dimnames = list(pair_genes, sample_manifest_rebuilt$sample_key)
)
for (name in names(cohort_expr)) {
  e <- cohort_expr[[name]]
  all_expression[rownames(e), colnames(e)] <- e
}
stopifnot(identical(colnames(all_expression), sample_manifest_rebuilt$sample_key))

if (use_final_711) {
  external_cohorts_711 <- c(
    "GSE53419", "GSE221786", "GSE59216_batch1",
    "GSE51404", "GSE44969", "GSE65391"
  )
  external_rows <- sample_manifest_rebuilt$cohort %in% external_cohorts_711
  external_manifest <- sample_manifest_rebuilt[external_rows, , drop = FALSE]
  external_expression <- all_expression[, external_rows, drop = FALSE]
  development_expression <- development_711_expression[pair_genes, , drop = FALSE]
  all_expression <- cbind(development_expression, external_expression)
  sample_manifest_rebuilt <- rbind(
    development_711_manifest[, names(external_manifest), drop = FALSE],
    external_manifest
  )
  rownames(sample_manifest_rebuilt) <- NULL
  colnames(all_expression) <- sample_manifest_rebuilt$sample_key
  stopifnot(
    ncol(development_expression) == 711L,
    identical(colnames(all_expression), sample_manifest_rebuilt$sample_key)
  )
}

original_cohorts <- c(
  "GTEx", "GSE232381", "GSE297198", "GSE264125", "GSE222408",
  "GSE175839", "GSE181228", "GSE199911", "GSE140713", "GSE59216_batch2"
)
external_cohorts <- c(
  "GSE53419", "GSE221786", "GSE59216_batch1", "GSE51404", "GSE44969", "GSE65391"
)
training <- if (use_final_711) seq_len(nrow(sample_manifest_rebuilt)) <= 711L else
  sample_manifest_rebuilt$cohort %in% c(original_cohorts, "GSE51405")
validation <- sample_manifest_rebuilt$cohort %in% external_cohorts
stopifnot(
  sum(training) == if (use_final_711) 711L else 682L,
  sum(validation) == 413L, !any(training & validation)
)
stopifnot(all(sample_manifest_rebuilt$label[training] %in% c(-1, 0, 1)))
training_y <- sample_manifest_rebuilt$label[training]

# Prepare three representations from the same gene-level matrix.
# 1) Absolute expression: all unique genes involved in the 390 pairs.
#    Expression is log2-transformed after a single non-negative shift estimated
#    from the training samples only.
# 2) Continuous log ratio: log2-transformed expression of gene A minus gene B.
# 3) Binary pair feature: 1 if A > B, otherwise 0 (the requested 0/1 coding).
training_values <- all_expression[, training, drop = FALSE]
finite_training_values <- training_values[is.finite(training_values)]
absolute_shift <- if (length(finite_training_values)) {
  max(0, -min(finite_training_values))
} else 0
absolute_input <- t(log2(pmax(all_expression + absolute_shift, 0) + 1))
dimnames(absolute_input) <- list(colnames(all_expression), pair_genes)
absolute_input <- as.matrix(absolute_input)

logratio_input <- matrix(NA_real_, nrow = nrow(absolute_input), ncol = length(pair_names),
                         dimnames = list(rownames(absolute_input), pair_names))
for (i in seq_along(pair_names)) {
  genes <- strsplit(pair_names[i], "_vs_", fixed = TRUE)[[1L]]
  a <- absolute_input[, genes[1L]]
  b <- absolute_input[, genes[2L]]
  logratio_input[, i] <- a - b
}
# Use the original 390-pair matrix for the 599 development samples, changing
# only its coding from -1/1 to 0/1. This preserves the original feature
# orientation and screening-time missing-pair handling. For all other samples,
# the pair order is computed from the reconstructed gene-level matrix below.
binary_input <- matrix(NA_real_, nrow = nrow(absolute_input), ncol = length(pair_names),
                       dimnames = list(rownames(absolute_input), pair_names))
if (all(original_binary_input %in% c(0, 1, NA))) {
  binary_input[seq_len(nrow(original_binary_input)), ] <- original_binary_input
} else {
  binary_input[seq_len(nrow(original_binary_input)), ] <- (original_binary_input + 1) / 2
}
for (i in seq_along(pair_names)) {
  genes <- strsplit(pair_names[i], "_vs_", fixed = TRUE)[[1L]]
  idx <- (nrow(original_binary_input) + 1L):nrow(binary_input)
  a <- absolute_input[idx, genes[1L]]
  b <- absolute_input[idx, genes[2L]]
  binary_input[idx, i] <- ifelse(is.finite(a) & is.finite(b) & a > b, 1, 0)
}
rownames(binary_input) <- sample_manifest_rebuilt$sample_key
colnames(binary_input) <- pair_names

# Impute missing continuous gene features using training medians only. Binary
# features use the training-set modal ordering for unavailable gene pairs.
fit_impute_continuous <- function(x_train, x_all) {
  med <- apply(x_train, 2L, function(v) {
    v <- v[is.finite(v)]
    if (length(v)) median(v) else 0
  })
  out <- x_all
  for (j in seq_len(ncol(out))) out[!is.finite(out[, j]), j] <- med[j]
  list(matrix = out, medians = med)
}
fit_impute_binary <- function(x_train, x_all) {
  modes <- apply(x_train, 2L, function(v) if (mean(v == 1, na.rm = TRUE) >= 0.5) 1 else 0)
  out <- x_all
  for (j in seq_len(ncol(out))) out[!is.finite(out[, j]), j] <- modes[j]
  list(matrix = out, modes = modes)
}

absolute_fit <- fit_impute_continuous(absolute_input[training, , drop = FALSE], absolute_input)
logratio_fit <- fit_impute_continuous(logratio_input[training, , drop = FALSE], logratio_input)
binary_fit <- fit_impute_binary(binary_input[training, , drop = FALSE], binary_input)

fit_svr <- function(x, y) {
  e1071::svm(
    x = x, y = y, type = "eps-regression", kernel = svr_params$kernel,
    degree = svr_params$degree, gamma = svr_params$gamma,
    coef0 = svr_params$coef0, cost = svr_params$cost, nu = svr_params$nu,
    epsilon = svr_params$epsilon, scale = TRUE
  )
}
models <- list(
  absolute_expression = fit_svr(absolute_fit$matrix[training, , drop = FALSE], training_y),
  continuous_log_ratio = fit_svr(logratio_fit$matrix[training, , drop = FALSE], training_y),
  binary_gene_pair_01 = fit_svr(binary_fit$matrix[training, , drop = FALSE], training_y)
)

predictions <- data.frame(
  sample_key = sample_manifest_rebuilt$sample_key,
  sample_id = sample_manifest_rebuilt$sample_id,
  cohort = sample_manifest_rebuilt$cohort,
  clinical_group = sample_manifest_rebuilt$clinical_group,
  true_label = sample_manifest_rebuilt$label,
  set = ifelse(training, "Training", ifelse(validation, "External validation", "Not evaluated")),
  stringsAsFactors = FALSE
)
input_matrices <- list(
  absolute_expression = absolute_fit$matrix,
  continuous_log_ratio = logratio_fit$matrix,
  binary_gene_pair_01 = binary_fit$matrix
)
for (method in names(models)) {
  predictions[[paste0("pred_", method)]] <- as.numeric(
    predict(models[[method]], input_matrices[[method]])
  )
}
write.csv(predictions, file.path(output_dir, "sample_level_predictions.csv"), row.names = FALSE)

safe_spearman <- function(x, y) {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}
safe_auc <- function(score, y, positive_label) {
  keep <- is.finite(score) & !is.na(y) & y %in% c(0, positive_label)
  if (sum(y[keep] == positive_label) < 2L || sum(y[keep] == 0) < 2L) return(c(auc = NA, ci_low = NA, ci_high = NA))
  roc_obj <- pROC::roc(as.integer(y[keep] == positive_label), score[keep],
                       levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  c(auc = as.numeric(pROC::auc(roc_obj)), ci_low = ci[1L], ci_high = ci[3L])
}

metric_rows <- list()
for (cohort in external_cohorts) {
  idx <- which(validation & sample_manifest_rebuilt$cohort == cohort)
  y_obs <- predictions$true_label[idx]
  groups <- predictions$clinical_group[idx]
  for (method in names(models)) {
    score <- predictions[idx, paste0("pred_", method)]
    # Overall ordinal performance against the predefined -1/0/1 endpoint.
    rmse <- sqrt(mean((score - y_obs)^2, na.rm = TRUE))
    rho <- safe_spearman(score, y_obs)
    # Disease-vs-healthy AUC: autoimmune is positive for autoimmune cohorts;
    # immunodeficiency is positive for immunodeficiency cohorts.
    disease_direction <- ifelse(any(y_obs == 1), 1, -1)
    auc_score <- if (disease_direction == -1) -score else score
    auc <- safe_auc(auc_score, y_obs, disease_direction)
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      cohort = cohort, comparison = ifelse(disease_direction == 1,
        "Autoimmune vs Healthy", "Immunodeficiency vs Healthy"),
      method = method, n = length(idx), n_case = sum(y_obs == disease_direction),
      n_healthy = sum(y_obs == 0), RMSE = rmse, Spearman_rho = rho,
      score_direction = ifelse(disease_direction == 1,
                               "Higher score indicates case",
                               "Lower score indicates case"),
      AUC = auc[["auc"]], AUC_CI_low = auc[["ci_low"]],
      AUC_CI_high = auc[["ci_high"]], stringsAsFactors = FALSE
    )
    # Subgroup AUCs for the multi-group CVID/XLA cohort.
    if (cohort == "GSE51404") {
      for (subgroup in setdiff(unique(groups), "Healthy")) {
        keep <- groups %in% c("Healthy", subgroup)
        sub_auc <- safe_auc(-score[keep], ifelse(groups[keep] == subgroup, -1, 0), -1)
        metric_rows[[length(metric_rows) + 1L]] <- data.frame(
          cohort = cohort, comparison = paste(subgroup, "vs Healthy"),
          method = method, n = sum(keep), n_case = sum(groups == subgroup),
          n_healthy = sum(groups == "Healthy"), RMSE = NA_real_,
          Spearman_rho = NA_real_,
          score_direction = "Lower score indicates case",
          AUC = sub_auc[["auc"]],
          AUC_CI_low = sub_auc[["ci_low"]], AUC_CI_high = sub_auc[["ci_high"]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
metrics <- do.call(rbind, metric_rows)
write.csv(metrics, file.path(output_dir, "external_validation_metrics.csv"), row.names = FALSE)
primary_metrics <- metrics[!is.na(metrics$RMSE), , drop = FALSE]
write.csv(primary_metrics, file.path(output_dir, "external_validation_primary_metrics.csv"), row.names = FALSE)
method_summary <- do.call(rbind, lapply(split(primary_metrics, primary_metrics$method), function(d) {
  data.frame(
    method = as.character(d$method[1L]),
    external_cohorts = nrow(d),
    mean_AUC = mean(d$AUC, na.rm = TRUE),
    median_AUC = median(d$AUC, na.rm = TRUE),
    mean_RMSE = mean(d$RMSE, na.rm = TRUE),
    mean_Spearman_rho = mean(d$Spearman_rho, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(method_summary, file.path(output_dir, "external_validation_method_summary.csv"), row.names = FALSE)

training_metrics <- do.call(rbind, lapply(names(models), function(method) {
  score <- predictions[training, paste0("pred_", method)]
  data.frame(
    cohort = "Training", comparison = "Ordinal endpoint", method = method,
    n = sum(training), n_case = NA_integer_, n_healthy = NA_integer_,
    RMSE = sqrt(mean((score - training_y)^2)),
    Spearman_rho = safe_spearman(score, training_y),
    AUC = NA_real_, AUC_CI_low = NA_real_, AUC_CI_high = NA_real_
  )
}))
write.csv(training_metrics, file.path(output_dir, "training_metrics.csv"), row.names = FALSE)

# One compact PDF per external cohort, showing the three prediction
# distributions on the same samples. The primary numerical results remain in
# external_validation_metrics.csv.
plot_data <- do.call(rbind, lapply(names(models), function(method) {
  data.frame(
    cohort = predictions$cohort[validation],
    clinical_group = predictions$clinical_group[validation],
    method = method,
    prediction = predictions[validation, paste0("pred_", method)]
  )
}))
for (cohort in external_cohorts) {
  d <- plot_data[plot_data$cohort == cohort, , drop = FALSE]
  d$method <- factor(d$method, levels = names(models), labels = c(
    "Absolute expression", "Continuous log ratio", "Binary gene pair (0/1)"
  ))
  d$clinical_group <- factor(d$clinical_group, levels = unique(d$clinical_group))
  p <- ggplot(d, aes(clinical_group, prediction, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.78), width = 0.65,
                 outlier.shape = NA, colour = "black", linewidth = 0.3) +
    geom_point(position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.78),
               size = 0.8, alpha = 0.55) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888", linewidth = 0.3) +
    scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B")) +
    labs(title = paste0("Three input representations: ", cohort),
         x = NULL, y = "Predicted immune balance score", fill = "Input") +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 9),
          axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "bottom")
  cairo_pdf(file.path(output_dir, paste0(cohort, "_three_input_comparison.pdf")),
            width = 150 / 25.4, height = 105 / 25.4, family = "Arial")
  print(p)
  dev.off()
}

write.csv(
  data.frame(
    item = c(
      "Feature set", "Feature count", "Absolute-expression dimension",
      "Absolute-expression definition",
      "Continuous-ratio definition", "Binary-pair definition", "Training cohorts",
      "Training samples", "Training label counts", "External validation cohorts",
      "SVR kernel", "SVR hyperparameters", "Continuous-feature imputation",
      "Binary-feature imputation", "Feature reselection", "Output"
    ),
    value = c(
      if (use_final_711)
        "Final locked 390-pair set selected from the 711-sample development set; no reselection"
      else "Original 390-pair set from best_svr0731/X0723; no reselection",
      "390 gene pairs",
      paste(length(pair_genes), "unique genes represented in the pair list"),
      "log2(expression after a training-derived non-negative shift + 1)",
      "log2-transformed gene A expression minus gene B expression for each pair",
      "1 if gene A > gene B, otherwise 0",
      if (use_final_711) "Final 711-sample development set across 11 cohorts" else
        "Original development cohorts (599) + complete GSE51405",
      as.character(sum(training)),
      paste(names(table(training_y)), as.integer(table(training_y)), collapse = "; "),
      paste(external_cohorts, collapse = "; "), svr_kernel,
      paste(names(svr_params)[-1], unlist(svr_params)[-1], sep = "=", collapse = "; "),
      "Training-set median per gene feature", "Training-set modal value per pair",
      "No; the previously generated top390 set was used unchanged", "PDF figures and CSV result tables"
    ), stringsAsFactors = FALSE
  ), file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE
)
write.csv(
  data.frame(cohort = names(cohort_expr), n = vapply(cohort_expr, ncol, integer(1)),
             clinical_groups = vapply(cohort_meta, function(x) paste(unique(x$clinical_group), collapse = "; "), character(1))),
  file.path(output_dir, "cohort_manifest.csv"), row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("Completed three-input SVR comparison:", output_dir, "\n")
