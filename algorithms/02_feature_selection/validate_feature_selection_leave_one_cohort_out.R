options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(randomForest)
  library(xgboost)
  library(e1071)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

# Leakage-controlled feature-discovery and leave-one-cohort-out validation for
# the final 711-sample development set. Every outer fold repeats gene coverage
# filtering, all-pairs screening, RF/XGBoost ranking, Top390 selection,
# training-mode estimation, and SVR fitting without the held-out cohort.
root <- Sys.getenv("IBI_TOOL_ROOT", unset = "C:/temp/ibi_tool")
dataset_dir <- file.path(root, "dataset")
input_dir <- file.path(root, "revision_supp_data")
output_dir <- file.path(
  root, "Result", "revision_711_top390_feature_selection_LOCO_20261109"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1L]) else "all"
valid_modes <- c("prepare", "full", "bootstrap", "loco", "summarize", "all")
if (!mode %in% valid_modes) stop("Mode must be one of: ", paste(valid_modes, collapse = ", "))

base_seed <- 20261021L
bootstrap_repeats <- if (length(args) >= 2L) as.integer(args[2L]) else 50L
bootstrap_start <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
bootstrap_end <- if (length(args) >= 4L) as.integer(args[4L]) else bootstrap_repeats
loco_start <- if (mode == "loco" && length(args) >= 2L) as.integer(args[2L]) else 1L
loco_end <- if (mode == "loco" && length(args) >= 3L) as.integer(args[3L]) else Inf
bootstrap_fraction <- 0.80
gene_sample_coverage_threshold <- 0.70
gene_cohort_coverage_fraction <- 0.50
pair_sample_coverage_threshold <- 0.70
minimum_order_fraction <- 0.05
prescreen_count <- 5000L
final_feature_count <- 390L
pair_block_gene_count <- 5L

cache_file <- file.path(output_dir, "development_711_expression_cache.rda")
full_result_file <- file.path(output_dir, "full711_selection_and_model.rda")

load_immune_genes <- function() {
  environment <- new.env(parent = emptyenv())
  load(file.path(root, "all_immune_gene.rda"), envir = environment)
  genes <- unique(toupper(trimws(as.character(environment$all_immune_gene))))
  genes <- genes[nzchar(genes)]
  stopifnot(length(genes) == 3198L)
  genes
}

prepare_expression <- function(expression, genes) {
  expression <- as.matrix(expression)
  storage.mode(expression) <- "double"
  symbols <- toupper(trimws(rownames(expression)))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols %in% genes
  if (!any(keep)) stop("No immune genes found in an expression matrix")
  expression <- expression[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (anyDuplicated(symbols)) {
    expression <- limma::avereps(expression, ID = symbols)
  } else {
    rownames(expression) <- symbols
  }
  expression[intersect(genes, rownames(expression)), , drop = FALSE]
}

read_series <- function(path, characteristic_keys) {
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
    if (is.na(line) || !length(line)) stop("Missing characteristic: ", key)
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    trimws(sub(
      paste0("^", key, ":"), "", sub('^"(.*)"$', "\\1", values),
      ignore.case = TRUE
    ))
  }
  metadata <- data.frame(
    sample_id = extract_field("^!Sample_geo_accession"),
    sample_title = extract_field("^!Sample_title"),
    stringsAsFactors = FALSE
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

map_probes <- function(expression, annotation, genes) {
  mapped_gene <- annotation$gene[match(trimws(rownames(expression)), annotation$probe)]
  keep <- !is.na(mapped_gene) & nzchar(mapped_gene) & mapped_gene %in% genes
  if (!any(keep)) stop("No probes mapped to immune genes")
  mapped <- limma::avereps(expression[keep, , drop = FALSE], ID = mapped_gene[keep])
  prepare_expression(mapped, genes)
}

read_gse59216_batch2 <- function(annotation, genes) {
  path <- file.path(
    dataset_dir, "GSE59216_normalizedData_log2transformed_HT12v3_21samples.csv.gz"
  )
  table <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  probe <- trimws(table[[1L]])
  sample_columns <- !grepl("^detectionPval\\.", names(table))
  sample_columns[1L] <- FALSE
  expression <- as.matrix(table[, sample_columns, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe
  map_probes(expression, annotation, genes)
}

build_development_cache_legacy_unused <- function() {
  immune_genes <- load_immune_genes()
  annotation_v3 <- read_annotation(
    file.path(input_dir, "GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx"),
    "Probe_Id", skip = 8L
  )

  base <- new.env()
  load(file.path(root, "gene_activator_inhibitor.rda"), envir = base)
  raw_original <- list(
    GTEx = base$whole_blood,
    GSE232381 = base$XLA_exp,
    GSE297198 = base$AIDS_exp
  )
  original_objects <- c(
    GSE264125 = "GSE264125_tpm_final",
    GSE222408 = "GSE222408_tpm_final",
    GSE175839 = "GSE175839_tpm_final",
    GSE181228 = "GSE181228_exp",
    GSE199911 = "GSE199911_exp_final",
    GSE140713 = "GSE140713_exp"
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
    environment <- new.env()
    load(file.path(dataset_dir, original_files[[name]]), envir = environment)
    raw_original[[name]] <- environment[[original_objects[[name]]]]
  }
  raw_original$GSE59216_batch2 <- read_gse59216_batch2(annotation_v3, immune_genes)

  label_environment <- new.env()
  load(file.path(root, "y0723.rda"), envir = label_environment)
  original_labels <- as.numeric(label_environment$y)
  original_counts <- vapply(raw_original, ncol, integer(1))
  stopifnot(sum(original_counts) == 599L, length(original_labels) == 599L)

  cohorts <- list()
  start <- 1L
  for (name in names(raw_original)) {
    expression <- prepare_expression(raw_original[[name]], immune_genes)
    indices <- start:(start + ncol(expression) - 1L)
    labels <- original_labels[indices]
    colnames(expression) <- paste(name, colnames(expression), sep = "::")
    cohorts[[name]] <- list(
      expression = expression,
      metadata = data.frame(
        sample_key = colnames(expression),
        sample_id = sub("^[^:]+::", "", colnames(expression)),
        cohort = name,
        clinical_group = c(`-1` = "Immunodeficiency", `0` = "Healthy", `1` = "Autoimmune")[
          as.character(labels)
        ],
        label = labels,
        stringsAsFactors = FALSE
      )
    )
    start <- max(indices) + 1L
  }

  gse51405 <- read_series(
    file.path(input_dir, "GSE51405_series_matrix.txt"), c("illness", "gender", "age")
  )
  expression_51405 <- map_probes(gse51405$expression, annotation_v3, immune_genes)
  group_51405 <- ifelse(
    tolower(gse51405$metadata$illness) == "healthy", "Healthy",
    gse51405$metadata$illness
  )
  labels_51405 <- ifelse(group_51405 == "Healthy", 0, -1)
  colnames(expression_51405) <- paste("GSE51405", colnames(expression_51405), sep = "::")
  cohorts$GSE51405 <- list(
    expression = expression_51405,
    metadata = data.frame(
      sample_key = colnames(expression_51405),
      sample_id = gse51405$metadata$sample_id,
      cohort = "GSE51405",
      clinical_group = group_51405,
      label = labels_51405,
      stringsAsFactors = FALSE
    )
  )

  sample_manifest <- do.call(rbind, lapply(cohorts, `[[`, "metadata"))
  rownames(sample_manifest) <- NULL
  all_expression <- matrix(
    NA_real_, nrow = length(immune_genes), ncol = nrow(sample_manifest),
    dimnames = list(immune_genes, sample_manifest$sample_key)
  )
  for (name in names(cohorts)) {
    expression <- cohorts[[name]]$expression
    all_expression[rownames(expression), colnames(expression)] <- expression
  }
  stopifnot(
    nrow(sample_manifest) == 682L,
    identical(colnames(all_expression), sample_manifest$sample_key),
    all(table(sample_manifest$label) == c(`-1` = 165L, `0` = 194L, `1` = 323L))
  )

  cohort_counts <- as.data.frame(table(
    cohort = sample_manifest$cohort,
    clinical_group = sample_manifest$clinical_group,
    label = sample_manifest$label
  ))
  cohort_counts <- cohort_counts[cohort_counts$Freq > 0L, , drop = FALSE]
  write.csv(sample_manifest, file.path(output_dir, "sample_manifest.csv"), row.names = FALSE)
  write.csv(cohort_counts, file.path(output_dir, "development_cohort_counts.csv"), row.names = FALSE)
  save(
    all_expression, sample_manifest, immune_genes,
    file = cache_file, compress = TRUE
  )
  invisible(cache_file)
}

build_development_cache <- function() {
  source_cache <- file.path(
    root, "Result", "revision_682_feature_stability_LOCO",
    "development_682_expression_cache.rda"
  )
  if (!file.exists(source_cache)) stop("Missing source expression cache: ", source_cache)
  source_environment <- new.env(parent = emptyenv())
  load(source_cache, envir = source_environment)
  required <- c("all_expression", "sample_manifest", "immune_genes")
  if (!all(required %in% ls(source_environment))) {
    stop("Source expression cache is incomplete")
  }

  immune_genes <- source_environment$immune_genes
  source_manifest <- source_environment$sample_manifest
  source_expression <- source_environment$all_expression

  # Replace the legacy 41-sample GTEx block with 70 independent GTEx v10
  # whole-blood donors aged 20-29. All non-GTEx development cohorts,
  # including GSE59216 batch 2, are retained unchanged.
  retain_source <- source_manifest$cohort != "GTEx"
  source_manifest <- source_manifest[retain_source, , drop = FALSE]
  source_expression <- source_expression[, source_manifest$sample_key, drop = FALSE]

  gtex_manifest_file <- file.path(
    input_dir, "GTEx_v10_Whole_Blood_age20_29_TPM_sample_manifest.tsv"
  )
  gtex_gct <- file.path(input_dir, "GTEx_v10_Whole_Blood_age20_29_TPM.gct")
  gtex_manifest <- data.table::fread(gtex_manifest_file, sep = "\t")
  required_manifest_columns <- c(
    "SAMPID", "SUBJID", "SEX", "AGE", "SMTSD", "ANALYTE_TYPE"
  )
  if (!all(required_manifest_columns %in% names(gtex_manifest))) {
    stop("GTEx manifest is missing required columns")
  }
  stopifnot(
    nrow(gtex_manifest) == 70L,
    data.table::uniqueN(gtex_manifest$SAMPID) == 70L,
    data.table::uniqueN(gtex_manifest$SUBJID) == 70L,
    all(gtex_manifest$AGE == "20-29"),
    all(gtex_manifest$SMTSD == "Whole Blood"),
    all(gtex_manifest$ANALYTE_TYPE == "RNA:Total RNA")
  )
  gtex_table <- data.table::fread(
    gtex_gct, sep = "\t", skip = 2L, check.names = FALSE, showProgress = TRUE
  )
  stopifnot(all(c("Name", "Description") %in% names(gtex_table)))
  gtex_samples <- setdiff(names(gtex_table), c("Name", "Description"))
  stopifnot(length(gtex_samples) == 70L, identical(gtex_samples, gtex_manifest$SAMPID))
  gene_symbols <- toupper(trimws(gtex_table$Description))
  keep_rows <- !is.na(gene_symbols) & nzchar(gene_symbols) & gene_symbols %in% immune_genes
  gtex_expression_observed <- as.matrix(gtex_table[keep_rows, ..gtex_samples])
  storage.mode(gtex_expression_observed) <- "double"
  rownames(gtex_expression_observed) <- gene_symbols[keep_rows]
  gtex_expression_observed <- limma::avereps(
    gtex_expression_observed, ID = rownames(gtex_expression_observed)
  )
  gtex_keys <- paste("GTEx", gtex_samples, sep = "::")
  colnames(gtex_expression_observed) <- gtex_keys
  gtex_expression <- matrix(
    NA_real_, nrow = length(immune_genes), ncol = length(gtex_samples),
    dimnames = list(immune_genes, gtex_keys)
  )
  shared_genes <- intersect(immune_genes, rownames(gtex_expression_observed))
  gtex_expression[shared_genes, ] <- gtex_expression_observed[shared_genes, , drop = FALSE]
  gtex_manifest_final <- data.frame(
    sample_key = gtex_keys,
    sample_id = gtex_samples,
    cohort = "GTEx",
    clinical_group = "Healthy",
    label = 0,
    stringsAsFactors = FALSE
  )

  sample_manifest <- rbind(gtex_manifest_final, source_manifest)
  all_expression <- cbind(gtex_expression, source_expression)
  development_cohorts <- c(
    "GTEx", "GSE232381", "GSE297198", "GSE264125", "GSE222408",
    "GSE175839", "GSE181228", "GSE199911", "GSE140713",
    "GSE59216_batch2", "GSE51405"
  )
  cohort_rank <- match(sample_manifest$cohort, development_cohorts)
  sample_manifest <- sample_manifest[order(cohort_rank), , drop = FALSE]
  rownames(sample_manifest) <- NULL
  all_expression <- all_expression[, sample_manifest$sample_key, drop = FALSE]
  stopifnot(
    nrow(sample_manifest) == 711L,
    identical(colnames(all_expression), sample_manifest$sample_key),
    identical(as.character(unique(sample_manifest$cohort)), development_cohorts),
    identical(
      as.integer(table(factor(sample_manifest$label, levels = c(-1, 0, 1)))),
      c(165L, 223L, 323L)
    )
  )

  cohort_counts <- as.data.frame(table(
    cohort = sample_manifest$cohort,
    clinical_group = sample_manifest$clinical_group,
    label = sample_manifest$label
  ))
  cohort_counts <- cohort_counts[cohort_counts$Freq > 0L, , drop = FALSE]
  write.csv(sample_manifest, file.path(output_dir, "sample_manifest.csv"), row.names = FALSE)
  write.csv(cohort_counts, file.path(output_dir, "development_cohort_counts.csv"), row.names = FALSE)
  save(
    all_expression, sample_manifest, immune_genes,
    file = cache_file, compress = TRUE
  )
  invisible(cache_file)
}

ensure_cache <- function() {
  if (!file.exists(cache_file)) build_development_cache()
  load(cache_file, envir = globalenv())
  stopifnot(
    exists("all_expression", envir = globalenv(), inherits = FALSE),
    exists("sample_manifest", envir = globalenv(), inherits = FALSE),
    exists("immune_genes", envir = globalenv(), inherits = FALSE),
    ncol(all_expression) == 711L
  )
}

retain_best <- function(current, incoming, n) {
  combined <- rbind(current, incoming)
  if (!nrow(combined)) return(combined)
  combined <- combined[order(
    -combined$prescreen_score, -combined$sample_coverage, combined$gene_pair
  ), , drop = FALSE]
  combined[seq_len(min(n, nrow(combined))), , drop = FALSE]
}

build_selected_features <- function(expression, candidates, training_modes = NULL) {
  result <- matrix(
    NA_real_, nrow = ncol(expression), ncol = nrow(candidates),
    dimnames = list(colnames(expression), candidates$gene_pair)
  )
  learned_modes <- numeric(nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    a <- expression[candidates$geneA[i], ]
    b <- expression[candidates$geneB[i], ]
    observed <- is.finite(a) & is.finite(b)
    values <- rep(NA_real_, length(a))
    values[observed] <- as.numeric(a[observed] > b[observed])
    if (is.null(training_modes)) {
      mode_value <- if (mean(values[observed] == 1, na.rm = TRUE) >= 0.5) 1 else 0
    } else {
      mode_value <- training_modes[i]
    }
    values[!observed] <- mode_value
    result[, i] <- values
    learned_modes[i] <- mode_value
  }
  list(matrix = result, modes = learned_modes)
}

select_top_pairs <- function(sample_indices, seed, verbose = TRUE) {
  expression <- all_expression[, sample_indices, drop = FALSE]
  metadata <- sample_manifest[sample_indices, , drop = FALSE]
  labels <- metadata$label
  cohort_factor <- factor(metadata$cohort)

  gene_present <- is.finite(expression)
  cohort_gene_presence <- vapply(
    levels(cohort_factor),
    function(name) rowSums(gene_present[, cohort_factor == name, drop = FALSE]) > 0L,
    logical(nrow(expression))
  )
  if (is.null(dim(cohort_gene_presence))) {
    cohort_gene_presence <- matrix(cohort_gene_presence, ncol = 1L)
  }
  minimum_cohorts <- ceiling(nlevels(cohort_factor) * gene_cohort_coverage_fraction)
  gene_coverage <- data.frame(
    gene = immune_genes,
    sample_coverage = rowMeans(gene_present),
    cohorts_present = rowSums(cohort_gene_presence),
    stringsAsFactors = FALSE
  )
  eligible <- gene_coverage$sample_coverage >= gene_sample_coverage_threshold &
    gene_coverage$cohorts_present >= minimum_cohorts
  eligible_genes <- gene_coverage$gene[eligible]
  if (length(eligible_genes) < 30L) stop("Too few genes passed coverage filtering")
  screen_expression <- expression[eligible_genes, , drop = FALSE]

  # Cohort centering makes the univariate prescreen depend on within-cohort
  # immune-state contrasts rather than platform-level mean shifts.
  centered_label <- labels - ave(labels, cohort_factor, FUN = mean)
  centered_label_squared <- centered_label^2
  best_candidates <- data.frame()
  total_pairs_evaluated <- 0L
  n_genes <- length(eligible_genes)
  block_starts <- seq.int(1L, n_genes - 1L, by = pair_block_gene_count)
  for (block_index in seq_along(block_starts)) {
    first <- block_starts[block_index]
    last <- min(first + pair_block_gene_count - 1L, n_genes - 1L)
    gene_a_index <- integer()
    gene_b_index <- integer()
    for (i in first:last) {
      gene_a_index <- c(gene_a_index, rep.int(i, n_genes - i))
      gene_b_index <- c(gene_b_index, seq.int(i + 1L, n_genes))
    }
    a <- screen_expression[gene_a_index, , drop = FALSE]
    b <- screen_expression[gene_b_index, , drop = FALSE]
    observed <- is.finite(a) & is.finite(b)
    z <- matrix(0, nrow = nrow(a), ncol = ncol(a))
    z[a > b & observed] <- 1
    n_observed <- rowSums(observed)
    sum_z <- rowSums(z)
    sum_y <- as.numeric(observed %*% centered_label)
    sum_y2 <- as.numeric(observed %*% centered_label_squared)
    sum_zy <- as.numeric(z %*% centered_label)
    numerator <- n_observed * sum_zy - sum_z * sum_y
    denominator <- sqrt(pmax(
      (n_observed * sum_z - sum_z^2) *
        (n_observed * sum_y2 - sum_y^2),
      0
    ))
    score <- abs(numerator / denominator)
    coverage <- n_observed / ncol(screen_expression)
    positive_fraction <- sum_z / pmax(n_observed, 1L)
    valid <- coverage >= pair_sample_coverage_threshold &
      positive_fraction >= minimum_order_fraction &
      positive_fraction <= (1 - minimum_order_fraction) &
      is.finite(score) & denominator > 0
    if (any(valid)) {
      incoming <- data.frame(
        geneA = eligible_genes[gene_a_index[valid]],
        geneB = eligible_genes[gene_b_index[valid]],
        gene_pair = paste(
          eligible_genes[gene_a_index[valid]], eligible_genes[gene_b_index[valid]],
          sep = "_vs_"
        ),
        prescreen_score = score[valid],
        sample_coverage = coverage[valid],
        positive_order_fraction = positive_fraction[valid],
        stringsAsFactors = FALSE
      )
      best_candidates <- retain_best(best_candidates, incoming, prescreen_count)
    }
    total_pairs_evaluated <- total_pairs_evaluated + length(gene_a_index)
    if (verbose && (block_index %% 100L == 0L || block_index == length(block_starts))) {
      cat(sprintf(
        "  prescreen %d/%d blocks; retained %d candidates\n",
        block_index, length(block_starts), nrow(best_candidates)
      ))
    }
  }
  best_candidates <- best_candidates[order(-best_candidates$prescreen_score), , drop = FALSE]
  best_candidates$prescreen_rank <- seq_len(nrow(best_candidates))

  built <- build_selected_features(expression, best_candidates)
  feature_matrix <- built$matrix
  label_factor <- droplevels(factor(labels, levels = c(-1, 0, 1)))
  class_counts <- table(label_factor)
  if (length(class_counts) < 2L) stop("Feature selection requires at least two labels")

  set.seed(seed)
  rf_model <- randomForest(
    x = feature_matrix, y = label_factor, ntree = 500L,
    mtry = max(1L, floor(sqrt(ncol(feature_matrix)))), importance = TRUE,
    sampsize = rep.int(min(class_counts), length(class_counts)),
    strata = label_factor, replace = TRUE
  )
  rf_importance <- randomForest::importance(rf_model, type = 1, scale = FALSE)
  if (is.matrix(rf_importance)) {
    rf_score <- rf_importance[, "MeanDecreaseAccuracy"]
    names(rf_score) <- rownames(rf_importance)
  } else {
    rf_score <- rf_importance
  }

  label_levels <- sort(unique(labels))
  xgb_label <- match(labels, label_levels) - 1L
  inverse_frequency <- 1 / as.numeric(table(factor(labels, levels = label_levels))[
    as.character(labels)
  ])
  inverse_frequency <- inverse_frequency / mean(inverse_frequency)
  dtrain <- xgboost::xgb.DMatrix(
    data = feature_matrix, label = xgb_label, weight = inverse_frequency
  )
  set.seed(seed)
  xgb_model <- xgboost::xgb.train(
    params = list(
      objective = "multi:softprob", num_class = length(label_levels),
      eval_metric = "mlogloss", eta = 0.05, max_depth = 4L,
      min_child_weight = 2, subsample = 0.8, colsample_bytree = 0.7,
      nthread = 4L
    ),
    data = dtrain, nrounds = 300L, verbose = 0L
  )
  xgb_table <- xgboost::xgb.importance(
    feature_names = colnames(feature_matrix), model = xgb_model
  )
  xgb_score <- setNames(rep(0, ncol(feature_matrix)), colnames(feature_matrix))
  if (nrow(xgb_table)) xgb_score[xgb_table$Feature] <- xgb_table$Gain

  ranking <- best_candidates
  ranking$rf_mean_decrease_accuracy <- rf_score[ranking$gene_pair]
  ranking$xgb_gain <- xgb_score[ranking$gene_pair]
  ranking$rf_rank <- rank(-ranking$rf_mean_decrease_accuracy, ties.method = "average")
  ranking$xgb_rank <- rank(-ranking$xgb_gain, ties.method = "average")
  ranking$integrated_mean_rank <- (ranking$rf_rank + ranking$xgb_rank) / 2
  ranking <- ranking[order(
    ranking$integrated_mean_rank, ranking$rf_rank, ranking$xgb_rank,
    -ranking$prescreen_score, ranking$gene_pair
  ), , drop = FALSE]
  ranking$integrated_rank <- seq_len(nrow(ranking))
  ranking <- ranking[, c(
    "integrated_rank", "gene_pair", "geneA", "geneB",
    "integrated_mean_rank", "rf_rank", "xgb_rank",
    "rf_mean_decrease_accuracy", "xgb_gain", "prescreen_rank",
    "prescreen_score", "sample_coverage", "positive_order_fraction"
  )]
  top390 <- ranking[seq_len(final_feature_count), , drop = FALSE]
  list(
    ranking = ranking,
    top390 = top390,
    training_modes = built$modes[match(top390$gene_pair, best_candidates$gene_pair)],
    eligible_gene_count = length(eligible_genes),
    evaluated_pair_count = total_pairs_evaluated,
    minimum_cohorts = minimum_cohorts,
    rank_correlation = suppressWarnings(cor(
      ranking$rf_rank, ranking$xgb_rank, method = "spearman", use = "complete.obs"
    ))
  )
}

svr_parameters <- function() {
  environment <- new.env()
  load(file.path(root, "best_svr0731.rda"), envir = environment)
  model <- environment$best_svr
  list(
    kernel = c("linear", "polynomial", "radial", "sigmoid")[as.integer(model$kernel) + 1L],
    cost = model$cost, gamma = model$gamma, degree = model$degree,
    coef0 = model$coef0, nu = model$nu, epsilon = model$epsilon
  )
}

train_svr <- function(x, y) {
  parameters <- svr_parameters()
  e1071::svm(
    x = x, y = y, type = "eps-regression", kernel = parameters$kernel,
    cost = parameters$cost, gamma = parameters$gamma,
    degree = parameters$degree, coef0 = parameters$coef0,
    nu = parameters$nu, epsilon = parameters$epsilon, scale = TRUE
  )
}

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

run_full_selection <- function() {
  ensure_cache()
  if (file.exists(full_result_file)) {
    cat("Full selection already exists; keeping checkpoint.\n")
    return(invisible(full_result_file))
  }
  started <- Sys.time()
  cat("Starting full 711-sample feature selection.\n")
  selection <- select_top_pairs(seq_len(nrow(sample_manifest)), base_seed)
  built <- build_selected_features(
    all_expression, selection$top390,
    training_modes = selection$training_modes
  )
  full_feature_matrix <- built$matrix
  full_svr_model <- train_svr(full_feature_matrix, sample_manifest$label)
  elapsed_minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  write.csv(
    selection$top390,
    file.path(output_dir, "full711_top390_gene_pairs.csv"), row.names = FALSE
  )
  write.csv(
    selection$ranking,
    file.path(output_dir, "full711_top5000_RF_XGBoost_ranking.csv"), row.names = FALSE
  )
  save(
    selection, full_feature_matrix, full_svr_model,
    file = full_result_file, compress = TRUE
  )
  metrics <- data.frame(
    samples = nrow(sample_manifest), cohorts = length(unique(sample_manifest$cohort)),
    initial_genes = length(immune_genes), eligible_genes = selection$eligible_gene_count,
    evaluated_pairs = selection$evaluated_pair_count,
    retained_pairs = nrow(selection$top390),
    rf_xgb_spearman = selection$rank_correlation,
    elapsed_minutes = elapsed_minutes
  )
  write.csv(metrics, file.path(output_dir, "full711_selection_metrics.csv"), row.names = FALSE)
  cat(sprintf("Full selection completed in %.2f minutes.\n", elapsed_minutes))
  invisible(full_result_file)
}

stratified_subsample <- function(fraction, seed) {
  set.seed(seed)
  strata <- interaction(sample_manifest$cohort, sample_manifest$label, drop = TRUE)
  unlist(lapply(split(seq_len(nrow(sample_manifest)), strata), function(index) {
    sample(index, max(1L, floor(length(index) * fraction)), replace = FALSE)
  }), use.names = FALSE)
}

run_bootstrap <- function(
    repeats = bootstrap_repeats, start_iteration = 1L, end_iteration = repeats) {
  ensure_cache()
  run_full_selection()
  load(full_result_file)
  bootstrap_dir <- file.path(output_dir, "bootstrap_replicates")
  dir.create(bootstrap_dir, recursive = TRUE, showWarnings = FALSE)
  iterations <- seq.int(max(1L, start_iteration), min(repeats, end_iteration))
  for (iteration in iterations) {
    result_path <- file.path(
      bootstrap_dir, sprintf("bootstrap_%03d_top390.csv", iteration)
    )
    metric_path <- file.path(
      bootstrap_dir, sprintf("bootstrap_%03d_metrics.csv", iteration)
    )
    if (file.exists(result_path) && file.exists(metric_path)) {
      cat(sprintf("Bootstrap %d/%d already exists; skipping.\n", iteration, repeats))
      next
    }
    indices <- stratified_subsample(bootstrap_fraction, base_seed + iteration)
    started <- Sys.time()
    cat(sprintf("Starting bootstrap %d/%d with %d samples.\n", iteration, repeats, length(indices)))
    current <- select_top_pairs(indices, base_seed + 1000L + iteration, verbose = FALSE)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    write.csv(current$top390, result_path, row.names = FALSE)
    write.csv(data.frame(
      iteration = iteration, seed = base_seed + iteration,
      samples = length(indices), cohorts = length(unique(sample_manifest$cohort[indices])),
      eligible_genes = current$eligible_gene_count,
      evaluated_pairs = current$evaluated_pair_count,
      jaccard_with_full = jaccard(current$top390$gene_pair, selection$top390$gene_pair),
      rf_xgb_spearman = current$rank_correlation,
      elapsed_minutes = elapsed
    ), metric_path, row.names = FALSE)
    cat(sprintf("Bootstrap %d completed in %.2f minutes.\n", iteration, elapsed))
  }
  invisible(NULL)
}

summarize_bootstrap <- function(repeats = bootstrap_repeats) {
  run_full_selection()
  load(full_result_file)
  bootstrap_dir <- file.path(output_dir, "bootstrap_replicates")
  result_files <- file.path(
    bootstrap_dir, sprintf("bootstrap_%03d_top390.csv", seq_len(repeats))
  )
  metric_files <- file.path(
    bootstrap_dir, sprintf("bootstrap_%03d_metrics.csv", seq_len(repeats))
  )
  complete <- file.exists(result_files) & file.exists(metric_files)
  if (!any(complete)) return(invisible(NULL))
  result_files <- result_files[complete]
  metric_files <- metric_files[complete]
  sets <- lapply(result_files, function(path) read.csv(path)$gene_pair)
  metrics <- do.call(rbind, lapply(metric_files, read.csv))
  write.csv(metrics, file.path(output_dir, "bootstrap_run_metrics.csv"), row.names = FALSE)

  frequency <- sort(table(unlist(sets)), decreasing = TRUE)
  frequency_table <- data.frame(
    gene_pair = names(frequency), selection_count = as.integer(frequency),
    selection_frequency = as.integer(frequency) / length(sets),
    in_full_top390 = names(frequency) %in% selection$top390$gene_pair,
    stringsAsFactors = FALSE
  )
  full_missing <- setdiff(selection$top390$gene_pair, frequency_table$gene_pair)
  if (length(full_missing)) {
    frequency_table <- rbind(frequency_table, data.frame(
      gene_pair = full_missing, selection_count = 0L, selection_frequency = 0,
      in_full_top390 = TRUE
    ))
  }
  frequency_table <- frequency_table[order(
    -frequency_table$selection_frequency, -frequency_table$in_full_top390,
    frequency_table$gene_pair
  ), , drop = FALSE]
  write.csv(
    frequency_table,
    file.path(output_dir, "bootstrap_feature_selection_frequency.csv"), row.names = FALSE
  )

  pairwise <- list()
  counter <- 1L
  if (length(sets) >= 2L) {
    for (i in seq_len(length(sets) - 1L)) {
      for (j in seq.int(i + 1L, length(sets))) {
        pairwise[[counter]] <- data.frame(
          iteration_1 = metrics$iteration[i], iteration_2 = metrics$iteration[j],
          jaccard = jaccard(sets[[i]], sets[[j]])
        )
        counter <- counter + 1L
      }
    }
  }
  pairwise <- if (length(pairwise)) do.call(rbind, pairwise) else data.frame(
    iteration_1 = integer(), iteration_2 = integer(), jaccard = numeric()
  )
  write.csv(
    pairwise, file.path(output_dir, "bootstrap_pairwise_jaccard.csv"), row.names = FALSE
  )
  summary_table <- data.frame(
    completed_repeats = length(sets),
    median_jaccard_with_full = median(metrics$jaccard_with_full),
    q1_jaccard_with_full = unname(quantile(metrics$jaccard_with_full, 0.25)),
    q3_jaccard_with_full = unname(quantile(metrics$jaccard_with_full, 0.75)),
    median_pairwise_jaccard = if (nrow(pairwise)) median(pairwise$jaccard) else NA_real_,
    q1_pairwise_jaccard = if (nrow(pairwise)) unname(quantile(pairwise$jaccard, 0.25)) else NA_real_,
    q3_pairwise_jaccard = if (nrow(pairwise)) unname(quantile(pairwise$jaccard, 0.75)) else NA_real_,
    features_frequency_ge_0_80 = sum(frequency_table$selection_frequency >= 0.80),
    full_features_frequency_ge_0_80 = sum(
      frequency_table$selection_frequency >= 0.80 & frequency_table$in_full_top390
    ),
    median_rf_xgb_spearman = median(metrics$rf_xgb_spearman, na.rm = TRUE)
  )
  write.csv(
    summary_table, file.path(output_dir, "bootstrap_jaccard_summary.csv"), row.names = FALSE
  )
  invisible(summary_table)
}

calculate_loco_metrics <- function(cohort, observed, predicted) {
  output <- data.frame(
    cohort = cohort, n = length(observed),
    n_immunodeficiency = sum(observed == -1),
    n_healthy = sum(observed == 0),
    n_autoimmune = sum(observed == 1),
    rmse = sqrt(mean((predicted - observed)^2)),
    mae = mean(abs(predicted - observed)),
    spearman_rho = if (length(unique(observed)) > 1L) {
      suppressWarnings(cor(observed, predicted, method = "spearman"))
    } else NA_real_,
    auc = NA_real_, auc_ci_low = NA_real_, auc_ci_high = NA_real_,
    auc_contrast = NA_character_, stringsAsFactors = FALSE
  )
  nonzero <- setdiff(unique(observed), 0)
  if (0 %in% observed && length(nonzero) == 1L) {
    disease <- observed != 0
    directional_score <- if (nonzero == 1) predicted else -predicted
    roc_object <- pROC::roc(
      response = disease, predictor = directional_score,
      levels = c(FALSE, TRUE), direction = "<", quiet = TRUE
    )
    interval <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
    output$auc <- as.numeric(pROC::auc(roc_object))
    output$auc_ci_low <- interval[1L]
    output$auc_ci_high <- interval[3L]
    output$auc_contrast <- if (nonzero == 1) {
      "Autoimmune vs healthy (higher IBI indicates disease)"
    } else {
      "Immunodeficiency vs healthy (lower IBI indicates disease)"
    }
  }
  output
}

run_loco <- function(start_fold = 1L, end_fold = Inf) {
  ensure_cache()
  run_full_selection()
  load(full_result_file)
  loco_dir <- file.path(output_dir, "loco_folds")
  dir.create(loco_dir, recursive = TRUE, showWarnings = FALSE)
  cohorts <- unique(sample_manifest$cohort)
  fold_indices <- seq.int(max(1L, start_fold), min(length(cohorts), end_fold))
  for (fold in fold_indices) {
    held_out <- cohorts[fold]
    metric_path <- file.path(loco_dir, paste0(held_out, "_metrics.csv"))
    prediction_path <- file.path(loco_dir, paste0(held_out, "_predictions.csv"))
    feature_path <- file.path(loco_dir, paste0(held_out, "_top390.csv"))
    if (all(file.exists(c(metric_path, prediction_path, feature_path)))) {
      cat("LOCO fold already exists; skipping ", held_out, ".\n", sep = "")
      next
    }
    test_indices <- which(sample_manifest$cohort == held_out)
    train_indices <- setdiff(seq_len(nrow(sample_manifest)), test_indices)
    started <- Sys.time()
    cat(sprintf(
      "Starting LOCO %d/%d: %s (%d training, %d held out).\n",
      fold, length(cohorts), held_out, length(train_indices), length(test_indices)
    ))
    current <- select_top_pairs(train_indices, base_seed + 2000L + fold, verbose = FALSE)
    training_built <- build_selected_features(
      all_expression[, train_indices, drop = FALSE], current$top390
    )
    test_built <- build_selected_features(
      all_expression[, test_indices, drop = FALSE], current$top390,
      training_modes = training_built$modes
    )
    model <- train_svr(training_built$matrix, sample_manifest$label[train_indices])
    predicted <- as.numeric(predict(model, test_built$matrix))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    metrics <- calculate_loco_metrics(
      held_out, sample_manifest$label[test_indices], predicted
    )
    metrics$training_n <- length(train_indices)
    metrics$eligible_genes <- current$eligible_gene_count
    metrics$jaccard_with_full <- jaccard(
      current$top390$gene_pair, selection$top390$gene_pair
    )
    metrics$rf_xgb_spearman <- current$rank_correlation
    metrics$elapsed_minutes <- elapsed
    predictions <- sample_manifest[test_indices, , drop = FALSE]
    predictions$IBI <- predicted
    predictions$residual <- predicted - predictions$label
    write.csv(metrics, metric_path, row.names = FALSE)
    write.csv(predictions, prediction_path, row.names = FALSE)
    write.csv(current$top390, feature_path, row.names = FALSE)
    cat(sprintf("LOCO %s completed in %.2f minutes.\n", held_out, elapsed))
  }
  invisible(summarize_loco())
}

summarize_loco <- function() {
  if (!file.exists(full_result_file)) return(invisible(NULL))
  load(full_result_file)
  loco_dir <- file.path(output_dir, "loco_folds")
  metric_files <- list.files(loco_dir, pattern = "_metrics\\.csv$", full.names = TRUE)
  prediction_files <- list.files(loco_dir, pattern = "_predictions\\.csv$", full.names = TRUE)
  feature_files <- list.files(loco_dir, pattern = "_top390\\.csv$", full.names = TRUE)
  if (length(metric_files)) {
    metrics <- do.call(rbind, lapply(metric_files, read.csv))
    write.csv(metrics, file.path(output_dir, "loco_performance_metrics.csv"), row.names = FALSE)
  }
  if (length(prediction_files)) {
    predictions <- do.call(rbind, lapply(prediction_files, read.csv))
    write.csv(predictions, file.path(output_dir, "loco_sample_predictions.csv"), row.names = FALSE)
  }
  if (length(feature_files)) {
    overlap <- do.call(rbind, lapply(feature_files, function(path) {
      features <- read.csv(path)$gene_pair
      data.frame(
        held_out_cohort = sub("_top390\\.csv$", "", basename(path)),
        shared_with_full = length(intersect(features, selection$top390$gene_pair)),
        jaccard_with_full = jaccard(features, selection$top390$gene_pair)
      )
    }))
    write.csv(overlap, file.path(output_dir, "loco_feature_overlap.csv"), row.names = FALSE)
  }
  invisible(NULL)
}

make_plots <- function() {
  bootstrap_metrics_path <- file.path(output_dir, "bootstrap_run_metrics.csv")
  frequency_path <- file.path(output_dir, "bootstrap_feature_selection_frequency.csv")
  pairwise_path <- file.path(output_dir, "bootstrap_pairwise_jaccard.csv")
  loco_metrics_path <- file.path(output_dir, "loco_performance_metrics.csv")
  predictions_path <- file.path(output_dir, "loco_sample_predictions.csv")
  qa_png_dir <- Sys.getenv("IBI_QA_PNG_DIR", unset = "")
  save_qa_png <- function(plot, filename, width_mm, height_mm) {
    if (!nzchar(qa_png_dir)) return(invisible(NULL))
    dir.create(qa_png_dir, recursive = TRUE, showWarnings = FALSE)
    grDevices::png(
      file.path(qa_png_dir, filename), width = width_mm / 25.4,
      height = height_mm / 25.4, units = "in", res = 150,
      type = "cairo", family = "Arial"
    )
    print(plot)
    dev.off()
    invisible(NULL)
  }

  theme_set(theme_classic(base_size = 8, base_family = "Arial") + theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    strip.text = element_text(size = 7, face = "bold"),
    plot.title = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 7), legend.text = element_text(size = 6.5)
  ))

  if (all(file.exists(c(bootstrap_metrics_path, frequency_path, pairwise_path)))) {
    metrics <- read.csv(bootstrap_metrics_path)
    frequency <- read.csv(frequency_path)
    pairwise <- read.csv(pairwise_path)
    frequency$rank <- seq_len(nrow(frequency))
    plot_a <- ggplot(frequency, aes(rank, selection_frequency, colour = in_full_top390)) +
      geom_point(size = 0.8, alpha = 0.75) +
      geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.35) +
      scale_colour_manual(values = c(`FALSE` = "#9A9A9A", `TRUE` = "#B24A4A")) +
      labs(
        title = "a  Feature-selection frequency", x = "Gene-pair rank by frequency",
        y = "Selection frequency", colour = "Full-data top390"
      )
    plot_b <- ggplot(metrics, aes(jaccard_with_full)) +
      geom_histogram(binwidth = 0.025, boundary = 0, fill = "#4C78A8", colour = "white") +
      labs(
        title = "b  Agreement with full-data top390", x = "Jaccard index", y = "Resamples"
      )
    plot_c <- ggplot(pairwise, aes(jaccard)) +
      geom_histogram(binwidth = 0.025, boundary = 0, fill = "#72A0A8", colour = "white") +
      labs(
        title = "c  Pairwise resample agreement", x = "Jaccard index", y = "Resample pairs"
      )
    combined <- plot_a / (plot_b | plot_c) + plot_layout(heights = c(1.15, 1))
    cairo_pdf(
      file.path(output_dir, "feature_stability_summary.pdf"),
      width = 180 / 25.4, height = 150 / 25.4, family = "Arial"
    )
    print(combined)
    dev.off()
    save_qa_png(combined, "feature_stability_summary.png", 180, 150)
  }

  if (file.exists(loco_metrics_path)) {
    metrics <- read.csv(loco_metrics_path)
    auc_data <- metrics[is.finite(metrics$auc), , drop = FALSE]
    if (nrow(auc_data)) {
      auc_data$cohort <- factor(
        auc_data$cohort, levels = auc_data$cohort[order(auc_data$auc)]
      )
      plot_auc <- ggplot(auc_data, aes(auc, cohort)) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#999999", linewidth = 0.35) +
        geom_errorbar(
          aes(xmin = auc_ci_low, xmax = auc_ci_high), width = 0.18,
          linewidth = 0.45, colour = "#4C78A8", orientation = "y"
        ) +
        geom_point(size = 2.1, colour = "#B24A4A") +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        labs(
          title = "LOCO discrimination by completely held-out cohort",
          x = "AUC (95% CI)", y = NULL
        )
      cairo_pdf(
        file.path(output_dir, "LOCO_AUC_summary.pdf"),
        width = 125 / 25.4, height = max(75, 15 * nrow(auc_data)) / 25.4,
        family = "Arial"
      )
      print(plot_auc)
      dev.off()
      save_qa_png(
        plot_auc, "LOCO_AUC_summary.png", 125, max(75, 15 * nrow(auc_data))
      )
    }
    metric_long <- rbind(
      data.frame(cohort = metrics$cohort, metric = "RMSE", value = metrics$rmse),
      data.frame(cohort = metrics$cohort, metric = "MAE", value = metrics$mae)
    )
    plot_error <- ggplot(metric_long, aes(reorder(cohort, value), value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.72), width = 0.64) +
      coord_flip() +
      scale_fill_manual(values = c(RMSE = "#4C78A8", MAE = "#72A0A8")) +
      labs(
        title = "LOCO prediction error by held-out cohort",
        x = NULL, y = "Error", fill = NULL
      )
    cairo_pdf(
      file.path(output_dir, "LOCO_error_summary.pdf"),
      width = 125 / 25.4, height = max(85, 10 * nrow(metrics)) / 25.4,
      family = "Arial"
    )
    print(plot_error)
    dev.off()
    save_qa_png(
      plot_error, "LOCO_error_summary.png", 125, max(85, 10 * nrow(metrics))
    )
  }

  if (file.exists(predictions_path)) {
    predictions <- read.csv(predictions_path)
    predictions$label_group <- factor(
      predictions$label, levels = c(-1, 0, 1),
      labels = c("Immunodeficiency", "Healthy", "Autoimmune")
    )
    plot_predictions <- ggplot(predictions, aes(label_group, IBI, colour = label_group)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "#999999", linewidth = 0.3) +
      geom_boxplot(
        width = 0.58, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.35
      ) +
      geom_point(
        position = position_jitter(width = 0.12, height = 0, seed = base_seed),
        size = 0.75, alpha = 0.65
      ) +
      facet_wrap(~ cohort, scales = "free_x", ncol = 4) +
      scale_colour_manual(
        values = c(Immunodeficiency = "#315A7D", Healthy = "#777777", Autoimmune = "#B24A4A"),
        drop = FALSE, guide = "none"
      ) +
      labs(
        title = "IBI predictions from leave-one-cohort-out validation",
        x = NULL, y = "Predicted IBI"
      ) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
    cairo_pdf(
      file.path(output_dir, "LOCO_prediction_distributions.pdf"),
      width = 183 / 25.4, height = 145 / 25.4, family = "Arial"
    )
    print(plot_predictions)
    dev.off()
    save_qa_png(plot_predictions, "LOCO_prediction_distributions.png", 183, 145)
  }
}

write_provenance <- function() {
  parameters <- svr_parameters()
  provenance <- data.frame(
    item = c(
      "Development samples", "Development cohorts", "Labels", "Initial immune genes",
      "Gene-pair encoding", "Gene sample coverage threshold",
      "Gene cohort coverage threshold", "Pair sample coverage threshold",
      "Minimum order fraction", "Univariate prescreen", "Feature ranking",
      "Final feature count", "Stability resampling", "LOCO definition",
      "Missing-pair handling", "SVR kernel", "SVR hyperparameters", "Random seed"
    ),
    value = c(
      "Final 711-sample development set: 165 immunodeficiency, 223 healthy, and 323 autoimmune samples",
      "11 cohorts; GSE59216 91-sample batch in development and 21-sample batch excluded from development",
      "-1 immunodeficiency; 0 healthy; 1 autoimmune",
      "3198 immune-related genes",
      "1 when gene A expression > gene B expression; otherwise 0",
      gene_sample_coverage_threshold,
      paste0("Present in at least ceiling(0.50 x number of training cohorts)"),
      pair_sample_coverage_threshold, minimum_order_fraction,
      "Top 5000 by absolute Pearson association with cohort-centered immune label",
      "Balanced Random Forest (500 trees) and weighted multiclass XGBoost (300 rounds); ascending mean rank",
      final_feature_count,
      paste0("", round(100 * bootstrap_fraction), "% within each cohort-by-label stratum; complete feature selection repeated"),
      "One entire study cohort excluded before coverage filtering, feature selection, imputation, and SVR fitting",
      "Training-set modal ordering; held-out data never determine imputation",
      parameters$kernel,
      paste0("cost=", parameters$cost, "; gamma=", parameters$gamma,
             "; degree=", parameters$degree, "; epsilon=", parameters$epsilon),
      base_seed
    ),
    stringsAsFactors = FALSE
  )
  write.csv(provenance, file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
}

if (mode == "prepare") {
  build_development_cache()
} else if (mode == "full") {
  ensure_cache()
  run_full_selection()
} else if (mode == "bootstrap") {
  ensure_cache()
  run_bootstrap(bootstrap_repeats, bootstrap_start, bootstrap_end)
} else if (mode == "loco") {
  ensure_cache()
  run_loco(loco_start, loco_end)
} else if (mode == "summarize") {
  ensure_cache()
  run_full_selection()
  summarize_bootstrap(bootstrap_repeats)
  summarize_loco()
  make_plots()
} else if (mode == "all") {
  ensure_cache()
  run_full_selection()
  run_bootstrap(bootstrap_repeats)
  summarize_bootstrap(bootstrap_repeats)
  run_loco()
  make_plots()
}
write_provenance()
cat("Completed mode:", mode, "\nOutput:", output_dir, "\n")
