options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(randomForest)
  library(xgboost)
})

root <- "C:/temp/ibi_tool"
dataset_dir <- file.path(root, "dataset")
input_dir <- file.path(root, "revision_supp_data")
output_dir <- file.path(root, "Result", "revision_all_samples_de_novo_top390")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261018L
gene_sample_coverage_threshold <- 0.70
gene_cohort_coverage_fraction <- 0.50
pair_sample_coverage_threshold <- 0.70
minimum_order_fraction <- 0.05
prescreen_count <- 5000L
final_feature_count <- 390L
pair_block_gene_count <- 5L

load(file.path(root, "all_immune_gene.rda"))
immune_genes <- unique(toupper(trimws(as.character(all_immune_gene))))
immune_genes <- immune_genes[nzchar(immune_genes)]
stopifnot(length(immune_genes) == 3198L)

prepare_expression <- function(expression, genes = immune_genes) {
  expression <- as.matrix(expression)
  storage.mode(expression) <- "double"
  symbols <- toupper(trimws(rownames(expression)))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols %in% genes
  expression <- expression[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (anyDuplicated(symbols)) expression <- avereps(expression, ID = symbols)
  else rownames(expression) <- symbols
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
    )]
    if (!length(line)) stop("Missing characteristic: ", key, " in ", basename(path))
    values <- strsplit(line[1L], "\t", fixed = TRUE)[[1L]][-1L]
    trimws(sub(
      paste0("^", key, ":"), "", sub('^"(.*)"$', "\\1", values),
      ignore.case = TRUE
    ))
  }
  metadata <- data.frame(
    sample_id = extract_field("^!Sample_geo_accession"),
    sample_title = extract_field("^!Sample_title")
  )
  for (key in characteristic_keys) metadata[[gsub(" ", "_", key, fixed = TRUE)]] <-
    extract_characteristic(key)
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
    gene = toupper(trimws(annotation[["Symbol"]]))
  )
}

map_probes <- function(expression, annotation) {
  mapped_gene <- annotation$gene[match(trimws(rownames(expression)), annotation$probe)]
  keep <- !is.na(mapped_gene) & nzchar(mapped_gene) & mapped_gene %in% immune_genes
  mapped <- avereps(expression[keep, , drop = FALSE], ID = mapped_gene[keep])
  prepare_expression(mapped)
}

read_gse59216_batch2 <- function(annotation) {
  path <- file.path(
    dataset_dir,
    "GSE59216_normalizedData_log2transformed_HT12v3_21samples.csv.gz"
  )
  table <- fread(path, data.table = FALSE, showProgress = FALSE)
  probe <- trimws(table[[1L]])
  sample_columns <- !grepl("^detectionPval\\.", names(table))
  sample_columns[1L] <- FALSE
  expression <- as.matrix(table[, sample_columns, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe
  map_probes(expression, annotation)
}

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
      paste0('^!Sample_characteristics_ch1\\t"', key, ':'), lines,
      ignore.case = TRUE
    )][1L]
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    trimws(sub(
      paste0("^", key, ":"), "", sub('^"(.*)"$', "\\1", values),
      ignore.case = TRUE
    ))
  }
  metadata <- data.frame(
    sample_id = extract_field("^!Sample_geo_accession"),
    subject = extract_characteristic("subject"),
    visit = suppressWarnings(as.numeric(extract_characteristic("visit"))),
    set = extract_characteristic("set"),
    disease_state = extract_characteristic("disease state")
  )
  metadata <- metadata[metadata$set != "Technical_Replicate", , drop = FALSE]
  metadata <- metadata[order(
    metadata$subject, ifelse(is.na(metadata$visit), Inf, metadata$visit),
    metadata$sample_id
  ), , drop = FALSE]
  metadata <- metadata[!duplicated(metadata$subject), , drop = FALSE]

  target_annotation <- annotation[
    annotation$gene %in% immune_genes & nzchar(annotation$probe), , drop = FALSE
  ]
  raw <- fread(
    path, skip = length(header) - 1L, sep = "\t", header = TRUE,
    quote = "\"", data.table = FALSE, showProgress = TRUE
  )
  names(raw)[1L] <- "probe_id"
  keep <- trimws(raw$probe_id) %in% target_annotation$probe
  probe <- trimws(raw$probe_id[keep])
  expression <- as.matrix(raw[keep, metadata$sample_id, drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- probe
  rm(raw)
  gc()
  list(expression = map_probes(expression, target_annotation), metadata = metadata)
}

cohorts <- list()
add_cohort <- function(name, expression, label, clinical_group = NULL) {
  expression <- prepare_expression(expression)
  stopifnot(ncol(expression) == length(label), all(label %in% c(-1, 0, 1)))
  if (is.null(clinical_group)) clinical_group <- as.character(label)
  stopifnot(length(clinical_group) == length(label))
  colnames(expression) <- paste(name, colnames(expression), sep = "::")
  cohorts[[name]] <<- list(
    expression = expression,
    metadata = data.frame(
      sample_key = colnames(expression), sample_id = sub("^[^:]+::", "", colnames(expression)),
      cohort = name, clinical_group = as.character(clinical_group),
      label = as.numeric(label)
    )
  )
}

annotation_v3 <- read_annotation(
  file.path(input_dir, "GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx"),
  "Probe_Id", skip = 8L
)
annotation_v4 <- read_annotation(
  file.path(input_dir, "GPL10558-50081.txt"), "ID"
)

# Reconstruct the original 599-sample development compendium from raw expression.
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
  environment <- new.env()
  load(file.path(dataset_dir, original_files[[name]]), envir = environment)
  raw_original[[name]] <- environment[[ls(environment)[1L]]]
}
raw_original$GSE59216_batch2 <- read_gse59216_batch2(annotation_v3)
load(file.path(root, "y0723.rda"))
original_counts <- vapply(raw_original, ncol, integer(1))
stopifnot(sum(original_counts) == length(y), length(y) == 599L)
label_start <- 1L
for (name in names(raw_original)) {
  indices <- label_start:(label_start + original_counts[[name]] - 1L)
  add_cohort(name, raw_original[[name]], as.numeric(y[indices]))
  label_start <- max(indices) + 1L
}
rm(raw_original, base)
gc()

# Add all currently available cohorts with unambiguous -1/0/1 labels.
gse51405 <- read_series(
  file.path(input_dir, "GSE51405_series_matrix.txt"), c("illness", "gender", "age")
)
gse51404 <- read_series(
  file.path(input_dir, "GSE51404_series_matrix.txt"), c("illness", "gender", "age")
)
gse44969 <- read_series(
  file.path(input_dir, "GSE44969_series_matrix.txt"), c("group", "age", "sex")
)
expression_51405 <- map_probes(gse51405$expression, annotation_v3)
expression_51404 <- map_probes(gse51404$expression, annotation_v4)
expression_44969 <- map_probes(gse44969$expression, annotation_v3)
group_51405 <- ifelse(tolower(gse51405$metadata$illness) == "healthy", "Healthy", "Immunodeficiency")
label_51405 <- ifelse(group_51405 == "Healthy", 0, -1)
group_51404 <- ifelse(
  grepl("_XLA$", gse51404$metadata$sample_title), "XLA",
  ifelse(tolower(gse51404$metadata$illness) == "healthy", "Healthy", gse51404$metadata$illness)
)
label_51404 <- ifelse(group_51404 == "Healthy", 0, -1)
group_44969 <- ifelse(
  gse44969$metadata$group == "Healthy control", "Healthy",
  ifelse(gse44969$metadata$group == "Chronic granulomatous patient", "CGD", "Carrier")
)
keep_44969 <- group_44969 != "Carrier"
add_cohort("GSE51405", expression_51405, label_51405, group_51405)
add_cohort("GSE51404", expression_51404, label_51404, group_51404)
add_cohort(
  "GSE44969", expression_44969[, keep_44969, drop = FALSE],
  ifelse(group_44969[keep_44969] == "Healthy", 0, -1), group_44969[keep_44969]
)
rm(gse51405, gse51404, gse44969, expression_51405, expression_51404, expression_44969)
gc()

gse65391 <- read_gse65391_targeted(annotation_v4)
label_65391 <- ifelse(gse65391$metadata$disease_state == "SLE", 1, 0)
add_cohort(
  "GSE65391", gse65391$expression, label_65391,
  gse65391$metadata$disease_state
)
rm(gse65391)
gc()

environment <- new.env()
load(file.path(dataset_dir, "GSE53419_tpm_final.rda"), envir = environment)
expression_53419 <- environment$GSE53419_tpm_final
group_53419 <- ifelse(grepl("^Control_", colnames(expression_53419)), "Healthy", "Autoimmune")
add_cohort("GSE53419", expression_53419, ifelse(group_53419 == "Healthy", 0, 1), group_53419)

environment <- new.env()
load(file.path(dataset_dir, "GSE59216_91_exp_final.rda"), envir = environment)
expression_59216 <- environment$GSE59216_91_exp_final
load(file.path(dataset_dir, "GSE59216_91_group_final.rda"))
group_59216 <- ifelse(as.character(group) == "Control", "Healthy", "22q11DS")
add_cohort("GSE59216_batch1", expression_59216, ifelse(group_59216 == "Healthy", 0, -1), group_59216)

environment <- new.env()
load(file.path(dataset_dir, "GSE221786_exp_tpm.txt"), envir = environment)
expression_221786 <- environment$GSE221786_exp_tpm
group_221786 <- ifelse(grepl("^Ctrl_", colnames(expression_221786)), "Healthy", "Autoimmune")
add_cohort("GSE221786", expression_221786, ifelse(group_221786 == "Healthy", 0, 1), group_221786)
rm(environment, expression_53419, expression_59216, expression_221786)
gc()

sample_manifest <- do.call(rbind, lapply(cohorts, `[[`, "metadata"))
rownames(sample_manifest) <- NULL
write.csv(sample_manifest, file.path(output_dir, "all_samples_manifest.csv"), row.names = FALSE)

all_expression <- matrix(
  NA_real_, nrow = length(immune_genes), ncol = nrow(sample_manifest),
  dimnames = list(immune_genes, sample_manifest$sample_key)
)
for (name in names(cohorts)) {
  expression <- cohorts[[name]]$expression
  all_expression[rownames(expression), colnames(expression)] <- expression
}
labels <- sample_manifest$label
cohort_factor <- factor(sample_manifest$cohort, levels = names(cohorts))
stopifnot(identical(colnames(all_expression), sample_manifest$sample_key))

gene_present <- is.finite(all_expression)
cohort_gene_presence <- vapply(
  levels(cohort_factor),
  function(name) rowSums(gene_present[, cohort_factor == name, drop = FALSE]) > 0L,
  logical(nrow(all_expression))
)
gene_coverage <- data.frame(
  gene = immune_genes,
  sample_coverage = rowMeans(gene_present),
  cohorts_present = rowSums(cohort_gene_presence),
  total_cohorts = nlevels(cohort_factor)
)
minimum_cohorts <- ceiling(nlevels(cohort_factor) * gene_cohort_coverage_fraction)
gene_coverage$eligible <- gene_coverage$sample_coverage >= gene_sample_coverage_threshold &
  gene_coverage$cohorts_present >= minimum_cohorts
write.csv(gene_coverage, file.path(output_dir, "immune_gene_coverage_filter.csv"), row.names = FALSE)
eligible_genes <- gene_coverage$gene[gene_coverage$eligible]
if (length(eligible_genes) < 30L) stop("Too few genes passed the coverage filter")
screen_expression <- all_expression[eligible_genes, , drop = FALSE]

# Center labels within cohorts so the prescreen emphasizes within-cohort biology.
within_cohort_label <- labels - ave(labels, cohort_factor, FUN = mean)
within_cohort_label_squared <- within_cohort_label^2

retain_best <- function(current, incoming, n) {
  combined <- rbind(current, incoming)
  if (!nrow(combined)) return(combined)
  combined <- combined[order(-combined$prescreen_score, -combined$sample_coverage), , drop = FALSE]
  combined[seq_len(min(n, nrow(combined))), , drop = FALSE]
}

best_candidates <- data.frame()
total_pairs_evaluated <- 0
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
  z <- matrix(-1, nrow = nrow(a), ncol = ncol(a))
  z[a > b & observed] <- 1
  z[!observed] <- 0
  n_observed <- rowSums(observed)
  positive_fraction <- rowSums(z == 1 & observed) / pmax(n_observed, 1L)
  sxy <- as.numeric(z %*% within_cohort_label)
  syy <- as.numeric(observed %*% within_cohort_label_squared)
  score <- abs(sxy) / sqrt(pmax(n_observed * syy, .Machine$double.eps))
  coverage <- n_observed / ncol(screen_expression)
  valid <- coverage >= pair_sample_coverage_threshold &
    positive_fraction >= minimum_order_fraction &
    positive_fraction <= (1 - minimum_order_fraction) &
    is.finite(score) & syy > 0
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
      positive_order_fraction = positive_fraction[valid]
    )
    best_candidates <- retain_best(best_candidates, incoming, prescreen_count)
  }
  total_pairs_evaluated <- total_pairs_evaluated + length(gene_a_index)
  if (block_index %% 25L == 0L || block_index == length(block_starts)) {
    cat(sprintf(
      "Prescreen block %d/%d; %.1f%%; retained %d candidates\n",
      block_index, length(block_starts), 100 * block_index / length(block_starts),
      nrow(best_candidates)
    ))
  }
  rm(a, b, observed, z)
}
best_candidates <- best_candidates[!duplicated(best_candidates$gene_pair), , drop = FALSE]
best_candidates <- best_candidates[order(-best_candidates$prescreen_score), , drop = FALSE]
best_candidates$prescreen_rank <- seq_len(nrow(best_candidates))
write.csv(
  best_candidates, file.path(output_dir, "top5000_univariate_prescreen.csv"),
  row.names = FALSE
)

build_selected_features <- function(expression, candidates) {
  result <- matrix(
    NA_integer_, nrow = ncol(expression), ncol = nrow(candidates),
    dimnames = list(colnames(expression), candidates$gene_pair)
  )
  for (i in seq_len(nrow(candidates))) {
    a <- expression[candidates$geneA[i], ]
    b <- expression[candidates$geneB[i], ]
    observed <- is.finite(a) & is.finite(b)
    values <- ifelse(a > b, 1L, -1L)
    modal <- ifelse(mean(values[observed] == 1L) >= 0.5, 1L, -1L)
    values[!observed] <- modal
    result[, i] <- values
  }
  result
}

feature_matrix <- build_selected_features(all_expression, best_candidates)
label_factor <- factor(labels, levels = c(-1, 0, 1))
class_counts <- table(label_factor)

set.seed(seed)
rf_model <- randomForest(
  x = feature_matrix, y = label_factor, ntree = 500L,
  mtry = max(1L, floor(sqrt(ncol(feature_matrix)))), importance = TRUE,
  sampsize = rep.int(min(class_counts), length(class_counts)),
  strata = label_factor, replace = TRUE
)
rf_importance <- importance(rf_model, type = 1, scale = FALSE)
rf_score <- rf_importance[, "MeanDecreaseAccuracy"]
names(rf_score) <- rownames(rf_importance)

xgb_label <- match(labels, c(-1, 0, 1)) - 1L
class_weight <- 1 / as.numeric(class_counts[as.character(label_factor)])
class_weight <- class_weight / mean(class_weight)
dtrain <- xgb.DMatrix(data = feature_matrix, label = xgb_label, weight = class_weight)
set.seed(seed)
xgb_model <- xgb.train(
  params = list(
    objective = "multi:softprob", num_class = 3L, eval_metric = "mlogloss",
    eta = 0.05, max_depth = 4L, min_child_weight = 2,
    subsample = 0.8, colsample_bytree = 0.7, nthread = 4L,
    seed = seed
  ),
  data = dtrain, nrounds = 300L, verbose = 1L, print_every_n = 25L
)
xgb_table <- xgb.importance(feature_names = colnames(feature_matrix), model = xgb_model)
xgb_score <- setNames(rep(0, ncol(feature_matrix)), colnames(feature_matrix))
xgb_score[xgb_table$Feature] <- xgb_table$Gain

ranking <- best_candidates
ranking$rf_mean_decrease_accuracy <- rf_score[ranking$gene_pair]
ranking$xgb_gain <- xgb_score[ranking$gene_pair]
ranking$rf_rank <- rank(-ranking$rf_mean_decrease_accuracy, ties.method = "average")
ranking$xgb_rank <- rank(-ranking$xgb_gain, ties.method = "average")
ranking$integrated_mean_rank <- (ranking$rf_rank + ranking$xgb_rank) / 2
ranking <- ranking[order(
  ranking$integrated_mean_rank, ranking$rf_rank, ranking$xgb_rank,
  -ranking$prescreen_score
), , drop = FALSE]
ranking$integrated_rank <- seq_len(nrow(ranking))
ranking <- ranking[, c(
  "integrated_rank", "gene_pair", "geneA", "geneB",
  "integrated_mean_rank", "rf_rank", "xgb_rank",
  "rf_mean_decrease_accuracy", "xgb_gain", "prescreen_rank",
  "prescreen_score", "sample_coverage", "positive_order_fraction"
)]
write.csv(ranking, file.path(output_dir, "all_5000_RF_XGBoost_importance.csv"), row.names = FALSE)

top390 <- ranking[seq_len(final_feature_count), , drop = FALSE]
write.csv(top390, file.path(output_dir, "de_novo_top390_gene_pairs.csv"), row.names = FALSE)
features_top390 <- top390$gene_pair
top390_feature_matrix <- feature_matrix[, features_top390, drop = FALSE]
save(
  features_top390, top390, top390_feature_matrix, sample_manifest,
  file = file.path(output_dir, "de_novo_top390_features_and_matrix.rda")
)
save(
  rf_model, xgb_model,
  file = file.path(output_dir, "feature_ranking_RF_XGBoost_models.rda")
)

provenance <- data.frame(
  item = c(
    "Analysis status", "Initial immune genes", "Theoretical all gene pairs",
    "Eligible genes after coverage filtering", "Pairs evaluated after gene coverage filtering",
    "Samples included", "Cohorts included", "Labels", "Gene-pair encoding",
    "Gene sample coverage threshold", "Gene cohort coverage threshold",
    "Pair sample coverage threshold", "Minimum order fraction",
    "Univariate prescreen", "Machine-learning ranking", "Integrated ranking",
    "Final feature count", "Data leakage caveat", "Random seed"
  ),
  value = c(
    "Exploratory de novo feature reselection using all currently available labelled samples",
    length(immune_genes), choose(length(immune_genes), 2), length(eligible_genes),
    total_pairs_evaluated, nrow(sample_manifest), nlevels(cohort_factor),
    "-1 immunodeficiency; 0 healthy; 1 autoimmune",
    "1 when gene A expression > gene B expression; otherwise -1",
    gene_sample_coverage_threshold,
    paste0("Present in at least ", minimum_cohorts, " of ", nlevels(cohort_factor), " cohorts"),
    pair_sample_coverage_threshold, minimum_order_fraction,
    "Top 5000 by absolute association with cohort-centered ordered immune label",
    "Balanced random forest (500 trees) and weighted multiclass XGBoost (300 rounds)",
    "Ascending mean of RF and XGBoost ranks", final_feature_count,
    "All included cohorts participated in feature selection and cannot be treated as independent validation",
    seed
  )
)
write.csv(provenance, file.path(output_dir, "analysis_provenance.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
cat("Completed. Output:", output_dir, "\n")
