suppressPackageStartupMessages({
  library(edgeR)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

# Independent validation of the manuscript-reported CAS biomarkers.
# GSE186143 is used only to freeze marker directions and the published-marker
# comparator. No feature selection, model fitting, or tuning uses either
# external cohort.

workspace <- normalizePath(
  Sys.getenv("IMMUNOCOMPASS_WORKSPACE", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
default_tool_root <- if (dir.exists(file.path(workspace, "Result"))) {
  workspace
} else {
  file.path(workspace, "返修", "免疫平衡评分工具")
}
tool_root <- gsub(
  "\\\\", "/",
  path.expand(Sys.getenv("IBI_TOOL_ROOT", unset = default_tool_root)),
  fixed = TRUE
)
if (!dir.exists(tool_root)) stop("Tool root does not exist: ", tool_root)
default_supplementary_data <- if (
  dir.exists(file.path(tool_root, "revision_supp_data"))
) {
  file.path(tool_root, "revision_supp_data")
} else {
  file.path(workspace, "返修", "补充数据")
}
supplementary_data <- gsub(
  "\\\\", "/",
  path.expand(Sys.getenv(
    "IRAE_SUPPLEMENTARY_DATA", unset = default_supplementary_data
  )),
  fixed = TRUE
)
if (!dir.exists(supplementary_data)) {
  stop("Supplementary-data directory does not exist: ", supplementary_data)
}
output_dir <- file.path(
  tool_root, "Result",
  "revision_irAE_external_validation_GSE319496_GSE307931_20261026"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

training_expression_file <- file.path(
  tool_root, "Result", "validation", "GSE186143_bulk_TPM.txt"
)
training_metadata_file <- file.path(
  tool_root, "Result", "validation", "GSE186143_sam_info.txt"
)
cas_statistics_file <- file.path(
  tool_root, "Result", "revision_BAS_CAS_IAS_sensitivity_top390",
  "standardized_BAS_CAS_IAS_gene_level.csv"
)
gse319496_counts_file <- file.path(
  supplementary_data, "GSE319496_GEO_raw_counts_SampleID.csv"
)
gse319496_metadata_file <- file.path(
  supplementary_data, "GSE319496_series_matrix.txt"
)
gse307931_metadata_file <- file.path(
  supplementary_data, "GSE307931_series_matrix.txt"
)
gse307931_raw_dir <- file.path(supplementary_data, "GSE307931_RAW")

required_files <- c(
  training_expression_file, training_metadata_file, cas_statistics_file,
  gse319496_counts_file, gse319496_metadata_file,
  gse307931_metadata_file
)
if (!all(file.exists(required_files))) {
  stop(
    "Missing input file(s): ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}
if (!dir.exists(gse307931_raw_dir)) {
  stop("Missing GSE307931 per-sample expression directory: ", gse307931_raw_dir)
}

master_seed <- 20261026L
n_bootstraps <- 2000L
cas_genes <- c("DTWD1", "BBS7", "CARF", "RABEP1", "DCK")
published_marker_candidates <- unique(c(
  "IL1B", "IL2", "CSF2", "CXCL8", "CSF3", "CCL2", "LEP", "CCL5",
  "IL6", "CRP", "IL17A", "ANGPT1", "CD40LG", "CXCL2", "CCL20",
  "CCL23", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "MICA", "CD163",
  "CX3CL1", "FGF2", "IFNA2", "IL12A", "IL1A", "IL1RN", "IL13",
  "IFNG", "CD8A", "GZMB", "PRF1", "STAT1", "IDO1"
))

strip_quotes <- function(x) {
  sub('^"', "", sub('"$', "", x))
}

parse_geo_metadata <- function(path) {
  lines <- readLines(path, warn = FALSE)
  read_geo_row <- function(prefix) {
    hit <- lines[startsWith(lines, prefix)]
    if (length(hit) != 1L) {
      stop("Expected one GEO row for ", prefix, " in ", path)
    }
    strip_quotes(strsplit(hit, "\t", fixed = TRUE)[[1L]][-1L])
  }

  accessions <- read_geo_row("!Sample_geo_accession")
  titles <- read_geo_row("!Sample_title")
  if (length(accessions) != length(titles)) {
    stop("GEO accession/title count mismatch in ", path)
  }
  metadata <- data.frame(
    geo_accession = accessions,
    sample_title = titles,
    stringsAsFactors = FALSE
  )

  characteristic_lines <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  for (line in characteristic_lines) {
    values <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1L]][-1L])
    if (length(values) != nrow(metadata)) {
      stop("GEO characteristic count mismatch in ", path)
    }
    key <- tolower(trimws(sub(":.*$", "", values[1L])))
    key <- gsub("[^a-z0-9]+", "_", key)
    key <- gsub("^_|_$", "", key)
    value <- trimws(sub("^[^:]+:\\s*", "", values))
    if (key %in% names(metadata)) {
      key <- make.unique(c(names(metadata), key))[length(names(metadata)) + 1L]
    }
    metadata[[key]] <- value
  }
  metadata
}

fixed_auc <- function(outcome, score) {
  if (length(unique(outcome)) != 2L || any(!is.finite(score))) {
    return(NA_real_)
  }
  as.numeric(auc(roc(
    outcome, score, levels = c(0, 1), direction = "<", quiet = TRUE
  )))
}

average_precision <- function(outcome, score) {
  ordering <- order(score, decreasing = TRUE)
  ordered_outcome <- outcome[ordering]
  if (sum(ordered_outcome == 1L) == 0L) return(NA_real_)
  cumulative_precision <- cumsum(ordered_outcome) / seq_along(ordered_outcome)
  mean(cumulative_precision[ordered_outcome == 1L])
}

standardize_rows <- function(expression_matrix) {
  means <- rowMeans(expression_matrix, na.rm = TRUE)
  standard_deviations <- apply(expression_matrix, 1L, sd, na.rm = TRUE)
  keep <- is.finite(standard_deviations) & standard_deviations > 0
  result <- sweep(expression_matrix[keep, , drop = FALSE], 1L, means[keep], "-")
  sweep(result, 1L, standard_deviations[keep], "/")
}

make_composite_score <- function(expression_matrix, genes, directions) {
  usable_genes <- intersect(genes, rownames(expression_matrix))
  if (length(usable_genes) == 0L) stop("No panel genes are evaluable.")
  z_matrix <- standardize_rows(expression_matrix[usable_genes, , drop = FALSE])
  usable_genes <- rownames(z_matrix)
  if (length(usable_genes) == 0L) stop("All evaluable panel genes are constant.")
  signed_z <- sweep(
    z_matrix, 1L, as.numeric(directions[usable_genes]), "*"
  )
  list(
    score = colMeans(signed_z),
    usable_genes = usable_genes,
    missing_genes = setdiff(genes, usable_genes)
  )
}

bootstrap_panel_metrics <- function(
  outcome, cas_score, published_score, n_iterations, seed
) {
  positive_index <- which(outcome == 1L)
  negative_index <- which(outcome == 0L)
  output <- vector("list", n_iterations * 2L)
  differences <- numeric(n_iterations)
  counter <- 0L
  set.seed(seed)
  for (iteration in seq_len(n_iterations)) {
    sampled_index <- c(
      sample(negative_index, length(negative_index), replace = TRUE),
      sample(positive_index, length(positive_index), replace = TRUE)
    )
    sampled_outcome <- outcome[sampled_index]
    cas_auc <- fixed_auc(sampled_outcome, cas_score[sampled_index])
    published_auc <- fixed_auc(
      sampled_outcome, published_score[sampled_index]
    )
    differences[iteration] <- cas_auc - published_auc
    for (method in c("Fixed CAS panel", "Fixed published-marker panel")) {
      score <- if (method == "Fixed CAS panel") cas_score else published_score
      counter <- counter + 1L
      output[[counter]] <- data.frame(
        iteration = iteration,
        method = method,
        metric = c("AUC", "Average_precision"),
        value = c(
          fixed_auc(sampled_outcome, score[sampled_index]),
          average_precision(sampled_outcome, score[sampled_index])
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  list(metrics = do.call(rbind, output), auc_difference = differences)
}

bootstrap_single_score <- function(outcome, score, n_iterations, seed) {
  positive_index <- which(outcome == 1L)
  negative_index <- which(outcome == 0L)
  values <- numeric(n_iterations)
  set.seed(seed)
  for (iteration in seq_len(n_iterations)) {
    sampled_index <- c(
      sample(negative_index, length(negative_index), replace = TRUE),
      sample(positive_index, length(positive_index), replace = TRUE)
    )
    values[iteration] <- fixed_auc(outcome[sampled_index], score[sampled_index])
  }
  values
}

training_expression <- as.matrix(read.delim(
  training_expression_file, row.names = 1L, check.names = FALSE
))
storage.mode(training_expression) <- "double"
training_metadata <- read.delim(
  training_metadata_file, row.names = 1L, check.names = FALSE
)
training_grade_text <- as.character(
  training_metadata[colnames(training_expression), "Grade"]
)
training_grade <- suppressWarnings(as.numeric(sub("^Grade", "", training_grade_text)))
training_keep <- is.finite(training_grade)
training_expression <- training_expression[, training_keep, drop = FALSE]
training_grade <- training_grade[training_keep]
training_outcome <- as.integer(training_grade >= 2L)
training_log_expression <- log2(training_expression + 1)

cas_statistics <- read.csv(
  cas_statistics_file, check.names = FALSE, stringsAsFactors = FALSE
)
cas_beta <- cas_statistics$beta_CAS[match(cas_genes, cas_statistics$gene)]
if (any(!is.finite(cas_beta)) || any(cas_beta == 0)) {
  stop("A CAS direction could not be frozen from the GSE186143 analysis.")
}
cas_directions <- setNames(sign(cas_beta), cas_genes)

available_published <- intersect(
  published_marker_candidates, rownames(training_log_expression)
)
published_training_results <- do.call(rbind, lapply(
  available_published,
  function(gene) {
    expression_value <- as.numeric(training_log_expression[gene, ])
    direction <- sign(
      mean(expression_value[training_outcome == 1L]) -
        mean(expression_value[training_outcome == 0L])
    )
    if (!is.finite(direction) || direction == 0) direction <- 1
    score <- direction * expression_value
    data.frame(
      gene = gene,
      direction = direction,
      training_apparent_auc = fixed_auc(training_outcome, score),
      stringsAsFactors = FALSE
    )
  }
))
published_training_results <- published_training_results[
  order(
    published_training_results$training_apparent_auc,
    decreasing = TRUE, na.last = NA
  ),
]
published_genes <- head(published_training_results$gene, length(cas_genes))
published_directions <- setNames(
  published_training_results$direction[
    match(published_genes, published_training_results$gene)
  ],
  published_genes
)

write.csv(
  data.frame(
    panel = "CAS", gene = cas_genes, direction = cas_directions,
    selection_source = "Manuscript-reported fixed panel",
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "fixed_CAS_panel.csv"), row.names = FALSE
)
write.csv(
  published_training_results,
  file.path(output_dir, "published_marker_training_ranking.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    panel = "Published marker", gene = published_genes,
    direction = published_directions,
    selection_source = "Frozen in GSE186143 before external validation",
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "fixed_published_marker_panel.csv"), row.names = FALSE
)

# GSE319496: raw counts -> TMM-normalized log2 CPM.
gse319496_metadata <- parse_geo_metadata(gse319496_metadata_file)
if (!"irae_status" %in% names(gse319496_metadata)) {
  stop("GSE319496 irAE status was not found in GEO metadata.")
}
gse319496_counts <- read.csv(
  gse319496_counts_file, row.names = 1L, check.names = FALSE
)
gse319496_counts <- as.matrix(gse319496_counts)
storage.mode(gse319496_counts) <- "double"
if (anyDuplicated(rownames(gse319496_counts))) {
  gse319496_counts <- rowsum(gse319496_counts, rownames(gse319496_counts))
}
gse319496_features_total <- nrow(gse319496_counts)
gse319496_complete_rows <- rowSums(is.na(gse319496_counts)) == 0L
gse319496_counts <- gse319496_counts[
  gse319496_complete_rows, , drop = FALSE
]
gse319496_features_excluded_missing <- sum(!gse319496_complete_rows)
gse319496_metadata$expression_id <- sub(
  ".*(Sample[0-9]+)$", "\\1", gse319496_metadata$sample_title
)
gse319496_metadata <- gse319496_metadata[
  match(colnames(gse319496_counts), gse319496_metadata$expression_id),
]
if (anyNA(gse319496_metadata$sample_title)) {
  stop("GSE319496 expression columns could not be matched to GEO metadata.")
}
gse319496_outcome <- as.integer(tolower(gse319496_metadata$irae_status) == "yes")
gse319496_dge <- DGEList(counts = gse319496_counts)
gse319496_dge <- calcNormFactors(gse319496_dge, method = "TMM")
gse319496_expression <- cpm(gse319496_dge, log = TRUE, prior.count = 0.5)

# GSE307931: assemble the deposited per-sample TPM files.
gse307931_metadata <- parse_geo_metadata(gse307931_metadata_file)
if (!"irae" %in% names(gse307931_metadata)) {
  stop("GSE307931 irAE status was not found in GEO metadata.")
}
gse307931_files <- list.files(
  gse307931_raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE
)
if (length(gse307931_files) != nrow(gse307931_metadata)) {
  stop(
    "GSE307931 file/sample mismatch: ", length(gse307931_files),
    " files versus ", nrow(gse307931_metadata), " metadata rows."
  )
}
target_genes <- unique(c(cas_genes, published_genes))
gse307931_features_deposited <- nrow(read.delim(
  gzfile(gse307931_files[1L]), check.names = FALSE
))
gse307931_parts <- lapply(gse307931_files, function(path) {
  table <- read.delim(gzfile(path), check.names = FALSE)
  gene <- as.character(table[[1L]])
  value <- as.numeric(table[[2L]])
  names(value) <- gene
  sample_geo <- sub("_.*$", "", basename(path))
  sample_title <- gse307931_metadata$sample_title[
    match(sample_geo, gse307931_metadata$geo_accession)
  ]
  if (length(sample_title) != 1L || is.na(sample_title)) {
    stop("Could not match GSE307931 file to metadata: ", basename(path))
  }
  output <- value[target_genes]
  names(output) <- target_genes
  list(sample_title = sample_title, expression = output)
})
gse307931_expression <- do.call(
  cbind, lapply(gse307931_parts, function(x) x$expression)
)
colnames(gse307931_expression) <- vapply(
  gse307931_parts, function(x) x$sample_title, character(1L)
)
gse307931_metadata <- gse307931_metadata[
  match(colnames(gse307931_expression), gse307931_metadata$sample_title),
]
gse307931_outcome <- as.integer(tolower(gse307931_metadata$irae) == "yes")
gse307931_expression <- log2(gse307931_expression + 1)

cohorts <- list(
  GSE319496 = list(
    expression = gse319496_expression,
    outcome = gse319496_outcome,
    sample_id = gse319496_metadata$geo_accession,
    endpoint = "CTCAE Grade 2-4 vs Grade 0-1",
    tissue = "Pretreatment whole blood",
    platform = "Bulk RNA-seq raw counts; TMM log2 CPM",
    expression_features_deposited = gse319496_features_total,
    expression_features_available_after_qc = nrow(gse319496_expression),
    target_panel_features_loaded = sum(
      target_genes %in% rownames(gse319496_expression)
    ),
    expression_features_excluded_missing = gse319496_features_excluded_missing
  ),
  GSE307931 = list(
    expression = gse307931_expression,
    outcome = gse307931_outcome,
    sample_id = gse307931_metadata$geo_accession,
    endpoint = "Any irAE vs no irAE",
    tissue = "Pretreatment whole blood",
    platform = "Bulk RNA-seq deposited TPM; log2(TPM + 1)",
    expression_features_deposited = gse307931_features_deposited,
    expression_features_available_after_qc = gse307931_features_deposited,
    target_panel_features_loaded = nrow(gse307931_expression),
    expression_features_excluded_missing = 0L
  )
)

cohort_summaries <- list()
coverage_outputs <- list()
sample_outputs <- list()
panel_metric_outputs <- list()
panel_bootstrap_outputs <- list()
panel_difference_outputs <- list()
individual_metric_outputs <- list()
individual_bootstrap_outputs <- list()

for (cohort_index in seq_along(cohorts)) {
  cohort_name <- names(cohorts)[cohort_index]
  cohort <- cohorts[[cohort_name]]
  outcome <- cohort$outcome
  expression_matrix <- cohort$expression
  if (length(outcome) != ncol(expression_matrix)) {
    stop("Outcome/expression mismatch in ", cohort_name)
  }
  if (!identical(sort(unique(outcome)), c(0L, 1L))) {
    stop("Unexpected outcome coding in ", cohort_name)
  }

  cas_panel <- make_composite_score(
    expression_matrix, cas_genes, cas_directions
  )
  published_panel <- make_composite_score(
    expression_matrix, published_genes, published_directions
  )

  cohort_summaries[[cohort_name]] <- data.frame(
    cohort = cohort_name,
    endpoint = cohort$endpoint,
    tissue = cohort$tissue,
    platform = cohort$platform,
    n_total = length(outcome),
    n_event = sum(outcome == 1L),
    n_nonevent = sum(outcome == 0L),
    cas_genes_requested = length(cas_genes),
    cas_genes_evaluable = length(cas_panel$usable_genes),
    published_genes_requested = length(published_genes),
    published_genes_evaluable = length(published_panel$usable_genes),
    expression_features_deposited = cohort$expression_features_deposited,
    expression_features_available_after_qc =
      cohort$expression_features_available_after_qc,
    target_panel_features_loaded = cohort$target_panel_features_loaded,
    expression_features_excluded_missing =
      cohort$expression_features_excluded_missing,
    stringsAsFactors = FALSE
  )

  coverage_outputs[[paste0(cohort_name, "_CAS")]] <- data.frame(
    cohort = cohort_name,
    panel = "Fixed CAS panel",
    gene = cas_genes,
    evaluable = cas_genes %in% cas_panel$usable_genes,
    reason = ifelse(
      cas_genes %in% cas_panel$usable_genes,
      "Evaluable", "Absent or invariant"
    ),
    stringsAsFactors = FALSE
  )
  coverage_outputs[[paste0(cohort_name, "_published")]] <- data.frame(
    cohort = cohort_name,
    panel = "Fixed published-marker panel",
    gene = published_genes,
    evaluable = published_genes %in% published_panel$usable_genes,
    reason = ifelse(
      published_genes %in% published_panel$usable_genes,
      "Evaluable", "Absent or invariant"
    ),
    stringsAsFactors = FALSE
  )

  sample_outputs[[cohort_name]] <- data.frame(
    cohort = cohort_name,
    sample_id = cohort$sample_id,
    outcome = outcome,
    outcome_label = ifelse(outcome == 1L, "irAE", "No irAE"),
    cas_score = cas_panel$score,
    published_marker_score = published_panel$score,
    stringsAsFactors = FALSE
  )

  central_panel_metrics <- do.call(rbind, lapply(
    c("Fixed CAS panel", "Fixed published-marker panel"),
    function(method) {
      score <- if (method == "Fixed CAS panel") {
        cas_panel$score
      } else {
        published_panel$score
      }
      data.frame(
        cohort = cohort_name,
        method = method,
        metric = c("AUC", "Average_precision"),
        estimate = c(
          fixed_auc(outcome, score),
          average_precision(outcome, score)
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  bootstrap_panel <- bootstrap_panel_metrics(
    outcome, cas_panel$score, published_panel$score,
    n_bootstraps, master_seed + cohort_index * 100000L
  )
  bootstrap_panel$metrics$cohort <- cohort_name
  panel_bootstrap_outputs[[cohort_name]] <- bootstrap_panel$metrics
  panel_metric_outputs[[cohort_name]] <- do.call(rbind, lapply(
    split(
      bootstrap_panel$metrics,
      list(bootstrap_panel$metrics$method, bootstrap_panel$metrics$metric),
      drop = TRUE
    ),
    function(values) {
      estimate <- central_panel_metrics$estimate[
        central_panel_metrics$method == unique(values$method) &
          central_panel_metrics$metric == unique(values$metric)
      ]
      data.frame(
        cohort = cohort_name,
        method = unique(values$method),
        metric = unique(values$metric),
        estimate = estimate,
        ci_lower = unname(quantile(values$value, 0.025, na.rm = TRUE)),
        ci_upper = unname(quantile(values$value, 0.975, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
  ))
  panel_difference_outputs[[cohort_name]] <- data.frame(
    cohort = cohort_name,
    contrast = "Fixed CAS panel - Fixed published-marker panel",
    estimate = central_panel_metrics$estimate[
      central_panel_metrics$method == "Fixed CAS panel" &
        central_panel_metrics$metric == "AUC"
    ] - central_panel_metrics$estimate[
      central_panel_metrics$method == "Fixed published-marker panel" &
        central_panel_metrics$metric == "AUC"
    ],
    ci_lower = unname(quantile(
      bootstrap_panel$auc_difference, 0.025, na.rm = TRUE
    )),
    ci_upper = unname(quantile(
      bootstrap_panel$auc_difference, 0.975, na.rm = TRUE
    )),
    stringsAsFactors = FALSE
  )

  individual_genes <- intersect(cas_genes, rownames(expression_matrix))
  for (gene_index in seq_along(individual_genes)) {
    gene <- individual_genes[gene_index]
    raw_score <- as.numeric(expression_matrix[gene, ])
    if (!is.finite(sd(raw_score)) || sd(raw_score) == 0) next
    score <- as.numeric(cas_directions[gene]) * raw_score
    bootstrap_auc <- bootstrap_single_score(
      outcome, score, n_bootstraps,
      master_seed + cohort_index * 100000L + gene_index * 1000L
    )
    individual_metric_outputs[[paste(cohort_name, gene, sep = "_")]] <- data.frame(
      cohort = cohort_name,
      gene = gene,
      direction = as.numeric(cas_directions[gene]),
      estimate = fixed_auc(outcome, score),
      ci_lower = unname(quantile(bootstrap_auc, 0.025, na.rm = TRUE)),
      ci_upper = unname(quantile(bootstrap_auc, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    individual_bootstrap_outputs[[paste(cohort_name, gene, sep = "_")]] <- data.frame(
      cohort = cohort_name,
      gene = gene,
      iteration = seq_len(n_bootstraps),
      AUC = bootstrap_auc,
      stringsAsFactors = FALSE
    )
  }
}

cohort_summary <- do.call(rbind, cohort_summaries)
coverage <- do.call(rbind, coverage_outputs)
sample_scores <- do.call(rbind, sample_outputs)
panel_metrics <- do.call(rbind, panel_metric_outputs)
panel_bootstrap <- do.call(rbind, panel_bootstrap_outputs)
panel_auc_differences <- do.call(rbind, panel_difference_outputs)
individual_cas_metrics <- do.call(rbind, individual_metric_outputs)
individual_cas_bootstrap <- do.call(rbind, individual_bootstrap_outputs)

write.csv(
  cohort_summary, file.path(output_dir, "external_cohort_summary.csv"),
  row.names = FALSE
)
write.csv(
  coverage, file.path(output_dir, "external_panel_gene_coverage.csv"),
  row.names = FALSE
)
write.csv(
  sample_scores, file.path(output_dir, "external_sample_scores.csv"),
  row.names = FALSE
)
write.csv(
  panel_metrics, file.path(output_dir, "external_panel_performance.csv"),
  row.names = FALSE
)
write.csv(
  panel_bootstrap,
  file.path(output_dir, "external_panel_bootstrap_metrics.csv"),
  row.names = FALSE
)
write.csv(
  panel_auc_differences,
  file.path(output_dir, "external_panel_auc_differences.csv"),
  row.names = FALSE
)
write.csv(
  individual_cas_metrics,
  file.path(output_dir, "external_individual_CAS_gene_AUC.csv"),
  row.names = FALSE
)
write.csv(
  individual_cas_bootstrap,
  file.path(output_dir, "external_individual_CAS_gene_bootstrap.csv"),
  row.names = FALSE
)

# Figure contract: independent, leakage-free external discrimination with
# endpoint-specific uncertainty. The external cohorts are never pooled because
# their irAE definitions differ.
theme_external <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 7.5),
    plot.tag = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

distribution_data <- sample_scores
distribution_data$cohort_label <- factor(
  distribution_data$cohort,
  levels = c("GSE319496", "GSE307931"),
  labels = c(
    "GSE319496: Grade 0-1 vs Grade 2-4",
    "GSE307931: no irAE vs any irAE"
  )
)
distribution_data$outcome_label <- factor(
  distribution_data$outcome_label, levels = c("No irAE", "irAE")
)
distribution_counts <- aggregate(
  sample_id ~ cohort_label + outcome_label,
  data = distribution_data, FUN = length
)
distribution_counts$display_label <- paste0(
  distribution_counts$outcome_label, "\n(n=", distribution_counts$sample_id, ")"
)
distribution_data$display_label <- distribution_counts$display_label[
  match(
    paste(distribution_data$cohort_label, distribution_data$outcome_label),
    paste(distribution_counts$cohort_label, distribution_counts$outcome_label)
  )
]
distribution_data$display_label <- factor(
  distribution_data$display_label,
  levels = unique(distribution_counts$display_label)
)

p_a <- ggplot(
  distribution_data,
  aes(display_label, cas_score, fill = outcome_label, colour = outcome_label)
) +
  geom_violin(trim = FALSE, alpha = 0.20, linewidth = 0.35) +
  geom_boxplot(
    width = 0.25, outlier.shape = NA, alpha = 0.65, linewidth = 0.35
  ) +
  geom_jitter(width = 0.10, size = 0.8, alpha = 0.65, stroke = 0) +
  facet_wrap(~ cohort_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("No irAE" = "#B8BDC5", "irAE" = "#C95C4A")) +
  scale_colour_manual(values = c("No irAE" = "#626970", "irAE" = "#943B2F")) +
  labs(
    title = "Direction-locked CAS score in independent cohorts",
    x = NULL, y = "Signed standardized CAS score"
  ) +
  theme_external +
  theme(legend.position = "none")

panel_auc_plot <- panel_metrics[panel_metrics$metric == "AUC", ]
panel_auc_plot$cohort_label <- factor(
  panel_auc_plot$cohort,
  levels = c("GSE319496", "GSE307931"),
  labels = c("GSE319496", "GSE307931")
)
panel_auc_plot$method <- factor(
  panel_auc_plot$method,
  levels = c("Fixed published-marker panel", "Fixed CAS panel")
)
p_b <- ggplot(
  panel_auc_plot,
  aes(estimate, method, colour = method)
) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#A5A5A5") +
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    orientation = "y", width = 0.13, linewidth = 0.45
  ) +
  geom_point(size = 2) +
  facet_wrap(~ cohort_label, ncol = 1) +
  scale_colour_manual(values = c(
    "Fixed CAS panel" = "#276A9B",
    "Fixed published-marker panel" = "#777D84"
  )) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    title = "External panel discrimination",
    x = "ROC AUC (95% bootstrap CI)", y = NULL
  ) +
  theme_external +
  theme(legend.position = "none")

individual_plot <- individual_cas_metrics
individual_plot$cohort_label <- factor(
  individual_plot$cohort,
  levels = c("GSE319496", "GSE307931")
)
individual_plot$gene <- factor(individual_plot$gene, levels = rev(cas_genes))
p_c <- ggplot(
  individual_plot,
  aes(estimate, gene, colour = cohort_label, shape = cohort_label)
) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#A5A5A5") +
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    orientation = "y", width = 0.13, linewidth = 0.4,
    position = position_dodge(width = 0.35)
  ) +
  geom_point(size = 1.8, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c("GSE319496" = "#276A9B", "GSE307931" = "#C95C4A")) +
  scale_shape_manual(values = c("GSE319496" = 16, "GSE307931" = 17)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    title = "Individual fixed CAS biomarkers",
    x = "Direction-locked ROC AUC (95% bootstrap CI)", y = NULL
  ) +
  theme_external +
  theme(legend.position = "bottom")

final_figure <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.05, 1), widths = c(0.9, 1.2)) +
  plot_annotation(tag_levels = "A")

figure_file <- file.path(
  output_dir, "irAE_independent_external_validation.pdf"
)
grDevices::cairo_pdf(
  figure_file, width = 183 / 25.4, height = 132 / 25.4,
  family = "Arial", onefile = TRUE
)
print(final_figure)
dev.off()

qa_preview_file <- Sys.getenv("IRAE_QA_PREVIEW", unset = "")
if (nzchar(qa_preview_file)) {
  grDevices::png(
    qa_preview_file, width = 183 / 25.4, height = 132 / 25.4,
    units = "in", res = 180, type = "cairo"
  )
  print(final_figure)
  dev.off()
}

provenance <- data.frame(
  item = c(
    "Discovery/internal cohort", "Fixed CAS panel",
    "External cohort 1", "External cohort 1 endpoint",
    "External cohort 2", "External cohort 2 endpoint",
    "External feature selection", "External model fitting",
    "Score construction", "Missing gene handling",
    "Confidence intervals", "Bootstrap seed", "Figure format"
  ),
  value = c(
    "GSE186143 (60 patients with complete CTCAE grade)",
    paste(cas_genes, collapse = ";"),
    "GSE319496", "CTCAE Grade 2-4 versus Grade 0-1",
    "GSE307931", "Any irAE versus no irAE",
    "None", "None",
    paste(
      "Outcome-blind within-cohort z-scores; training-direction-locked",
      "genes averaged with equal weights"
    ),
    "No imputation; absent or invariant genes excluded and reported",
    paste0(n_bootstraps, " stratified patient-level bootstrap resamples"),
    master_seed, "PDF only"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  provenance, file.path(output_dir, "analysis_provenance.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    cohort_summary = cohort_summary,
    coverage = coverage,
    sample_scores = sample_scores,
    panel_metrics = panel_metrics,
    panel_auc_differences = panel_auc_differences,
    individual_cas_metrics = individual_cas_metrics,
    fixed_cas_panel = data.frame(
      gene = cas_genes, direction = cas_directions
    ),
    fixed_published_panel = data.frame(
      gene = published_genes, direction = published_directions
    ),
    provenance = provenance
  ),
  file.path(output_dir, "external_validation_objects.rds")
)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))

message("Independent external validation complete: ", output_dir)
