options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(
    paste(
      "Usage: Rscript new_data_prediction.r",
      "<gene_by_sample_expression.tsv> <output.csv>",
      "[development_mode|zero|error]"
    ),
    call. = FALSE
  )
}

input_file <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_file <- args[[2L]]
missing_rule <- if (length(args) == 3L) args[[3L]] else "development_mode"
if (!missing_rule %in% c("zero", "development_mode", "error")) {
  stop("Unknown missing-pair rule: ", missing_rule)
}

suppressPackageStartupMessages(library(ImmunoCompass))

expression <- read.delim(
  input_file,
  row.names = 1L,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
expression <- as.matrix(expression)
storage.mode(expression) <- "double"

result <- calculate_ibi(
  expression,
  missing = missing_rule,
  return_features = TRUE
)

pair_genes <- strsplit(gene_pairs, "_vs_", fixed = TRUE)
fully_evaluable <- vapply(
  pair_genes,
  function(genes) all(genes %in% rownames(expression)),
  logical(1)
)

output <- data.frame(
  sample_id = names(result$IBI),
  IBI = as.numeric(result$IBI),
  evaluable_gene_pairs = sum(fully_evaluable),
  total_gene_pairs = length(gene_pairs),
  evaluable_percentage = 100 * mean(fully_evaluable),
  missing_pair_rule = missing_rule,
  stringsAsFactors = FALSE
)
write.csv(output, output_file, row.names = FALSE)

cat("ImmunoCompass prediction complete\n")
cat("  samples:", nrow(output), "\n")
cat(
  "  evaluable pairs:", sum(fully_evaluable), "/", length(gene_pairs),
  sprintf("(%.1f%%)\n", 100 * mean(fully_evaluable))
)
cat("  output:", output_file, "\n")
