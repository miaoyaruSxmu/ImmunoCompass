library(ImmunoCompass)

data(example_expression)
data(example_metadata)
data(gene_pairs)
data(gene_pair_modes)
data(model_metadata)

features <- build_pair_features(example_expression)
stopifnot(
  identical(dim(features), c(6L, 390L)),
  all(features %in% c(0, 1)),
  identical(colnames(features), gene_pairs),
  identical(names(gene_pair_modes), gene_pairs),
  all(gene_pair_modes %in% c(0, 1)),
  identical(model_metadata$training_samples, 807L)
)

ibi <- calculate_ibi(example_expression)
stopifnot(length(ibi) == 6L, all(is.finite(ibi)))

effects <- decompose_gene_effect(
  example_expression, ibi, example_metadata$condition
)
stopifnot(nrow(effects) == nrow(example_expression))
