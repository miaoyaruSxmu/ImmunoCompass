# Avoid false-positive global-variable notes for lazy-loaded package data.
utils::globalVariables(c("gene_pairs", "gene_pair_modes", "svr_model"))

#' Convert an expression input to a numeric gene-by-sample matrix.
#'
#' @param expr A gene-by-sample matrix, data frame, or tab-delimited file path.
#'   Gene identifiers must be row names or the first column of a file.
#' @return A numeric matrix with genes in rows and samples in columns.
#' @keywords internal
as_expression_matrix <- function(expr) {
  if (is.character(expr) && length(expr) == 1L) {
    expr <- utils::read.delim(
      expr, row.names = 1L, check.names = FALSE, quote = "",
      comment.char = "", stringsAsFactors = FALSE
    )
  }
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr)) stop("expr must be a matrix, data frame, or file path")
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("expr must have gene row names and sample column names")
  }
  if (anyDuplicated(rownames(expr))) {
    stop("expr contains duplicated gene identifiers; aggregate them before scoring")
  }
  suppressWarnings(storage.mode(expr) <- "double")
  expr
}

#' Parse the two gene symbols in a locked gene-pair identifier.
#' @keywords internal
parse_gene_pair <- function(pair) {
  parts <- strsplit(pair, "_vs_", fixed = TRUE)[[1L]]
  if (length(parts) != 2L || any(!nzchar(parts))) {
    stop("Invalid gene-pair identifier: ", pair)
  }
  parts
}

#' Build the fixed binary gene-pair feature matrix.
#'
#' A pair is encoded as 1 when GeneA expression is greater than GeneB
#' expression and 0 otherwise, including ties. Missing or non-finite pair
#' values are assigned the corresponding binary ordering mode estimated from
#' the 807 model-development samples by default.
#'
#' @param expr A gene-by-sample matrix, data frame, or tab-delimited file path.
#' @param pairs Character vector of `GeneA_vs_GeneB` identifiers. Defaults to
#'   the bundled 390-pair panel.
#' @param missing How missing pair values are handled: `"development_mode"`
#'   uses the bundled pair-specific development modes, `"zero"` assigns 0,
#'   and `"error"` stops.
#' @param missing_state Optional numeric replacement value for missing pairs.
#'   A scalar is recycled across pairs, or a named or ordered vector can be
#'   supplied for pair-specific states. When `NULL`, the bundled development
#'   modes are used for `missing = "development_mode"`.
#' @return A sample-by-pair numeric matrix.
#' @export
build_pair_features <- function(
    expr, pairs = gene_pairs,
    missing = c("development_mode", "zero", "error"),
    missing_state = NULL) {
  missing <- match.arg(missing)
  expr <- as_expression_matrix(expr)
  pairs <- as.character(pairs)
  if (!length(pairs)) stop("pairs must contain at least one gene pair")

  if (missing == "development_mode" && is.null(missing_state)) {
    mode_index <- match(pairs, names(gene_pair_modes))
    if (anyNA(mode_index)) {
      stop("Bundled development modes are unavailable for one or more requested pairs")
    }
    missing_state <- unname(gene_pair_modes[mode_index])
  }
  if (is.null(missing_state)) missing_state <- 0
  if (length(missing_state) == 1L) missing_state <- rep(missing_state, length(pairs))
  if (length(missing_state) != length(pairs)) {
    stop("missing_state must be scalar or have one value per pair")
  }
  if (missing == "zero") missing_state[] <- 0
  if (missing == "error") missing_state[] <- NA_real_

  out <- matrix(
    NA_real_, nrow = ncol(expr), ncol = length(pairs),
    dimnames = list(colnames(expr), pairs)
  )
  for (j in seq_along(pairs)) {
    genes <- parse_gene_pair(pairs[j])
    if (all(genes %in% rownames(expr))) {
      a <- as.numeric(expr[genes[1L], ])
      b <- as.numeric(expr[genes[2L], ])
      observed <- is.finite(a) & is.finite(b)
      out[observed, j] <- as.numeric(a[observed] > b[observed])
      if (any(!observed)) out[!observed, j] <- missing_state[j]
    } else {
      out[, j] <- missing_state[j]
    }
    if (missing == "error" && any(!is.finite(out[, j]))) {
      stop("Missing expression prevents computation of pair ", pairs[j])
    }
  }
  out
}

#' Calculate the Immune Balance Index (IBI).
#'
#' @param expr A gene-by-sample expression matrix, data frame, or file path.
#' @param model A fitted e1071 SVR object. Defaults to the bundled model.
#' @param pairs The fixed gene-pair panel. Defaults to the bundled 390 pairs.
#' @param missing Missing-pair handling passed to [build_pair_features()].
#' @param missing_state Replacement value or pair-specific development modes.
#' @param return_features If `TRUE`, return IBI and the generated feature matrix.
#' @return A named numeric IBI vector, or a list containing `IBI` and `features`.
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @export
calculate_ibi <- function(
    expr, model = svr_model, pairs = gene_pairs,
    missing = c("development_mode", "zero", "error"),
    missing_state = NULL, return_features = FALSE) {
  features <- build_pair_features(
    expr, pairs = pairs, missing = missing, missing_state = missing_state
  )
  prediction <- as.numeric(predict(model, newdata = as.data.frame(features)))
  names(prediction) <- rownames(features)
  if (return_features) return(list(IBI = prediction, features = features))
  prediction
}

#' Standardize a numeric vector without silently accepting constants.
#' @keywords internal
z_score <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("standardization variables must be finite")
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) stop("cannot standardize a constant variable")
  (x - mean(x)) / s
}

#' Decompose gene-level expression effects into BAS, CAS and IAS.
#'
#' The fitted model is `gene expression ~ standardized IBI * standardized
#' condition`. BAS is the absolute standardized IBI coefficient, CAS is the
#' absolute standardized condition coefficient, and IAS is the absolute
#' standardized interaction coefficient. The returned signed coefficients are
#' retained for biological interpretation.
#'
#' @param expression A gene-by-sample numeric expression matrix or data frame.
#'   Rows are genes and columns are samples. Expression is standardized within
#'   gene before fitting unless `standardize_expression = FALSE`.
#' @param ibi Numeric IBI values, one per sample.
#' @param condition Numeric clinical condition or grade, one per sample.
#'   Factors are converted to their underlying integer level order.
#' @param dominance Fold threshold for assigning a dominant component. The
#'   default is 1.5; a component must exceed both alternatives by this factor.
#' @param standardize_expression Whether to z-standardize each gene first.
#' @return A data frame with signed coefficients, absolute component scores,
#'   nominal P values, BH-adjusted FDR values and dominant-component labels.
#' @export
decompose_gene_effect <- function(
    expression, ibi, condition, dominance = 1.5,
    standardize_expression = TRUE) {
  expression <- as_expression_matrix(expression)
  if (is.factor(condition)) condition <- as.numeric(condition)
  ibi <- as.numeric(ibi)
  condition <- as.numeric(condition)
  if (length(ibi) != ncol(expression) || length(condition) != ncol(expression)) {
    stop("ibi and condition must have one value per expression column")
  }
  complete <- is.finite(ibi) & is.finite(condition)
  if (sum(complete) < 5L) stop("at least five complete samples are required")
  expression <- expression[, complete, drop = FALSE]
  ibi <- ibi[complete]
  condition <- condition[complete]
  ibi_z <- z_score(ibi)
  condition_z <- z_score(condition)

  gene_matrix <- expression
  if (standardize_expression) {
    gene_matrix <- t(apply(gene_matrix, 1L, function(y) {
      if (any(!is.finite(y)) || stats::sd(y) == 0) return(rep(NA_real_, length(y)))
      as.numeric(scale(y))
    }))
    rownames(gene_matrix) <- rownames(expression)
    colnames(gene_matrix) <- colnames(expression)
  }
  fit_one <- function(y) {
    ok <- is.finite(y)
    if (sum(ok) < 5L) return(rep(NA_real_, 6L))
    fit_data <- data.frame(
      y = y[ok], IBI = ibi_z[ok], Condition = condition_z[ok]
    )
    fit <- stats::lm(y ~ IBI * Condition, data = fit_data)
    sm <- summary(fit)$coefficients
    coef_names <- c("IBI", "Condition", "IBI:Condition")
    if (!all(coef_names %in% rownames(sm))) return(rep(NA_real_, 6L))
    c(sm[coef_names, 1L], sm[coef_names, 4L])
  }
  fitted <- t(apply(gene_matrix, 1L, fit_one))
  colnames(fitted) <- c("beta_BAS", "beta_CAS", "beta_IAS",
                        "p_BAS", "p_CAS", "p_IAS")
  result <- data.frame(gene = rownames(gene_matrix), fitted, row.names = NULL,
                       check.names = FALSE)
  result$BAS <- abs(result$beta_BAS)
  result$CAS <- abs(result$beta_CAS)
  result$IAS <- abs(result$beta_IAS)
  result$FDR_BAS <- stats::p.adjust(result$p_BAS, method = "BH")
  result$FDR_CAS <- stats::p.adjust(result$p_CAS, method = "BH")
  result$FDR_IAS <- stats::p.adjust(result$p_IAS, method = "BH")
  scores <- as.matrix(result[, c("BAS", "CAS", "IAS")])
  result$dominant_component <- "None"
  for (component in colnames(scores)) {
    others <- setdiff(colnames(scores), component)
    keep <- is.finite(scores[, component]) &
      scores[, component] > dominance * scores[, others[1L]] &
      scores[, component] > dominance * scores[, others[2L]]
    result$dominant_component[keep] <- component
  }
  result[, c("gene", "beta_BAS", "beta_CAS", "beta_IAS", "BAS", "CAS", "IAS",
             "p_BAS", "p_CAS", "p_IAS", "FDR_BAS", "FDR_CAS", "FDR_IAS",
             "dominant_component")]
}

#' Backward-compatible alias for [calculate_ibi()].
#' @param expr_loc Expression matrix, data frame or tab-delimited file path.
#' @return A named numeric IBI vector.
#' @export
ImmunoCompass <- function(expr_loc) calculate_ibi(expr_loc)

#' Development-set modal states for the locked gene-pair panel
#'
#' A named integer vector containing the modal binary ordering state for each
#' of the 390 locked gene pairs across the 807 model-development samples. These
#' values are used by default when a pair cannot be computed in new data.
#'
#' @format A named integer vector of length 390 with values 0 or 1.
#' @source The Table-S2-aligned 807-sample model-development matrix.
"gene_pair_modes"
