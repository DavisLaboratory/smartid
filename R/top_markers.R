#' compute group summarized score and order genes based on processed scores
#'
#' @inheritParams top_markers_abs
#' @inheritParams top_markers_glm
#' @param use.glm logical, if to use [stats::glm()] to compute group mean score,
#'     if TRUE, also compute mean score difference as output
#' @param ... params for [top_markers_abs()] or [top_markers_glm()]
#'
#' @return a tibble with feature names, group labels and ordered processed scores
#' @export
#'
#' @examples
#' data <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
#' top_markers_init(data, label = rep(c("A", "B"), 5))
top_markers_init <- function(data, label, n = 10,
                             use.glm = TRUE,
                             batch = NULL,
                             scale = TRUE,
                             use.mgm = TRUE,
                             softmax = TRUE,
                             ...) {
  if (use.glm == TRUE) {
    data <- top_markers_glm(
      data = data,
      label = label,
      batch = batch,
      n = n,
      scale = scale,
      use.mgm = use.mgm,
      softmax = softmax,
      ...
    )
  } else {
    data <- top_markers_abs(
      data = data,
      label = label,
      n = n,
      scale = scale,
      use.mgm = use.mgm,
      softmax = softmax,
      ...
    )
  }

  return(data)
}

#' calculate group median, MAD or mean score and order genes based on scores
#'
#' @inheritParams scale_mgm
#' @inheritParams top_markers_glm
#' @param method character, specify metric to compute, can be one of "median",
#'     "mad", "mean"
#'
#' @return a tibble with feature names, group labels and ordered processed scores
#' @export
#'
#' @examples
#' data <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
#' top_markers_abs(data, label = rep(c("A", "B"), 5))
top_markers_abs <- function(data, label, n = 10,
                            pooled.sd = FALSE,
                            method = c("median", "mad", "mean"),
                            scale = TRUE, use.mgm = TRUE,
                            softmax = TRUE,
                            tau = 1) {
  method <- match.arg(method)
  data <- apply_row_scaling(data, label, scale, use.mgm, pooled.sd)

  ## G x K aggregation goes straight to a long data.frame; skips the
  ## legacy `t() |> as.data.frame() |> summarise_all()` path which
  ## materialised a dense N x G frame (tens of GB on large inputs).
  agg  <- aggregate_rows_by_group(data, label, method)
  long <- long_format_from_group_matrix(agg)
  finalize_top_markers(long, n = n, softmax = softmax, tau = tau)
}

#' calculate group mean score using glm and order genes based on scores difference
#'
#' @inheritParams scale_mgm
#' @param data matrix, features in row and samples in column
#' @param n integer, number of returned top genes for each group
#' @param family family for glm, details in [stats::glm()]
#' @param batch a vector of batch labels, default NULL
#' @param scale logical, if to scale data by row
#' @param use.mgm logical, if to scale data using [scale_mgm()]
#' @param softmax logical, if to apply softmax transformation on output
#' @param tau numeric, hyper parameter for softmax
#'
#' @return a tibble with feature names, group labels and ordered processed scores
#' @export
#'
#' @examples
#' data <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
#' top_markers_glm(data, label = rep(c("A", "B"), 5))
top_markers_glm <- function(data, label, n = 10,
                            family = gaussian(), # score are continuous non-negative, can use gamma or inverse.gaussian, if continuous and unbounded use gaussian, if discrete use poisson, if binary or proportions between [0,1] or binary freq counts use binomial
                            batch = NULL,
                            scale = TRUE, use.mgm = TRUE,
                            pooled.sd = FALSE,
                            # log = TRUE,
                            softmax = TRUE,
                            tau = 1) {
  label <- factor(label)
  if (!is.null(batch)) batch <- factor(batch)

  data  <- apply_row_scaling(data, label, scale, use.mgm, pooled.sd)
  # ## log score
  # if(log == TRUE) {
  #   data <- log(data + 1e-8)
  # }
  betas <- fit_label_betas(data, label, batch, family)      # K x G
  betas <- betas_to_logfc_1v_max(betas)                     # K x G
  rownames(betas) <- levels(label)

  long <- data.frame(.dot = rownames(betas), betas,
                     check.names = FALSE, stringsAsFactors = FALSE) |>
    tidyr::pivot_longer(-`.dot`, names_to = "Genes",
                        values_to = "Scores") |>
    dplyr::group_by(`.dot`)

  finalize_top_markers(long, n = n, softmax = softmax, tau = tau)
}

## sigmoid: [0, 1], multi-label, no need to sum to 1
sigmoid <- function(x) {
  x <- x / max(abs(x))
  1 / (1 + exp(-x))
}

## softmax: [0, 1], one-label, multi-class, sum to 1
softmax <- function(x, tau = 1) {
  x <- x / tau
  exp(x) / sum(exp(x), na.rm = TRUE)
}

## tanh: [-1, 1], similar to sigmoid, no need to sum 1
tanh <- function(x) 2 / (1 + exp(-2 * x)) - 1

#################################################
# Internal helpers for top_markers_abs/glm
#################################################

## Row-wise z-score without densifying via `scale(t(data))`. Rows with
## zero or NA SD are collapsed to zero, matching the
## `data[is.na(data)] <- 0` guard from the legacy path.
row_scale_zmean <- function(data) {
  mu <- sparseMatrixStats::rowMeans2(data, na.rm = TRUE)
  sd <- sparseMatrixStats::rowSds(data,   na.rm = TRUE)
  sd[sd == 0 | is.na(sd)] <- 1
  out <- (data - mu) / sd
  out[is.na(out)] <- 0
  out
}

## Single entry point for the three `scale` / `use.mgm` branches shared
## between `top_markers_abs()` and `top_markers_glm()`.
apply_row_scaling <- function(data, label, scale, use.mgm, pooled.sd) {
  if (!isTRUE(scale)) return(data)
  if (isTRUE(use.mgm)) {
    return(scale_mgm(expr = data, label = label, pooled.sd = pooled.sd))
  }
  row_scale_zmean(data)
}

## G x K matrix of per-group row statistics. Uses `sparseMatrixStats` so
## the path stays sparse-friendly for dgCMatrix inputs.
aggregate_rows_by_group <- function(data, label, method) {
  groups <- unique(as.character(label))
  fn <- switch(method,
    mean   = sparseMatrixStats::rowMeans2,
    median = sparseMatrixStats::rowMedians,
    mad    = sparseMatrixStats::rowMads,
    stop("Unknown aggregation method: ", method, call. = FALSE)
  )
  label_ch <- as.character(label)
  agg <- vapply(groups, function(g)
    fn(data[, label_ch == g, drop = FALSE], na.rm = TRUE),
    numeric(nrow(data)))
  colnames(agg) <- groups
  rownames(agg) <- rownames(data)
  agg
}

## Transform G x K group-statistic matrix into the grouped long
## data.frame expected by downstream `markers_*()` consumers.
long_format_from_group_matrix <- function(agg) {
  out <- data.frame(
    .dot   = rep(colnames(agg), each = nrow(agg)),
    Genes  = rep(rownames(agg), times = ncol(agg)),
    Scores = as.vector(agg),
    stringsAsFactors = FALSE
  )
  dplyr::group_by(out, `.dot`)
}

## Softmax (optional) + top-n slice, shared finalisation.
finalize_top_markers <- function(long, n, softmax, tau) {
  if (isTRUE(softmax)) {
    # long <- dplyr::mutate(Scores = Scores / sd(Scores, na.rm = TRUE)) |> # norm by sd
    #  dplyr::mutate(Scores = sigmoid(Scores)) |> # sigmoid
    #  dplyr::mutate(Scores = tanh(Scores)) |> # tanh
    long <- dplyr::mutate(long, Scores = softmax(Scores, tau = tau))
  }
  dplyr::slice_max(long, Scores, n = n)
}

## Estimate per-gene label coefficients. Uses the closed-form solution
## whenever `family` is gaussian-identity and the design is full-rank;
## otherwise falls back to the legacy per-gene `glm()` loop so that
## users passing non-gaussian families keep the same behaviour.
fit_label_betas <- function(data, label, batch = NULL,
                            family = stats::gaussian()) {
  is_gauss_identity <- identical(family$family, "gaussian") &&
    identical(family$link, "identity")
  if (isTRUE(is_gauss_identity)) {
    res <- try(
      fit_label_betas_closed_form(data, label, batch),
      silent = TRUE
    )
    if (!inherits(res, "try-error") && !is.null(res)) return(res)
  }
  fit_label_betas_glm_loop(data, label, batch, family)
}

## Closed-form ordinary least squares for gaussian identity link.
## Returns K x G matrix of label coefficients; NULL signals a rank-
## deficient design so the caller can fall back to `glm()`.
fit_label_betas_closed_form <- function(data, label, batch = NULL) {
  X <- if (is.null(batch)) {
    Matrix::sparse.model.matrix(~ 0 + label)
  } else {
    Matrix::sparse.model.matrix(~ 0 + label + batch)
  }
  XtX <- Matrix::crossprod(X)
  if (Matrix::rankMatrix(XtX)[1] < ncol(X)) return(NULL)
  ## betas_all: K_total x G  (K_total = n label levels [+ batch levels])
  betas_all <- as.matrix(
    Matrix::solve(XtX, Matrix::crossprod(X, Matrix::t(data)))
  )
  rownames(betas_all) <- colnames(X)
  keep  <- grep("^label", rownames(betas_all))
  betas <- betas_all[keep, , drop = FALSE]
  rownames(betas) <- sub("^label", "", rownames(betas))
  betas
}

## Legacy per-gene glm loop; retained for non-gaussian families.
fit_label_betas_glm_loop <- function(data, label, batch, family) {
  # Build design matrix ONCE (biggest single win: avoids G formula parses)
  design <- if (is.null(batch)) {
    stats::model.matrix(~ 0 + label)
  } else {
    stats::model.matrix(~ 0 + label + batch)
  }
  
  # Loop over genes; tryCatch to assign NA on failure (e.g. perfect separation)
  betas <- apply(data, 1, function(y)
    tryCatch(
      stats::glm.fit(x = design, y = y, family = family,
                     intercept = FALSE)$coefficients,
      error = function(e) rep(NA_real_, ncol(design))
    ))
  # betas is K_total x G; we want K x G, so subset to label rows
  betas <- betas[grep("^label", rownames(betas)), , drop = FALSE]
  # rownames(betas) are "labelX" where X is the group; strip the "label" prefix
  rownames(betas) <- sub("^label", "", rownames(betas))
  betas
}

## 1-vs-max(other) log fold-change contrast used by top_markers_glm().
## Input is a K x G matrix; output has the same shape.
betas_to_logfc_1v_max <- function(betas) {
  vapply(
    seq_len(nrow(betas)), function(i)
      betas[i, ] - sparseMatrixStats::colMaxs(
        betas[-i, , drop = FALSE]
      ),
    numeric(ncol(betas))
  ) |> t()
}

utils::globalVariables(c(".dot", "Scores"))
