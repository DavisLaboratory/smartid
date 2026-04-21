#' @include tf_idf_iae_wrappers.R
NULL

#' Calculate scores of each cell on given features
#'
#' @param score matrix, features in row and samples in column
#' @param features vector, feature names to compute score
#'
#' @return a vector of score
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' gs_score_init(data, 1:5)
gs_score_init <- function(score, features = NULL) {
  if (is.null(features)) features <- rownames(score)

  ## check features
  if (!all(features %in% rownames(score))) {
    warning(sprintf(
      "Feature %s is not in score!\n",
      setdiff(features, rownames(score))
    ))
  }
  features <- intersect(rownames(score), features)
  stopifnot("less than 2 features are in score rows!" = length(features) > 1)

  ## calculate mean score of features
  m_score <- colMeans(score[features, , drop = FALSE], na.rm = TRUE)
  return(m_score)
}

#' @title Get names of available IDF and IAE methods
#'
#' @description Returns a named vector of IDF/IAE methods
#' @return names of methods implemented
#' @export
#'
#' @examples
#' idf_iae_methods()
idf_iae_methods <- function() {
  return(sort(c(
    "label probability" = "prob", "label relative frequency" = "rf",
    "label IGM" = "igm", "null" = "null",
    "unlabel max" = "m", "unlabel SD" = "sd",
    "unlabel HDBSCAN" = "hdb", "unlabel standard" = "standard"
  )))
}

#' Calculate score for each feature in each cell
#'
#' @param expr a count matrix, features in row and cells in column
#' @param tf a character, specify the TF method to use, can be "tf" or "logtf"
#' @param idf a character, specify the IDF method to use. Available methods can
#'     be accessed using [idf_iae_methods()]
#' @param iae a character, specify the IAE method to use. Available methods can
#'     be accessed using [idf_iae_methods()]
#' @param par.idf other parameters for specified IDF methods
#' @param par.iae other parameters for specified IAE methods
#' @param return.intermediate logical, if TRUE the returned list also contains
#'     the intermediate `tf`, `idf` and `iae` objects. Default `FALSE` keeps
#'     only the combined `score` to avoid the memory overhead of three extra
#'     feature-by-cell matrices on large inputs.
#'
#' @return a list always containing `score`; when `return.intermediate = TRUE`
#'     the list additionally contains `tf`, `idf` and `iae`.
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' label <- sample(c("A", "B"), 10, replace = TRUE)
#' smartid:::cal_score_init(data,
#'   par.idf = list(label = label),
#'   par.iae = list(label = label)
#' )
cal_score_init <- function(expr, tf = c("logtf", "tf"),
                           idf = "prob", iae = "prob",
                           par.idf = NULL, par.iae = NULL,
                           return.intermediate = FALSE) {
  ## check
  tf <- match.arg(tf)
  idf <- match.arg(idf, choices = idf_iae_methods())
  iae <- match.arg(iae, choices = idf_iae_methods())
  stopifnot(
    "par.idf must be a named list or NULL" = is.null(par.idf) | is.list(par.idf),
    "par.iae must be a named list or NULL" = is.null(par.iae) | is.list(par.iae),
    "return.intermediate must be a single logical" =
      is.logical(return.intermediate) && length(return.intermediate) == 1L
  )

  ## compute tf
  tf <- tf(expr, log = (tf == "logtf"))

  ## compute idf
  if (idf == "null") {
    idf <- 1
  } else {
    idf <- ifelse(idf == "standard", "idf", paste0("idf_", idf))
    idf <- do.call(idf, c(list(expr = expr), par.idf))
  }

  ## compute iae
  if (iae == "null") {
    iae <- 1
  } else {
    iae <- ifelse(iae == "standard", "iae", paste0("iae_", iae))
    iae <- do.call(iae, c(list(expr = expr), par.iae))
  }

  ## combined score via per-group column-block broadcast; avoids
  ## materialising the full G x N copy that a naive `tf * idf * iae`
  ## would trigger when either factor is a G x K compact matrix.
  score <- combine_tf_idf_iae(tf, idf, iae,
                              label_idf = par.idf$label,
                              label_iae = par.iae$label)

  if (isTRUE(return.intermediate)) {
    return(list(score = score, tf = tf, idf = idf, iae = iae))
  }
  list(score = score)
}

## Per-group column-block composition of score = tf * idf * iae.
##
## Each factor can be one of:
##   * scalar 1 (the "null" path),
##   * a G-vector (cell-independent, e.g. idf/iae/idf_sd/iae_sd/idf_igm/iae_igm),
##   * a G x K compact matrix with colnames = unique labels (idf_prob, idf_rf,
##     iae_prob, iae_rf after Phase B),
##   * a full G x N matrix (idf_m, iae_m, idf_hdb, iae_hdb — the latter two
##     expand internally to preserve their legacy contract).
##
## When at least one factor is compact we loop over groups and broadcast
## only the active slice into the corresponding columns of `score`. When
## both factors are cell-independent or full-cell we fall back to the
## direct algebraic form.
combine_tf_idf_iae <- function(tf_mat, idf_obj, iae_obj,
                               label_idf = NULL, label_iae = NULL) {
  N <- ncol(tf_mat)
  is_compact <- function(x) {
    if (is.null(dim(x))) return(FALSE)
    !is.null(colnames(x)) && ncol(x) < N
  }
  idf_is_gk <- is_compact(idf_obj)
  iae_is_gk <- is_compact(iae_obj)

  if (!idf_is_gk && !iae_is_gk) {
    ## no compact factor; direct algebra preserves sparsity for scalar /
    ## G-vector factors and only densifies when a factor is already G x N.
    return(tf_mat * idf_obj * iae_obj)
  }

  label <- label_idf %||% label_iae
  stopifnot(
    "par.idf$label or par.iae$label is required when idf/iae return compact G x K matrices" =
      !is.null(label),
    "length(label) must equal ncol(expr)" = length(label) == N
  )
  label_ch <- as.character(label)

  ## Slice a factor for a given (column subset, group name) pair.
  slice_factor <- function(x, cols, group_name) {
    if (is.null(dim(x))) return(x)             # scalar or G-vector
    if (ncol(x) == N) return(x[, cols, drop = FALSE])  # G x N full
    x[, group_name]                             # G x K compact
  }

  score <- tf_mat
  for (g in unique(label_ch)) {
    cols <- which(label_ch == g)
    if (!length(cols)) next
    idf_g <- slice_factor(idf_obj, cols, g)
    iae_g <- slice_factor(iae_obj, cols, g)
    score[, cols] <- score[, cols, drop = FALSE] * idf_g * iae_g
  }
  score
}

## Null-coalescing helper (kept local to avoid a new Imports).
`%||%` <- function(x, y) if (is.null(x)) y else x
