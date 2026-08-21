#################################################
#-----------------Helpers----------------------#
#################################################

## For each row i and each column k of `mat`, return the maximum of
## `mat[i, -k]`. Vectorised O(G * K) replacement for the nested
## `apply(mat, 1, function(x) max(x[names(x) != type]))` pattern used in
## the multi-class branches of the labelled IDF/IAE helpers.
rowwise_notin_max <- function(mat) {
  stopifnot(!is.null(dim(mat)))
  G <- nrow(mat)
  K <- ncol(mat)
  if (K == 1L) {
    ## no "other" column exists; return -Inf to mark degenerate input
    out <- matrix(-Inf, nrow = G, ncol = 1L, dimnames = dimnames(mat))
    return(out)
  }
  row_max  <- sparseMatrixStats::rowMaxs(mat, na.rm = TRUE)
  argmax   <- max.col(mat, ties.method = "first")
  masked   <- mat
  masked[cbind(seq_len(G), argmax)] <- -Inf
  row_max2 <- sparseMatrixStats::rowMaxs(masked, na.rm = TRUE)
  ## broadcast row_max across all K columns, then overwrite with row_max2
  ## wherever column index == argmax (i.e. the excluded-self case).
  notin <- matrix(row_max, nrow = G, ncol = K)
  for (k in seq_len(K)) {
    hit <- argmax == k
    if (any(hit)) notin[hit, k] <- row_max2[hit]
  }
  dimnames(notin) <- dimnames(mat)
  notin
}

## Sparse-preserving equivalent of `pmax(x - thres, 0)`. For dgCMatrix
## we mutate the non-zero slot in place and drop structural zeros,
## avoiding the dense allocation that `x[x < 0] <- 0` would trigger.
## `thres == 0` short-circuits because scRNA-seq counts are already
## non-negative, which is the common default path.
pmax0_offset <- function(x, thres = 0) {
  if (thres == 0) return(x)
  if (methods::is(x, "sparseMatrix")) {
    x@x <- pmax(x@x - thres, 0)
    return(Matrix::drop0(x))
  }
  out <- x - thres
  out[out < 0] <- 0
  out
}

#################################################
#-----------------TF variants------------------#
#################################################

### $$\mathbf{TF_{i,j}}=\frac{N_{i,j}}{\sum_j{N_{i,j}}}$$

## term frequency
#' compute term/feature frequency within each cell
#'
#' @details
#' \deqn{\mathbf{TF_{i,j}}=\frac{N_{i,j}}{\sum_j{N_{i,j}}}}
#' where \eqn{N_{i,j}} is the counts of feature i in cell j.
#'
#' @param expr a count matrix, features in row and cells in column
#' @param log logical, if to do log-transformation
#'
#' @return a matrix of term/gene frequency. For a `dgCMatrix` input the
#'     returned object preserves sparsity (`log1p(0) == 0`); dense input
#'     returns a dense matrix.
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::tf(data)
tf <- function(expr, log = FALSE) {
  ## Column sums: use Matrix-aware accessor so dgCMatrix stays sparse.
  cs <- Matrix::colSums(expr, na.rm = TRUE) + 0.01
  if (methods::is(expr, "sparseMatrix")) {
    ## Right-multiplying by a diagonal with 1/cs scales each column
    ## without materialising a dense copy. log1p preserves sparsity
    ## because log1p(0) == 0.
    out <- expr %*% Matrix::Diagonal(x = 1 / cs)
    dimnames(out) <- dimnames(expr)
  } else {
    out <- sweep(expr, 2, cs, FUN = "/")
  }
  if (log) out <- log1p(out)
  out
}

#################################################
#-----------------IDF variants------------------#
#################################################

## -----------------unlabeled-------------------##

## inverse document frequency
### $$\mathbf{IDF_i} = log(1+\frac{n}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' standard inverse cell frequency
#'
#' @details
#' \deqn{\mathbf{IDF_i} = log(1+\frac{n}{n_i+1})}
#' where \eqn{n} is the total number of cells, \eqn{n_i} is the number of cells
#' containing feature i.
#'
#' @inheritParams idf_rf
#'
#' @return a vector of inverse cell frequency score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf(data)
idf <- function(expr, features = NULL, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr) ## number of total obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  ## `expr > thres` is an lgCMatrix for sparse input, which
  ## `base::rowSums()` rejects; the re-exported generic handles both
  ## backends. Mirrors the accessor already used by `iae_sd()`.
  n_sub <- sparseMatrixStats::rowSums2(
    expr[features, , drop = FALSE] > thres
  ) ## num of obs contain feature i > thres

  idf <- log1p(n_obs / (n_sub + 1))
  return(idf)
}

## inverse document frequency max
### $$\mathbf{IDF_{i,j}} = log(\frac{max_{\{i^{'}\in j\}}(n_{i^{'}})}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse cell frequency: max
#'
#' @details
#' \deqn{\mathbf{IDF_{i,j}} = log(\frac{max_{\{i^{'}\in j\}}(n_{i^{'}})}{n_i+1})}
#' where \eqn{i} is the feature \eqn{i} and \eqn{i^{'}} is the feature except
#' \eqn{i}, \eqn{n_i} is the number of cells containing feature i, and
#' \eqn{n_{i^{'}}} is the number of cells containing feature \eqn{i^{'}}.
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse cell frequency score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_m(data)
idf_m <- function(expr, features = NULL, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr, probs = 0.25, na.rm = TRUE)
  above <- expr > thres
  n_sub <- sparseMatrixStats::rowSums2(above) ## num of obs contain feature i > thres

  ## For each cell: max over the features expressed in that cell of their
  ## document frequency. The arithmetic mask replaces `ifelse()`, which
  ## silently coerces an lgCMatrix to a plain logical vector and so drops
  ## `dim`. Multiplying instead keeps both the dimensions and the sparsity.
  ## `n_sub` is a count, hence non-negative, so the two forms agree
  ## elementwise. Same form as `iae_m()`.
  n_max <- sparseMatrixStats::colMaxs(above * n_sub)

  idf <- matrix(1 / (1 + n_sub), ncol = 1) %*% matrix(n_max, nrow = 1)
  dimnames(idf) <- dimnames(expr)

  idf <- log(idf[features, , drop = FALSE])
  return(idf)
}

## inverse document frequency sd
### $$\mathbf{IDF_i} = log(1+sd(tf_{i})*\frac{n}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse cell frequency using standard deviation (SD)
#'
#' @details
#' \deqn{\mathbf{IDF_i} = log(1+sd(tf_{i})*\frac{n}{n_i+1})}
#' where \eqn{tf_i} is the term frequency of feature \eqn{i}, see details in
#' [tf()], \eqn{n} is the total number of cells and \eqn{n_i} is the number of
#' cells containing feature \eqn{i}.
#'
#' @inheritParams idf_rf
#' @inheritParams tf
#'
#' @return a vector of inverse cell frequency score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_sd(data)
idf_sd <- function(expr, features = NULL, log = FALSE, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr) ## number of total obs

  tfs <- tf(expr, log = log)

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  ## Matrix-aware accessor: see the note in `idf()`.
  n_sub <- sparseMatrixStats::rowSums2(
    expr[features, , drop = FALSE] > thres
  ) ## num of obs contain feature i > thres
  sd_row <- sparseMatrixStats::rowSds(tfs[features, , drop = FALSE], na.rm = TRUE)

  idf <- log1p(sd_row * n_obs / (n_sub + 1))
  return(idf)
}

## inverse document frequency using hdbscan cluster as label
#' inverse document frequency using hdbscan cluster as label
#'
#' @details
#' Details as [idf_prob()].
#'
#' @inheritParams idf_rf
#' @param minPts integer, minimum size of clusters, default 2.
#'     Details in [dbscan::hdbscan()].
#' @param ... parameters for [dbscan::hdbscan()]
#'
#' @return a matrix of inverse cell frequency score
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_hdb(data)
idf_hdb <- function(expr, features = NULL, multi = TRUE,
                    thres = 0, minPts = 2, ...) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  ## initially compute naive tf-idf (sparse-preserving via tf() helper)
  tfidf <- tf(expr)[features, , drop = FALSE] *
    idf(expr, features = features, thres = thres)

  ## cluster obs based on given features
  cluster <- dbscan::hdbscan(Matrix::t(tfidf), minPts = minPts, ...)$cluster
  ## factor cluster
  cluster <- factor(cluster)

  ## `idf_prob()` now returns a compact G x K matrix; `idf_hdb()` owns the
  ## cluster labels and is not reachable from `cal_score_init()`, so we
  ## expand here to keep the G x N external contract.
  idf_compact <- idf_prob(
    expr = expr, features = features,
    label = cluster, multi = multi,
    thres = thres
  )
  idf_compact[, as.character(cluster), drop = FALSE]
}

## ------------------labeled--------------------##

## labeled inverse document frequency: relative frequency
### $$\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}})$$

#' labeled inverse cell frequency: relative frequency
#'
#' @details
#' \deqn{\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}})}
#' where \eqn{n_{i,j\in D}} is the number of cells containing feature \eqn{i} in
#' class \eqn{D}, \eqn{n_{j\in D}} is the total number of cells in class \eqn{D},
#' \eqn{\hat D} is the class except \eqn{D}.
#'
#' @param expr a matrix, features in row and cells in column
#' @param features vector, feature names or indexes to compute
#' @param label vector, group label of each cell
#' @param multi logical, if to compute based on binary (FALSE) or multi-class (TRUE)
#' @param thres numeric, cell only counts when expr > threshold, default 0
#'
#' @return a matrix of inverse cell frequency score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_rf(data, label = sample(c("A", "B"), 10, replace = TRUE))
idf_rf <- function(expr, features = NULL, label,
                   multi = TRUE, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  df_n <- expr[features, , drop = FALSE] > thres ## if contain feature i > thres

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- vapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(df_n[, label == type, drop = FALSE],
      na.rm = TRUE
    )
  }, rep(1, nrow(df_n))) ## mean counts for each gene in the group
  if (multi == TRUE) {
    ## G x K: row-wise max over "other" columns; O(G * K) vectorisation
    mean_row_notin <- rowwise_notin_max(mean_row_in)
  } else {
    mean_row_notin <- vapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(df_n[, label != type, drop = FALSE],
        na.rm = TRUE
      )
    }, rep(1, nrow(df_n))) ## doc freq for each gene not in group for bi-class
  }

  ## Return compact G x K matrix; `cal_score_init()` handles per-group
  ## broadcast so we never materialise a G x N copy via `[, label]`.
  idf <- log1p(mean_row_in / (mean_row_notin + 1e-8))
  colnames(idf) <- colnames(mean_row_in)
  idf
}

## labeled inverse document frequency: probability based
### $$\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}}\frac{n_{i,j\in D}}{n_{j\in D}})$$
### modified from
### $$\mathbf{IDF_{i,j}} = log(1+\frac{A}{max(B)+1}*\frac{A}{C+1})$$
### A denotes the number of cells belonging to category D where the gene i occurs at least once; B denotes the number of cells not belonging to category D where the gene i occurs at least once; C denotes the number of cells belonging to category D where the gene i does not occur; D denotes the number of cells not belonging to category D where the gene i does not occur.

#' labeled inverse cell frequency: probability based
#'
#' @details
#' \deqn{\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}}\frac{n_{i,j\in D}}{n_{j\in D}})}
#' where \eqn{n_{i,j\in D}} is the number of cells containing feature \eqn{i} in
#' class \eqn{D}, \eqn{n_{j\in D}} is the total number of cells in class \eqn{D},
#' \eqn{\hat D} is the class except \eqn{D}.
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse cell frequency score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_prob(data, label = sample(c("A", "B"), 10, replace = TRUE))
idf_prob <- function(expr, features = NULL, label,
                     multi = TRUE, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  df_n <- expr[features, , drop = FALSE] > thres ## if contain feature i > thres

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- vapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(df_n[, label == type, drop = FALSE],
      na.rm = TRUE
    )
  }, rep(1, nrow(df_n))) ## mean counts for each gene in the group
  if (multi == TRUE) {
    ## G x K: row-wise max over "other" columns; O(G * K) vectorisation
    mean_row_notin <- rowwise_notin_max(mean_row_in)
  } else {
    mean_row_notin <- vapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(df_n[, label != type, drop = FALSE],
        na.rm = TRUE
      )
    }, rep(1, nrow(df_n))) ## doc freq for each gene not in group for bi-class
  }

  ## Return compact G x K; `cal_score_init()` broadcasts per group.
  idf <- log1p(mean_row_in^2 / (mean_row_notin + 1e-8))
  colnames(idf) <- colnames(mean_row_in)
  idf

  return(idf)
}

## labeled inverse document frequency IGM
### $$\mathbf{IGM_i} = log(1+\lambda\frac{max(n_{i,j\in D})_{k}}{\sum_{k}^{K}((n_{i,j\in D})_{k}*r_{k})+e^{-8}})$$
### $$\mathbf{k}: type\ in\ total\ group\ K$$
### $$\mathbf{r_{k}}: rank\ of\ n_{i,j\in D})_{k}\ in\ total\ group\ K$$

#' labeled inverse cell frequency: IGM
#'
#' @details
#' \deqn{\mathbf{IGM_i} = log(1+\lambda\frac{max(n_{i,j\in D})_{k}}{\sum_{k}^{K}((n_{i,j\in D})_{k}*r_{k})+e^{-8}})}
#' where \eqn{\lambda} is the hyper parameter, \eqn{n_{i,j\in D}} is the number
#' of cells containing feature \eqn{i} in class \eqn{D}, \eqn{r_k} is the rank
#' of \eqn{n_{i,j\in D}}.
#'
#' @inheritParams idf_rf
#' @param lambda numeric, hyperparameter for IGM
#'
#' @return a vector of inverse gravity moment score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_igm(data, label = sample(c("A", "B"), 10, replace = TRUE))
idf_igm <- function(expr, features = NULL, label, lambda = 7, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  df_n <- expr[features, , drop = FALSE] > thres ## if contain feature i > thres

  mean_row_in <- vapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(df_n[, label == type, drop = FALSE],
      na.rm = TRUE
    )
  }, rep(1, nrow(df_n))) ## mean counts for each gene in the group
  ## calculate igm for each gene
  igm <- apply(mean_row_in, 1, function(x) {
    max(x) / (sum(x * rank(-x), na.rm = TRUE) + 1e-8)
  })

  idf <- log1p(lambda * igm) ## TF-IDF scores
  return(idf)
}

#################################################
#-----------------IAE variants------------------#
#################################################

## -----------------unlabeled-------------------##

## inverse average expression
### $$\mathbf{IAE_i} = log(1+\frac{n}{\sum_j^n\hat N_{i,j}+1})$$
### $$n: total\ counts\ of\ documents;\ \hat N_{i,j}: max(0, N_{i,j} - threshold)$$

#' standard inverse average expression
#'
#' @details
#' \deqn{\mathbf{IAE_i} = log(1+\frac{n}{\hat N_{i,j}+1})}
#' where \eqn{n} is the total number of cells, \eqn{N_{i,j}} is the counts of
#' feature \eqn{i} in cell \eqn{j}.
#'
#' @inheritParams idf_rf
#'
#' @return a vector of inverse average expression score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae(data)
iae <- function(expr, features = NULL, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr) ## number of total obs

  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  expr_offset <- pmax0_offset(expr[features, , drop = FALSE], thres)
  s_row <- Matrix::rowSums(expr_offset, na.rm = TRUE) ## per-feature total counts

  iae <- log1p(n_obs / (s_row + 1))
  return(iae)
}

## inverse average expression max
### $$\mathbf{IAE_{i,d}} = log(1+\frac{max_{\{i^{'}\in d\}}(n_{i^{'}})}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse average expression: max
#'
#' @details
#' \deqn{\mathbf{IAE_{i,j}} = log(1+\frac{max_{\{i^{'}\in j\}}(n_{i^{'}})}{\sum_{j = 1}^{n} max(0, N_{i,j} - threshold)+1})}
#' where \eqn{i} is the feature \eqn{i} and \eqn{i^{'}} is the feature except
#' \eqn{i}, \eqn{N_{i,j}} is the counts of feature \eqn{i} in cell \eqn{j}, and
#' \eqn{n_{i^{'}}} is \eqn{\sum_{j = 1}^{n} sign(N_{i,j} > threshold)}.
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse average expression score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_m(data)
iae_m <- function(expr, features = NULL, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  expr_offset <- pmax0_offset(expr, thres)
  s_row <- Matrix::rowSums(expr_offset, na.rm = TRUE) ## per-feature total counts
  ## For each cell: max over features of s_row restricted to features
  ## that are expressed > thres in that cell. `ifelse()` on the logical
  ## mask preserves sparsity for dgCMatrix.
  nonzero <- expr_offset > 0
  s_max <- sparseMatrixStats::colMaxs(nonzero * s_row, na.rm = TRUE)

  iae <- matrix(1 / (1 + s_row), ncol = 1) %*% matrix(s_max, nrow = 1)
  dimnames(iae) <- dimnames(expr)
  iae <- log1p(iae[features, , drop = FALSE])
  return(iae)
}

## inverse average expression sd
### $$\mathbf{IAE} = log(1+sd(tf_{i})*\frac{n}{\sum_{j=1}^{n}N_{i,j}+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse average expression using standard deviation (SD)
#'
#' @details
#' \deqn{\mathbf{IAE} = log(1+sd(tf_{i})*\frac{n}{\sum_{j=1}^{n}max(0,N_{i,j})+1})}
#' where \eqn{tf_i} is the term frequency of feature \eqn{i}, see details in
#' [tf()], \eqn{n} is the total number of cells and \eqn{N_{i,j}} is the counts
#' of feature \eqn{i} in cell \eqn{j}.
#'
#' @inheritParams idf_rf
#' @inheritParams tf
#'
#' @return a vector of inverse average expression score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_sd(data)
iae_sd <- function(expr, features = NULL, log = FALSE, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr) ## number of obs

  tfs <- tf(expr, log = log)

  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  expr_offset <- pmax0_offset(expr[features, , drop = FALSE], thres)
  s_row <- sparseMatrixStats::rowSums2(expr_offset, na.rm = TRUE) ## summed counts for each gene
  sd_row <- sparseMatrixStats::rowSds(
    tfs[features, , drop = FALSE],
    na.rm = TRUE
  )

  iae <- log1p(sd_row * n_obs / (s_row + 1)) ## IDF scores
  return(iae)
}

## inverse average expression using hdbscan cluster as label
#' inverse average expression using hdbscan cluster as label
#'
#' @details
#' Details as [iae_prob()].
#'
#' @inheritParams idf_rf
#' @param minPts integer, minimum size of clusters, default 2.
#'     Details in [dbscan::hdbscan()].
#' @param ... parameters for [dbscan::hdbscan()]
#'
#' @return a matrix of inverse average expression score
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_hdb(data)
iae_hdb <- function(expr, features = NULL, multi = TRUE,
                    thres = 0, minPts = 2, ...) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  ## initially compute naive tf-iae (sparse-preserving via tf() helper)
  tfidf <- tf(expr)[features, , drop = FALSE] *
    iae(expr, features = features, thres = thres)

  ## cluster obs based on given features
  cluster <- dbscan::hdbscan(Matrix::t(tfidf), minPts = minPts, ...)$cluster
  ## factor cluster
  cluster <- factor(cluster)

  ## Same rationale as `idf_hdb()`: expand G x K compact to G x N here
  ## since cluster labels are not visible to `cal_score_init()`.
  iae_compact <- iae_prob(
    expr = expr, features = features,
    label = cluster, multi = multi,
    thres = thres
  )
  iae_compact[, as.character(cluster), drop = FALSE]
}

## ------------------labeled--------------------##

## labeled inverse average expression: relative frequency
### $$\mathbf{IAE} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}})$$
###   or
### $$\mathbf{IAE} = log(1+\frac{mean(N_{i,j\in D})}{mean(N_{i,j\notin D})+ e^{-8}})$$

#' labeled inverse average expression: relative frequency
#'
#' @details
#' \deqn{\mathbf{IAE} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}})}
#' where \eqn{N_{i,j\in D}} is the counts of feature \eqn{i} in cell \eqn{j}
#' within class \eqn{D}, and \eqn{\hat D} is the class except \eqn{D}.
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse average expression score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_rf(data, label = sample(c("A", "B"), 10, replace = TRUE))
iae_rf <- function(expr, features = NULL, label,
                   multi = TRUE, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  expr_offset <- pmax0_offset(expr[features, , drop = FALSE], thres)

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- vapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(expr_offset[, label == type, drop = FALSE],
      na.rm = TRUE
    )
  }, rep(1, nrow(expr_offset))) ## mean counts for each gene in the group
  if (multi == TRUE) {
    mean_row_notin <- rowwise_notin_max(mean_row_in)
  } else {
    mean_row_notin <- vapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(expr_offset[, label != type, drop = FALSE],
        na.rm = TRUE
      )
    }, rep(1, nrow(expr_offset))) ## mean counts for each gene not in group for bi-class
  }

  ## Return compact G x K; `cal_score_init()` broadcasts per group.
  iae <- log1p(mean_row_in / (mean_row_notin + 1e-8))
  colnames(iae) <- colnames(mean_row_in)
  iae
}

## labeled inverse average expression: probability based
### $$\mathbf{IAE_{i,j}} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}}*mean(N_{i,j\in D}))$$
### modified from
### $$\mathbf{IDF_{i,j}} = log(1+\frac{A}{max(B)+1}*\frac{A}{C+1})$$
### A denotes the number of cells belonging to category D where the gene i occurs at least once; B denotes the number of cells not belonging to category D where the gene i occurs at least once; C denotes the number of cells belonging to category D where the gene i does not occur; D denotes the number of cells not belonging to category D where the gene i does not occur.

#' labeled inverse average expression: probability based
#'
#' @details
#' \deqn{\mathbf{IAE_{i,j}} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}}*mean(N_{i,j\in D}))}
#' where \eqn{N_{i,j\in D}} is the counts of feature \eqn{i} in cell \eqn{j}
#' within class \eqn{D}, and \eqn{\hat D} is the class except \eqn{D}.
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse average expression score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_prob(data, label = sample(c("A", "B"), 10, replace = TRUE))
iae_prob <- function(expr, features = NULL, label,
                     multi = TRUE, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  expr_offset <- pmax0_offset(expr[features, , drop = FALSE], thres)

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- vapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(expr_offset[, label == type, drop = FALSE],
      na.rm = TRUE
    )
  }, rep(1, nrow(expr_offset))) ## mean counts for each gene in the group
  if (multi == TRUE) {
    mean_row_notin <- rowwise_notin_max(mean_row_in)
  } else {
    mean_row_notin <- vapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(expr_offset[, label != type, drop = FALSE],
        na.rm = TRUE
      )
    }, rep(1, nrow(expr_offset))) ## mean counts for each gene not in group for bi-class
  }

  ## Return compact G x K; `cal_score_init()` broadcasts per group.
  iae <- log1p(mean_row_in^2 / (mean_row_notin + 1e-8))
  colnames(iae) <- colnames(mean_row_in)
  iae
}

## labeled inverse average expression IGM
### $$\mathbf{IGM_i} = log(1+\lambda\frac{max(mean(N_{i,j\in D})_{k})}{\sum_{k}^{K}(mean(N_{i,j\in D})_{k}*r_{k})+e^{-8}})$$
### $$\mathbf{k}: type\ in\ total\ group\ K$$
### $$\mathbf{r_{k}}: rank\ of\ mean(N_{i,j\in D})_{k}\ in\ total\ group\ K$$

#' labeled inverse average expression: IGM
#'
#' @details
#' \deqn{\mathbf{IGM_i} = log(1+\lambda\frac{max(mean(N_{i,j\in D})_{k})}{\sum_{k}^{K}(mean(N_{i,j\in D})_{k}*r_{k})+e^{-8}})}
#' where \eqn{\lambda} is the hyper parameter, \eqn{N_{i,j\in D}} is the counts
#' of feature \eqn{i} in cell \eqn{j} within class \eqn{D}, and \eqn{r_k} is the
#' rank of \eqn{mean(N_{i,j\in D})}.
#'
#' @inheritParams idf_rf
#' @param lambda numeric, hyperparameter for IGM
#'
#' @return a vector of inverse gravity moment score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_igm(data, label = sample(c("A", "B"), 10, replace = TRUE))
iae_igm <- function(expr, features = NULL, label, lambda = 7, thres = 0) {
  if (is.null(features)) features <- seq_len(nrow(expr))

  # thres <- sparseMatrixStats::rowQuantiles(expr[features, , drop = FALSE], probs = 0.25, na.rm = TRUE)
  expr_offset <- pmax0_offset(expr[features, , drop = FALSE], thres)

  mean_row_in <- vapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(expr_offset[, label == type, drop = FALSE],
      na.rm = TRUE
    )
  }, rep(1, nrow(expr_offset))) ## mean counts for each gene in the group
  ## calculate igm for each gene
  igm <- apply(mean_row_in, 1, function(x) {
    max(x) / (sum(x * rank(-x), na.rm = TRUE) + 1e-8)
  })

  iae <- log1p(lambda * igm) ## TF-IDF scores
  return(iae)
}
