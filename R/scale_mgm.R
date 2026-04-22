#' scale by mean of group mean for imbalanced data
#'
#' @details
#' \deqn{z=\frac{x-\frac{\sum_k^{n_D}(\mu_k)}{n_D}}{s}}
#' where \eqn{\mu_k} is the mean of x in \eqn{k^{th}} class, and \eqn{n_D} is
#' the number of classes, \eqn{s} is the standard deviation of x,
#' when `pooled.sd` is set to be TRUE, \eqn{s} will be replaced with
#' \eqn{s_{pooled}}, \eqn{s_{pooled}=\sqrt{\frac{\sum_k^{n_D}{(n_k-1){s_k}^2}}{\sum_k^{n_D}{n_k}-k}}}
#'
#' @param expr matrix
#' @param label a vector of group label
#' @param pooled.sd logical, if to use pooled SD for scaling
#'
#' @return scaled matrix
#' @export
#'
#' @examples
#' scale_mgm(matrix(rnorm(100), 10), label = rep(letters[1:2], 5))
scale_mgm <- function(expr, label, pooled.sd = FALSE) {
  ## Cache column indices per group once; the group-mean and pooled-SD
  ## paths previously recomputed `label == i` inside every `vapply`
  ## iteration, which is O(K * N) in scan cost.
  idx_by_grp <- split(seq_len(ncol(expr)), label)

  sds <- if (isTRUE(pooled.sd)) {
    row_pool_sds_from_idx(expr, idx_by_grp)
  } else {
    sparseMatrixStats::rowSds(expr, na.rm = TRUE)
  }

  # Compute per-group means, then average them to get the mean of group means (MGM).
  group_means <- vapply(idx_by_grp, function(cols)
    sparseMatrixStats::rowMeans2(expr[, cols, drop = FALSE], na.rm = TRUE),
    numeric(nrow(expr)))
  mgm <- rowMeans(group_means, na.rm = TRUE) # mean of per-group means

  # scale
  ## Single broadcast + single allocation: `(expr - mgm) * inv_sd`
  ## collapses the prior two temporaries ((expr - mgm), then divide) into
  ## one. Division-by-zero is guarded by the additive epsilon.
  inv_sd <- 1 / (sds + 1e-8)
  (expr - mgm) * inv_sd
}


## Row-wise pooled SDs given pre-computed per-group column indices.
row_pool_sds_from_idx <- function(expr, idx_by_grp) {
  # Compute per-group variances, then combine them with the group sizes to get the pooled SD.
  vars <- vapply(idx_by_grp, function(cols)
    sparseMatrixStats::rowVars(expr[, cols, drop = FALSE], na.rm = TRUE),
    numeric(nrow(expr)))
  # Group sizes: number of columns in each group.
  ng <- lengths(idx_by_grp)
  as.numeric(sqrt((vars %*% cbind(ng - 1)) / sum(ng - 1)))
}

## Back-compat shim: preserves the old `row_pool_sds(expr, label)`
## internal call signature in case any caller still uses it.
row_pool_sds <- function(expr, label) {
  row_pool_sds_from_idx(expr, split(seq_len(ncol(expr)), label))
}
