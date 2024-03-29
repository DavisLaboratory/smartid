#' scale by mean of group mean for imbalanced data
#'
#' @details
#' \deqn{z=\frac{x-\frac{\sum_k^{n_D}(\mu_k)}{n_D}}{s_{pooled}}}
#' where \eqn{s_{pooled}=\sqrt{\frac{\sum_k^{n_D}{(n_k-1){s_k}^2}}{\sum_k^{n_D}{n_k}-k}}}
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
  if (pooled.sd) {
    ## compute pooled sds
    sds <- vapply(
      unique(label), \(i)
      sparseMatrixStats::rowVars(expr[, label == i, drop = FALSE],
        na.rm = TRUE
      ),
      rep(1, nrow(expr))
    ) # get vars of each group
    ng <- table(label)[unique(label)] # get group sizes in the same order
    sds <- sds %*% cbind(ng - 1)
    sds <- as.numeric(sqrt(sds / sum(ng - 1)))
  } else {
    ## compute overall sds
    sds <- sparseMatrixStats::rowSds(expr, na.rm = TRUE)

    # ## compute group sds
    # sds <- vapply(unique(label), \(i)
    #               sparseMatrixStats::rowSds(expr[, label == i, drop = FALSE],
    #                                         na.rm = TRUE),
    #               rep(1, nrow(expr))
    #        ) # get sds of each group
    # sds <- sparseMatrixStats::rowMeans2(sds)
  }

  ## compute group means
  mgm <- vapply(
    unique(label), \(i)
    sparseMatrixStats::rowMeans2(expr[, label == i, drop = FALSE],
      na.rm = TRUE
    ),
    rep(1, nrow(expr))
  ) |> # get mean of each group
    rowMeans(na.rm = TRUE) # get mean of group mean

  ## scale
  expr <- (expr - mgm) / (sds + 1e-8)

  # expr[is.na(expr)] <- 0 # assign 0 to NA when sd = 0

  return(expr)
}
