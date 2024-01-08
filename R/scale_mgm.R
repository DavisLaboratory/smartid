#' scale by mean of group mean in case extreme unbalanced data
#'
#' @param expr matrix
#' @param label a vector of group label
#'
#' @return scaled matrix
#' @export
#'
#' @examples
#' scale_mgm(matrix(rnorm(100), 10), label = rep(letters[1:2], 5))
scale_mgm <- function(expr, label) {
  ## compute sds
  sds <- sparseMatrixStats::rowSds(expr, na.rm = TRUE)
  # sds <- sapply(unique(label), \(i)
  #               sparseMatrixStats::rowSds(expr[, label == i], na.rm = TRUE)
  #        ) # get mean of each group
  # colnames(sds) <- unique(label)

  ## compute group means
  mgm <- sapply(unique(label), \(i)
                sparseMatrixStats::rowMeans2(expr[, label == i], na.rm = TRUE)) |>  # get mean of each group
    rowMeans(na.rm = TRUE) # get mean of group mean

  ## scale
  expr <- (expr - mgm) / (sds + 1e-8)

  # expr[is.na(expr)] <- 0 # assign 0 to NA when sd = 0

  return(expr)
}
