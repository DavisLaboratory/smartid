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
                             scale = TRUE,
                             use.mgm = TRUE,
                             softmax = TRUE,
                             ...) {
  if (use.glm == TRUE) {
    data <- top_markers_glm(
      data = data,
      label = label,
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
#' @param data matrix, features in row and samples in column
#' @param label vector, group labels
#' @param n integer, number of returned top genes for each group
#' @param method character, specify metric to compute, can be one of "median",
#'     "mad", "mean"
#' @param scale logical, if to scale data by row
#' @param use.mgm logical, if to scale data using [scale_mgm()]
#' @param softmax logical, if to apply softmax transformation on output
#'
#' @return a tibble with feature names, group labels and ordered processed scores
#' @export
#'
#' @examples
#' data <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
#' top_markers_abs(data, label = rep(c("A", "B"), 5))
top_markers_abs <- function(data, label, n = 10,
                            method = c("median", "mad", "mean"),
                            scale = TRUE, use.mgm = TRUE,
                            softmax = TRUE) {
  method <- match.arg(method)
  if (scale & use.mgm) {
    data <- scale_mgm(expr = data, label = label)
  } else if (scale & !use.mgm) {
    ## scale scores on rows
    # mu_s <- sparseMatrixStats::rowMeans2(data, na.rm = TRUE)
    # sd_s <- sparseMatrixStats::rowSds(data, na.rm = TRUE)
    # data <- (data - mu_s) / sd_s

    data <- t(scale(t(data)))
    data[is.na(data)] <- 0 # assign 0 to NA when sd = 0
  }

  data <- data |>
    t() |>
    as.data.frame() |>
    dplyr::group_by(.dot = label) |> ## group by label
    dplyr::summarise_all(method, na.rm = TRUE) |> ## aggregate scores
    tidyr::gather("Genes", "Scores", -`.dot`) |> ## transform into long data
    dplyr::group_by(`.dot`) ## group by label again

  if (softmax == TRUE) {
    data <- data |>
      # dplyr::mutate(Scores = Scores / sd(Scores, na.rm = TRUE)) |> # norm by sd
      # dplyr::mutate(Scores = sigmoid(Scores)) |> # sigmoid
      # dplyr::mutate(Scores = tanh(Scores)) |> # tanh
      dplyr::mutate(Scores = softmax(Scores)) # softmax
  }

  data <- dplyr::slice_max(data, Scores, n = n) ## extract top n markers

  # ## softmax
  # if(softmax == TRUE)
  #   data <- dplyr::mutate(data, Scores = softmax(Scores))

  return(data)
}

#' calculate group mean score using glm and order genes based on scores difference
#'
#' @param data matrix, features in row and samples in column
#' @param label vector, group labels
#' @param n integer, number of returned top genes for each group
#' @param family family for glm, details in [stats::glm()]
#' @param scale logical, if to scale data by row
#' @param use.mgm logical, if to scale data using [scale_mgm()]
#' @param softmax logical, if to apply softmax transformation on output
#'
#' @return a tibble with feature names, group labels and ordered processed scores
#' @export
#'
#' @examples
#' data <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
#' top_markers_glm(data, label = rep(c("A", "B"), 5))
top_markers_glm <- function(data, label, n = 10,
                            family = gaussian(), # score are continuous non-negative, can use gamma or inverse.gaussian, if continuous and unbounded use gaussian, if discrete use poisson, if binary or proportions between [0,1] or binary freq counts use binomial
                            scale = TRUE, use.mgm = TRUE,
                            # log = TRUE,
                            softmax = TRUE) {
  label <- factor(label) # factorize label

  ## scale
  if (scale & !use.mgm) {
    ## scale scores on rows
    # mu_s <- sparseMatrixStats::rowMeans2(data, na.rm = TRUE)
    # sd_s <- sparseMatrixStats::rowSds(data, na.rm = TRUE)
    # data <- (data - mu_s) / sd_s

    data <- t(scale(t(data)))
    data[is.na(data)] <- 0 # assign 0 to NA when sd = 0
  } else if (scale & use.mgm) {
    data <- scale_mgm(expr = data, label = label)
  }

  # ## log score
  # if(log == TRUE) {
  #   data <- log(data + 1e-8)
  # }

  ## estimate betas based on given group
  betas <- apply(data, 1, \(s) glm(s ~ 0 + label, family = family)$coef)
  rownames(betas) <- gsub("label", "", rownames(betas))

  # ## compute logFC (1 vs all mean) for each group
  # betas <- apply(betas, 2, \(x) x - (sum(x) - x)/(length(x) - 1))

  ## compute logFC (1 vs max excluding self) for each group
  betas <- vapply(
    seq_len(nrow(betas)), \(i)
    betas[i, ] - sparseMatrixStats::colMaxs(betas[-i, , drop = FALSE]),
    rep(1, ncol(betas))
  ) |>
    t()
  rownames(betas) <- levels(label)

  data <- data.frame(.dot = rownames(betas), betas) |>
    tidyr::pivot_longer(-`.dot`, names_to = "Genes", values_to = "Scores") |>
    dplyr::group_by(`.dot`) ## group by label again

  if (softmax == TRUE) {
    data <- data |>
      # dplyr::mutate(Scores = Scores / sd(Scores, na.rm = TRUE)) |> # norm by sd
      # dplyr::mutate(Scores = sigmoid(Scores)) |> # sigmoid
      # dplyr::mutate(Scores = tanh(Scores)) |> # tanh
      dplyr::mutate(Scores = softmax(Scores)) # softmax
  }

  data <- dplyr::slice_max(data, Scores, n = n) ## extract top n markers

  # ## softmax
  # if(softmax == TRUE)
  #   data <- dplyr::mutate(data, Scores = softmax(Scores))

  return(data)
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

utils::globalVariables(c(".dot", "Scores"))
