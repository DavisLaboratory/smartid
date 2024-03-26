#' @include plot.R
NULL

#' select markers using mixtools EM method
#'
#' @param top_markers output of [top_markers()]
#' @param column character, specify which column used as group label
#' @param prob numeric, probability cutoff for 1st component classification
#' @param k integer, number of components of mixtures
#' @param ratio numeric, ratio cutoff of 1st component mu to 2nd component mu,
#'     only when ratio > cutoff will return markers for the group
#' @param dist can be one of "norm" and "gamma", specify if to use
#'     [mixtools::normalmixEM()] or [mixtools::gammamixEM()]
#' @param maxit integer, maximum number of iterations for EM
#' @param plot logical, if to plot mixture density and hist
#' @param ... other params for [mixtools::normalmixEM()] or [mixtools::gammamixEM()]
#'
#' @return a list of markers for each group
#' @export
#'
#' @examples
#' set.seed(1000)
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' top_n <- top_markers(data, label = rep(c("A", "B"), 5))
#' markers_mixmdl(top_n, k = 3)
markers_mixmdl <- function(top_markers, column = ".dot",
                           prob = 0.99, k = 3, ratio = 2,
                           dist = c("norm", "gamma"), # gamma takes much longer time, when data is extreme imbalanced using norm
                           maxit = 1e5,
                           plot = FALSE,
                           ...) {
  dist <- match.arg(dist)
  stopifnot(
    "Please provide a valid column name!" = column %in% colnames(top_markers)
  )
  ## split tibble into list
  top_markers <- top_markers |>
    dplyr::group_by(!!ggplot2::sym(column)) |>
    dplyr::group_split() |>
    setNames(unique(top_markers[[column]]) |> sort())

  markers_list <- lapply(top_markers, function(m) {
    if (dist == "norm") {
      mixmdl <- mixtools::normalmixEM(m$Scores, k = k, maxit = maxit, ...)
    } else {
      mixmdl <- mixtools::gammamixEM(m$Scores, k = k, maxit = maxit, ...)
      ## compute mu
      mixmdl$mu <- mixmdl$gamma.pars[1, ] * mixmdl$gamma.pars[2, ] # alpha * beta
    }


    if (plot == TRUE) {
      # # base R plot, easy
      # show(plot(mixmdl, which = 2))
      # abline(v = mixmdl$mu, col = seq_along(mixmdl$mu) + 1)

      ## ggplot
      p <- plot_mm(mixmdl, dist = dist)
      show(p)
    }
    ## get order of mu
    ord <- order(mixmdl$mu, decreasing = TRUE)
    ## return null if ratio of 1st highest and 2nd highest mu is < thres
    if (mixmdl$mu[ord[1]] / mixmdl$mu[ord[2]] < ratio) {
      return(NULL)
    }

    idx <- which(mixmdl$posterior[, ord[1]] > prob)
    m$Genes[idx]
  })

  return(markers_list)
}

#' select markers using mclust EM method
#'
#' @inheritParams markers_mixmdl
#' @param s_thres NULL or numeric, only features with score > threshold will be
#'     returned, if NULL will use 2 * average probability as threshold
#' @param method can be "max.one" or "remove.min", if to only keep features in
#'     1st component or return features not in the last component
#' @param ... other params for [mclust::densityMclust()]
#'
#' @return a list of markers for each group
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' top_n <- top_markers(data, label = rep(c("A", "B"), 5))
#' markers_mclust(top_n)
markers_mclust <- function(top_markers,
                           column = ".dot",
                           prob = 0.99,
                           s_thres = NULL,
                           method = c("max.one", "remove.min"),
                           plot = FALSE,
                           ...) {
  method <- match.arg(method)
  if (is.null(s_thres)) {
    s_thres <- 2 / table(top_markers[[column]])[1]
  }

  stopifnot(
    "Please provide a valid column name!" = column %in% colnames(top_markers)
  )
  ## split tibble into list
  top_markers <- top_markers |>
    dplyr::group_by(!!ggplot2::sym(column)) |>
    dplyr::group_split() |>
    setNames(unique(top_markers[[column]]) |> sort())

  markers_list <- lapply(names(top_markers), function(g) {
    m <- top_markers[[g]]
    s <- m$Scores

    ## cluster cell types based on each gene profile
    cl <- mclust::densityMclust(s, plot = FALSE, ...)

    if (plot == TRUE) {
      ## bootstrap
      boot <- mclust::MclustBootstrap(cl)
      ## set plot par
      par(mfrow = c(2, cl$G))
      ## plot
      plot(boot, what = "mean")
      plot(boot, what = "pro")
      title(g, outer = TRUE, line = -2)
      # reset plot par
      par(mfrow = c(1, 1))

      # show(plot_mm_clust(cl$data, cl$classification))
    }
    # plot(cbind(s, s), col=alpha(cl$classification, 0.4), pch=20)

    if (method == "max.one") {
      flag <- max(cl$classification) ## which cluster has the max expression
      m$Genes[cl$classification == flag & s > s_thres & cl$z[, flag] > prob]
    } else {
      flag <- min(cl$classification) ## which cluster has the min expression
      if (length(unique(cl$classification)) == 1 && min(s, na.rm = TRUE) > s_thres) {
        m$Genes
      } else {
        m$Genes[cl$classification != flag & s > s_thres]
      }
    }
  }) |>
    setNames(names(top_markers))

  return(markers_list)
}

#' select markers using HDBSCAN method
#'
#' @inheritParams markers_mclust
#' @param minPts integer, minimum size of clusters for [dbscan::hdbscan()]
#' @param ... other params for [dbscan::hdbscan()]
#'
#' @return a list of markers for each group
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' top_n <- top_markers(data, label = rep(c("A", "B"), 5))
#' markers_hdbscan(top_n, minPts = 2)
markers_hdbscan <- function(top_markers,
                            column = ".dot",
                            s_thres = NULL,
                            method = c("max.one", "remove.min"),
                            minPts = 5,
                            plot = FALSE,
                            ...) {
  method <- match.arg(method)
  if (is.null(s_thres)) {
    s_thres <- 2 / table(top_markers[[column]])[1]
  }

  stopifnot(
    "Please provide a valid column name!" = column %in% colnames(top_markers)
  )
  ## split tibble into list
  top_markers <- top_markers |>
    dplyr::group_by(!!sym(column)) |>
    dplyr::group_split() |>
    setNames(unique(top_markers[[column]]) |> sort())

  markers_list <- lapply(top_markers, function(m) {
    s <- m$Scores

    ## cluster cell types based on each gene profile
    cl <- dbscan::hdbscan(cbind(s, s), minPts = minPts, ...)

    if (plot == TRUE) {
      ## plot hist and density with mu
      show(plot_mm_clust(s, cl$cluster))
      ## plot hc with stable clusters
      plot(cl, show_flat = TRUE)
    }
    # plot(cbind(s, s), col=alpha(cl$cluster+1, 0.4), pch=20)

    if (method == "max.one") {
      flag <- cl$cluster[which.max(s)] ## which cluster has the max expression
      m$Genes[cl$cluster == flag & s > s_thres]
    } else {
      flag <- cl$cluster[which.min(s)] ## which cluster has the min expression
      if (length(unique(cl$cluster)) == 1 && min(s, na.rm = TRUE) > s_thres) {
        m$Genes
      } else {
        m$Genes[cl$cluster != flag & s > s_thres]
      }
    }
  })

  return(markers_list)
}
