#--------------plot for marker selection---------------#

## plot density and hist for EM mixture model in mixtools
plot_mm <- function(mixmdl, dist = c("norm", "gamma")) {
  dist <- match.arg(dist)

  p <- ggplot(data.frame(x = mixmdl$x), aes(x = x)) +
    geom_histogram()
  if (dist == "norm") {
    p <- p +
      mapply(
        function(mean, sd, lambda, col) {
          stat_function(
            fun = function(x) {
              (dnorm(x, mean = mean, sd = sd)) * lambda
            },
            col = col
          )
        },
        mean = mixmdl$mu,
        sd = mixmdl$sigma,
        lambda = mixmdl$lambda,
        col = seq_along(mixmdl$mu) + 1
      ) +
      geom_vline(xintercept = mixmdl$mu, col = seq_along(mixmdl$mu) + 1) +
      labs(x = "score", y = "density") +
      theme_bw()
  } else {
    ## compute mu
    mixmdl$mu <- mixmdl$gamma.pars[1, ] * mixmdl$gamma.pars[2, ] # alpha * beta

    p <- p +
      mapply(
        function(mean, alpha, beta, lambda, col) {
          stat_function(
            fun = function(x) {
              (dgamma(x, shape = alpha, rate = beta)) * lambda
            },
            col = col
          )
        },
        mean = mixmdl$mu,
        alpha = mixmdl$gamma.pars["alpha", ],
        beta = mixmdl$gamma.pars["beta", ],
        lambda = mixmdl$lambda,
        col = seq_along(mixmdl$mu) + 1
      ) +
      geom_vline(xintercept = mixmdl$mu, col = seq_along(mixmdl$mu) + 1) +
      labs(x = "score", y = "density") +
      theme_bw()
  }

  p
}

## plot density and hist for EM mixture model in mclust or clusters in hdbscan
plot_mm_clust <- function(score, clust) {
  clust <- factor(clust)
  p <- data.frame(Scores = score, Comp = clust) |>
    ggplot(aes(x = Scores, group = Comp, col = Comp)) +
    geom_histogram(aes(y = ..density.., fill = Comp), alpha = 0.5) +
    geom_density() +
    geom_vline(
      data = stack(tapply(score, clust, mean)),
      aes(xintercept = values, col = ind)
    ) +
    labs(col = "Component", fill = "Component") +
    theme_classic()

  return(p)
}


#--------------plot for features---------------#
#' boxplot of split single feature score
#'
#' @param data matrix, features in row and samples in column
#' @param features vector, feature names to plot
#' @param ref.group character, reference group name
#' @param label vector, group labels
#' @param method character, statistical test to use,
#'               details in [ggpubr::stat_compare_means()]
#'
#' @return faceted ggplot object
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' sin_score_boxplot(data, 1:2, ref.group = "A", label = rep(c("A", "B"), 5))
sin_score_boxplot <- function(data, features = NULL,
                              ref.group, label,
                              method = "t.test") {
  if (is.null(features)) features <- rownames(data)
  data[features, ] |>
    as.matrix() |>
    as.data.frame() |>
    dplyr::add_rownames("Gene") |>
    setNames(c("Gene", as.character(label))) |>
    tidyr::pivot_longer(-"Gene", names_to = "Type", values_to = "Score") |>
    dplyr::mutate(Type = gsub("\\.\\d+$", "", Type)) |>
    ggplot(aes(x = Type, y = Score, fill = Type)) +
    # geom_boxplot() +
    geom_violin() +
    stat_summary(fun = median, geom = "crossbar") +
    facet_wrap(~Gene, scales = "free") +
    ggpubr::stat_compare_means(
      label = "p.signif", method = method,
      ref.group = ref.group, label.y.npc = 1
    ) +
    theme_classic() +
    theme(axis.text.x = element_blank())
}

#' boxplot of features overall score
#'
#' @inheritParams sin_score_boxplot
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' ova_score_boxplot(data, 1:5, ref.group = "A", label = rep(c("A", "B"), 5))
ova_score_boxplot <- function(data, features,
                              ref.group, label,
                              method = "t.test") {
  data.frame(
    Score = gs_score_init(data, features = features),
    Group = label
  ) |>
    ggplot(aes(x = Group, y = Score, col = Group)) +
    geom_boxplot() +
    # geom_jitter(aes(col = Group), alpha = 0.2) +
    ggpubr::stat_compare_means(
      label = "p.signif", method = method,
      ref.group = ref.group, label.y.npc = 1
    ) +
    theme_classic() +
    theme(axis.text.x = element_blank())
}

#' barplot of processed score
#'
#' @param top_markers output of [top_markers()]
#' @param column character, specify which column used as group label
#' @param f_list a named list of markers
#' @param n numeric, number of returned top genes for each group
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' top_n <- top_markers(data, label = rep(c("A", "B"), 5))
#' score_barplot(top_n)
score_barplot <- function(top_markers, column = ".dot", f_list, n = 30) {
  ## set f_list as features label for color
  if (missing(f_list)) {
    f_list <- list(Features = top_markers$Genes)
  }

  # ## get top n markers
  # top_markers <- top_markers(
  #   data[features, ],
  #   label = label,
  #   n = n,
  #   ...
  # )

  ## extract top n markers
  top_markers <- dplyr::slice_max(top_markers, Scores, n = n)

  ## add markers type
  top_markers <- merge(top_markers, stack(f_list),
    by.x = "Genes",
    by.y = "values", all.x = TRUE
  )

  ## plot
  ggplot(
    top_markers,
    aes(
      y = tidytext::reorder_within(Genes, Scores, !!ggplot2::sym(column)),
      x = Scores,
      fill = ind
    )
  ) +
    geom_bar(stat = "identity") +
    facet_wrap(ggplot2::sym(column), scales = "free") +
    labs(y = "Genes", fill = "Markers Type") +
    tidytext::scale_y_reordered() +
    theme_classic()
}

utils::globalVariables(c(
  "Group", "Score", "x", "Scores", "Comp", "Genes",
  "..density..", "stack", "values", "ind", "Type"
))
