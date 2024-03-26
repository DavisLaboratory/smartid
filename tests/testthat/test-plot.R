test_that("plot works", {
  ## sim data
  data <- matrix(rnorm(100), 10, dimnames = list(1:10))
  label <- rep(c("A", "B"), 5)

  ## individual feature score plot
  p <- sin_score_boxplot(
    data = data,
    features = seq_len(2),
    ref.group = "A",
    label = label
  )
  expect_true(is.ggplot(p))

  ## boxplot of overall score for given features
  p <- ova_score_boxplot(
    data = data,
    features = seq_len(5),
    ref.group = "A",
    label = label,
    method = "t.test"
  )
  expect_true(is.ggplot(p))

  ## boxplot of overall score for given features
  top_n <- top_markers(data, label = label)
  p <- score_barplot(
    top_markers = top_n,
    f_list = list(A = paste0("X", seq_len(3)))
  )
  expect_true(is.ggplot(p))

  ## test plot_mm for gaussian
  set.seed(123)
  mixmdl <- mixtools::normalmixEM(rnorm(50), k = 2)
  p <- plot_mm(mixmdl = mixmdl, dist = "norm")
  expect_true(is.ggplot(p))

  ## test plot_mm for gamma
  set.seed(123)
  mixmdl <- mixtools::gammamixEM(rgamma(50, 1), k = 2)
  p <- plot_mm(mixmdl = mixmdl, dist = "gamma")
  expect_true(is.ggplot(p))

  ## test plot_mm_clust for gamma
  set.seed(123)
  p <- plot_mm_clust(score = rnorm(50), clust = rep(c("A", "B"), 25))
  expect_true(is.ggplot(p))
})
