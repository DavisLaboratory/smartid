test_that("marker selection works", {
  ## sim data
  set.seed(1000)
  data <- matrix(rnorm(100), 10, dimnames = list(1:10))
  top_n <- top_markers(data, label = rep(c("A", "B"), 5))

  ## test mixmdl gaussian
  res <- markers_mixmdl(top_n, k = 3)
  expect_type(res, "list")

  ## test mixmdl gamma
  res <- markers_mixmdl(top_n, k = 3, dist = "gamma")
  expect_type(res, "list")

  ## test mclust
  res <- markers_mclust(top_n)
  expect_type(res, "list")

  ## test hdbscan
  res <- markers_hdbscan(top_n, minPts = 2)
  expect_type(res, "list")
})
