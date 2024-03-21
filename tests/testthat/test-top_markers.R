test_that("top markers works", {
  ## sim data
  expr <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
  label <- rep(c("A", "B"), 5)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = expr),
    colData = data.frame(group = label)
  )

  ## test matrix
  ### glm
  res <- top_markers(
    data = expr,
    label = label,
    n = Inf
  )
  expect_s3_class(res, "data.frame")
  ### abs
  res <- top_markers(
    data = expr,
    label = label,
    n = 3,
    use.glm = FALSE
  )
  expect_s3_class(res, "data.frame")
  ### no scale, softmax
  res <- top_markers(
    data = expr,
    label = label,
    n = 3,
    scale = FALSE,
    softmax = FALSE
  )
  expect_s3_class(res, "data.frame")
  ### standard scale
  res <- top_markers(
    data = expr,
    label = label,
    n = 3,
    use.mgm = FALSE
  )
  expect_s3_class(res, "data.frame")

  ## test se
  res <- top_markers(
    data = se,
    label = "group",
    n = 3,
    slot = "counts"
  )
  expect_s3_class(res, "data.frame")
})
