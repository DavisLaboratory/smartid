test_that("gs_score works", {
  ## sim data
  set.seed(123)
  data <- matrix(rnorm(100), 10, dimnames = list(seq_len(10)))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = data)
  )

  ## test matrix
  res <- gs_score(data, features = seq_len(3))
  expect_type(res, "double")

  ## test matrix, feature list
  res <- gs_score(data, features = list(a = seq_len(3), b = seq_len(5)))
  expect_type(res, "list")

  ## test se
  res <- gs_score(se, features = seq_len(3), slot = "counts")
  expect_type(res$score, "double")

  ## test se, feature list
  res <- gs_score(se,
    features = list(a = seq_len(3), b = seq_len(5)),
    slot = "counts"
  )
  expect_type(colnames(SummarizedExperiment::colData(res)), "character")
})
