test_that("score works", {
  ## sim data
  set.seed(123)
  expr <- matrix(rpois(100, 2), 10, dimnames = list(1:10, letters[1:10]))
  label <- sample(c("A", "B"), 10, replace = TRUE)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = expr),
    colData = data.frame(group = label)
  )

  ## test matrix
  ### label
  res <- cal_score(expr,
    par.idf = list(label = label),
    par.iae = list(label = label)
  )
  expect_type(res, "list")
  res <- cal_score(expr,
    tf = "tf", par.idf = list(label = label),
    par.iae = list(label = label)
  )
  expect_type(res, "list")
  ### unlabel
  res <- cal_score(expr, idf = "sd", iae = "sd")
  expect_type(res, "list")

  ## test sce
  ### label
  res <- cal_score(se,
    par.idf = list(label = "group"),
    par.iae = list(label = "group")
  )
  expect_s4_class(res, "SummarizedExperiment")
  ### unlabel
  res <- cal_score(se, idf = "sd", iae = "sd")
  expect_s4_class(res, "SummarizedExperiment")
})
