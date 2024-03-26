test_that("wrappers work", {
  ## sim data
  set.seed(123)
  expr <- matrix(rpois(100, 2), 10, dimnames = list(1:10, letters[1:10]))
  label <- sample(c("A", "B"), 10, replace = TRUE)

  # test TF
  #--------------
  res <- tf(expr, log = TRUE)
  expect_type(res, "double")

  # test IDF
  #--------------
  ## test idf
  res <- idf(expr)
  expect_type(res, "double")

  ## test idf_m
  res <- idf_m(expr)
  expect_type(res, "double")

  ## test idf_sd
  res <- idf_sd(expr)
  expect_type(res, "double")

  ## test idf_hdb
  res <- idf_hdb(expr)
  expect_type(res, "double")

  ## test idf_rf
  res <- idf_rf(expr, label = label)
  expect_type(res, "double")

  ## test idf_prob
  res <- idf_prob(expr, label = label)
  expect_type(res, "double")

  ## test idf_igm
  res <- idf_igm(expr, label = label)
  expect_type(res, "double")

  # test IAE
  #--------------
  ## test iae
  res <- iae(expr)
  expect_type(res, "double")

  ## test iae_m
  res <- iae_m(expr)
  expect_type(res, "double")

  ## test iae_sd
  res <- iae_sd(expr)
  expect_type(res, "double")

  ## test iae_hdb
  res <- iae_hdb(expr)
  expect_type(res, "double")

  ## test iae_rf
  res <- iae_rf(expr, label = label)
  expect_type(res, "double")

  ## test iae_prob
  res <- iae_prob(expr, label = label)
  expect_type(res, "double")

  ## test iae_igm
  res <- iae_igm(expr, label = label)
  expect_type(res, "double")
})
