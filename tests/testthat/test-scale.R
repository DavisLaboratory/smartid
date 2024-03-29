test_that("scale works", {
  ## test scale using overall sd
  res <- scale_mgm(matrix(rnorm(100), 10), label = rep(letters[1:2], 5))
  expect_type(res, "double")

  ## test scale using pooled sd
  res <- scale_mgm(matrix(rnorm(100), 10),
    label = rep(letters[1:2], 5),
    pooled.sd = TRUE
  )
  expect_type(res, "double")
})
