## Independent formula oracles for the IDF/IAE helpers.
##
## The dense code path is the known-good contract: it is what every
## published result was produced with. The sparse conversion replaces
## `base::rowSums()` / `base::colMeans()` with the MatrixGenerics
## generics re-exported by sparseMatrixStats, and an `ifelse()` mask with
## an arithmetic mask, none of which may move a dense number.
##
## Rather than snapshotting current output to an .rds (a golden master
## only restates whatever the code does today, and the package
## deliberately carries no test data), each oracle below re-implements
## the documented formula from scratch in plain base R. It shares no code
## with the implementation, so a disagreement localises the defect to one
## side or the other.
##
## `iae()`, `iae_m()` and `iae_sd()` are already Matrix-aware; including
## them validates the oracles themselves rather than the implementation.

## ---------------------------------------------------------------------
## Oracles
## ---------------------------------------------------------------------

## TF: column-normalise by total counts + 0.01, optional log1p.
## `t(t(x) / cs)` is the transpose-divide-transpose idiom, independent of
## the `sweep()` / `%*% Diagonal()` forms used in the package.
tf_oracle <- function(expr, log = FALSE) {
  cs  <- colSums(expr) + 0.01
  out <- t(t(expr) / cs)
  if (log) log1p(out) else out
}

## Number of cells in which each feature exceeds `thres`.
n_above <- function(expr, thres) {
  vapply(seq_len(nrow(expr)),
         function(i) sum(expr[i, ] > thres), numeric(1))
}

## Per-cell maximum of `stat` over the features expressed in that cell,
## with 0 as the floor -- the explicit form of
## `colMaxs(ifelse(expr > thres, stat, 0))`.
per_cell_max <- function(expr, thres, stat) {
  vapply(seq_len(ncol(expr)),
         function(j) max(c(0, stat[expr[, j] > thres])), numeric(1))
}

## IDF_i = log(1 + n / (n_i + 1))
idf_oracle <- function(expr, features, thres = 0) {
  n_sub <- n_above(expr[features, , drop = FALSE], thres)
  log1p(ncol(expr) / (n_sub + 1))
}

## IDF_{i,j} = log( max_{i' in j}(n_{i'}) / (n_i + 1) )
idf_m_oracle <- function(expr, features, thres = 0) {
  n_sub <- n_above(expr, thres)
  n_max <- per_cell_max(expr, thres, n_sub)
  out   <- outer(1 / (1 + n_sub), n_max)
  dimnames(out) <- dimnames(expr)
  log(out[features, , drop = FALSE])
}

## IDF_i = log(1 + sd(tf_i) * n / (n_i + 1))
idf_sd_oracle <- function(expr, features, log = FALSE, thres = 0) {
  tfs    <- tf_oracle(expr, log = log)
  sd_row <- apply(tfs[features, , drop = FALSE], 1, sd)
  n_sub  <- n_above(expr[features, , drop = FALSE], thres)
  log1p(sd_row * ncol(expr) / (n_sub + 1))
}

## Counts clipped at `thres`, matching pmax0_offset().
offset_oracle <- function(expr, thres) {
  if (thres == 0) expr else pmax(expr - thres, 0)
}

## IAE_i = log(1 + n / (sum_j max(0, N_ij - thres) + 1))
iae_oracle <- function(expr, features, thres = 0) {
  s_row <- rowSums(offset_oracle(expr[features, , drop = FALSE], thres))
  log1p(ncol(expr) / (s_row + 1))
}

## IAE_{i,j} = log(1 + max_{i' in j}(s_{i'}) / (s_i + 1))
iae_m_oracle <- function(expr, features, thres = 0) {
  off   <- offset_oracle(expr, thres)
  s_row <- rowSums(off)
  s_max <- per_cell_max(off, 0, s_row)
  out   <- outer(1 / (1 + s_row), s_max)
  dimnames(out) <- dimnames(expr)
  log1p(out[features, , drop = FALSE])
}

## IAE_i = log(1 + sd(tf_i) * n / (s_i + 1)); note the SD is taken over
## the *requested* features, exactly as idf_sd() does.
iae_sd_oracle <- function(expr, features, log = FALSE, thres = 0) {
  tfs    <- tf_oracle(expr, log = log)
  sd_row <- apply(tfs[features, , drop = FALSE], 1, sd)
  s_row  <- rowSums(offset_oracle(expr[features, , drop = FALSE], thres))
  log1p(sd_row * ncol(expr) / (s_row + 1))
}

## ---------------------------------------------------------------------
## Fixture
## ---------------------------------------------------------------------

## Non-negative counts so `pmax0_offset()`'s `thres == 0` short-circuit
## and an explicit `pmax(x, 0)` cannot diverge.
oracle_fixture <- function(G = 24L, N = 18L, seed = 7L) {
  set.seed(seed)
  dense <- matrix(rpois(G * N, 1.5), G, N,
                  dimnames = list(paste0("g", seq_len(G)),
                                  paste0("c", seq_len(N))))
  list(dense = dense, sparse = as(dense, "CsparseMatrix"))
}

## ---------------------------------------------------------------------
## Unlabelled IDF / IAE helpers vs their oracles
## ---------------------------------------------------------------------

ORACLE_CASES <- list(
  list(name = "idf",    impl = "idf",    oracle = idf_oracle),
  list(name = "idf_m",  impl = "idf_m",  oracle = idf_m_oracle),
  list(name = "idf_sd", impl = "idf_sd", oracle = idf_sd_oracle),
  list(name = "iae",    impl = "iae",    oracle = iae_oracle),
  list(name = "iae_m",  impl = "iae_m",  oracle = iae_m_oracle),
  list(name = "iae_sd", impl = "iae_sd", oracle = iae_sd_oracle)
)

## `thres = 0` is the default; `thres = 1` exercises the offset branch.
invisible(lapply(ORACLE_CASES, function(case) {
  lapply(c(0, 1), function(thres) {
    test_that(paste0(case$name, "() on dense input matches its formula",
                     " (thres = ", thres, ")"), {
      fx <- oracle_fixture()
      fs <- rownames(fx$dense)
      impl <- get(case$impl, envir = asNamespace("smartid"))
      expect_equal(
        unname(as.matrix(impl(fx$dense, features = fs, thres = thres))),
        unname(as.matrix(case$oracle(fx$dense, features = fs, thres = thres))),
        tolerance = 1e-10
      )
    })
  })
}))

## The same equivalence must hold once the helpers accept sparse input,
## which is what pins the conversion to a value-preserving one.
invisible(lapply(ORACLE_CASES, function(case) {
  test_that(paste0(case$name, "() on dgCMatrix matches its formula"), {
    fx <- oracle_fixture()
    fs <- rownames(fx$dense)
    impl <- get(case$impl, envir = asNamespace("smartid"))
    expect_equal(
      unname(as.matrix(impl(fx$sparse, features = fs, thres = 0))),
      unname(as.matrix(case$oracle(fx$dense, features = fs, thres = 0))),
      tolerance = 1e-10
    )
  })
}))

## ---------------------------------------------------------------------
## The two rewrites that must be value-neutral, stated directly
## ---------------------------------------------------------------------

## The sparse conversion swaps `base::rowSums()` for
## `sparseMatrixStats::rowSums2()` and `base::colMeans()` for
## `colMeans2()`. Those are the MatrixGenerics generics re-exported by
## sparseMatrixStats, and for an ordinary matrix they resolve to the
## matrixStats backend, so dense values must be bit-identical.
test_that("rowSums2/colMeans2 agree with base on an ordinary matrix", {
  fx   <- oracle_fixture()
  mask <- fx$dense > 0
  expect_identical(sparseMatrixStats::rowSums2(mask), base::rowSums(mask))
  expect_identical(sparseMatrixStats::colMeans2(fx$dense, na.rm = TRUE),
                   base::colMeans(fx$dense, na.rm = TRUE))
})

## `idf_m()`'s `ifelse(expr > thres, n_sub, 0)` becomes `(expr > thres) *
## n_sub`, mirroring the form `iae_m()` already uses. Both recycle
## `n_sub` down the columns; the arithmetic form additionally preserves
## sparsity. Equivalence needs `n_sub >= 0`, which holds because it is a
## count.
test_that("arithmetic mask equals ifelse mask for a non-negative statistic", {
  fx    <- oracle_fixture()
  n_sub <- base::rowSums(fx$dense > 0)
  expect_identical(
    (fx$dense > 0) * n_sub,
    ifelse(fx$dense > 0, n_sub, 0)
  )
})

## ---------------------------------------------------------------------
## cal_score dense output must not move
## ---------------------------------------------------------------------

## Composition check: for the unlabelled methods the score is simply
## tf * idf * iae, so an independent recomposition catches any drift in
## how `combine_tf_idf_iae()` broadcasts the factors.
invisible(lapply(c("standard", "m", "sd"), function(method) {
  test_that(paste0("cal_score dense output recomposes for idf/iae = ",
                   method), {
    fx <- oracle_fixture()
    fs <- rownames(fx$dense)

    got <- cal_score(fx$dense, tf = "logtf", idf = method, iae = method)$score

    idf_fn <- switch(method, standard = idf_oracle, m = idf_m_oracle,
                     sd = idf_sd_oracle)
    iae_fn <- switch(method, standard = iae_oracle, m = iae_m_oracle,
                     sd = iae_sd_oracle)
    want <- tf_oracle(fx$dense, log = TRUE) *
      idf_fn(fx$dense, features = fs, thres = 0) *
      iae_fn(fx$dense, features = fs, thres = 0)

    expect_equal(unname(as.matrix(got)), unname(as.matrix(want)),
                 tolerance = 1e-10)
  })
}))
