## Dense/sparse parity contract for smartid.
##
## Every user-facing entry point must accept a sparse matrix and return
## numerically identical results to the dense equivalent. `cal_score()`
## already had a single sparse test covering `idf = "prob", iae = "prob"`,
## which is exactly the one branch that was sparse-aware; the branches
## below were unprotected and are the reason the gap went unnoticed.
##
## Parameterised blocks are expanded with `lapply()` so each case becomes
## its own named `test_that()` and a failure names the offending method
## instead of aborting the rest of the file.

## ---------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------

## Small Poisson counts with a genuinely sparse structure (~30% zeros)
## so `dgCMatrix` is a faithful representation rather than a dense matrix
## in sparse clothing.
sparse_fixture <- function(G = 30L, N = 24L, seed = 42L) {
  set.seed(seed)
  dense <- matrix(rpois(G * N, 1.2), G, N,
                  dimnames = list(paste0("g", seq_len(G)),
                                  paste0("c", seq_len(N))))
  list(
    dense  = dense,
    sparse = as(dense, "CsparseMatrix"),
    label  = rep(c("A", "B", "C"), length.out = N)
  )
}

## Numeric view of any smartid return value, so dense and sparse results
## can be compared with a single `expect_equal()`.
as_num <- function(x) {
  if (is.list(x) && !is.null(x$score)) x <- x$score
  unname(as.matrix(x))
}

## `idf`/`iae` method catalogue: which methods need a group label, and
## which can mathematically preserve sparsity. `m` and `hdb` both build a
## G x N product, so a dense result is the correct expectation there.
IDF_IAE_CASES <- list(
  list(method = "null",     needs_label = FALSE, keeps_sparse = TRUE),
  list(method = "standard", needs_label = FALSE, keeps_sparse = TRUE),
  list(method = "sd",       needs_label = FALSE, keeps_sparse = TRUE),
  list(method = "igm",      needs_label = TRUE,  keeps_sparse = TRUE),
  list(method = "prob",     needs_label = TRUE,  keeps_sparse = TRUE),
  list(method = "rf",       needs_label = TRUE,  keeps_sparse = TRUE),
  list(method = "m",        needs_label = FALSE, keeps_sparse = FALSE),
  list(method = "hdb",      needs_label = FALSE, keeps_sparse = FALSE)
)

## ---------------------------------------------------------------------
## 1. cal_score(): every idf x iae method, dense vs dgCMatrix
## ---------------------------------------------------------------------

invisible(lapply(IDF_IAE_CASES, function(case) {
  test_that(paste0("cal_score matches dense on dgCMatrix for method ",
                   case$method), {
    if (case$method == "hdb") skip_if_not_installed("dbscan")
    fx  <- sparse_fixture()
    par <- if (case$needs_label) list(label = fx$label) else NULL

    r_sparse <- cal_score(fx$sparse, idf = case$method, iae = case$method,
                          par.idf = par, par.iae = par)
    r_dense  <- cal_score(fx$dense,  idf = case$method, iae = case$method,
                          par.idf = par, par.iae = par)

    expect_equal(as_num(r_sparse), as_num(r_dense), tolerance = 1e-10)
    if (case$keeps_sparse) expect_s4_class(r_sparse$score, "CsparseMatrix")
  })
}))

invisible(lapply(
  Filter(function(x) !x$method %in% c("null", "hdb"), IDF_IAE_CASES),
  function(case) {
    test_that(paste0("cal_score intermediates match dense for method ",
                     case$method), {
      fx  <- sparse_fixture()
      par <- if (case$needs_label) list(label = fx$label) else NULL

      r_sparse <- cal_score(fx$sparse, idf = case$method, iae = case$method,
                            par.idf = par, par.iae = par,
                            return.intermediate = TRUE)
      r_dense  <- cal_score(fx$dense,  idf = case$method, iae = case$method,
                            par.idf = par, par.iae = par,
                            return.intermediate = TRUE)

      expect_equal(unname(as.matrix(r_sparse$tf)),
                   unname(as.matrix(r_dense$tf)),  tolerance = 1e-10)
      expect_equal(unname(as.matrix(r_sparse$idf)),
                   unname(as.matrix(r_dense$idf)), tolerance = 1e-10)
      expect_equal(unname(as.matrix(r_sparse$iae)),
                   unname(as.matrix(r_dense$iae)), tolerance = 1e-10)
    })
  }
))

## ---------------------------------------------------------------------
## 2. Internal idf_* / iae_* helpers
## ---------------------------------------------------------------------

## The IDF side is the half that was never converted to Matrix-aware
## accessors; the IAE side already uses `Matrix::rowSums()` /
## `sparseMatrixStats::rowSums2()`.
IDF_IAE_HELPERS <- function(label) list(
  tf       = function(x) smartid:::tf(x, log = TRUE),
  idf      = function(x) smartid:::idf(x),
  idf_m    = function(x) smartid:::idf_m(x),
  idf_sd   = function(x) smartid:::idf_sd(x),
  idf_prob = function(x) smartid:::idf_prob(x, label = label),
  idf_rf   = function(x) smartid:::idf_rf(x, label = label),
  idf_igm  = function(x) smartid:::idf_igm(x, label = label),
  iae      = function(x) smartid:::iae(x),
  iae_m    = function(x) smartid:::iae_m(x),
  iae_sd   = function(x) smartid:::iae_sd(x),
  iae_prob = function(x) smartid:::iae_prob(x, label = label),
  iae_rf   = function(x) smartid:::iae_rf(x, label = label),
  iae_igm  = function(x) smartid:::iae_igm(x, label = label)
)

invisible(lapply(names(IDF_IAE_HELPERS("A")), function(nm) {
  test_that(paste0(nm, "() accepts dgCMatrix and matches dense"), {
    fx <- sparse_fixture()
    f  <- IDF_IAE_HELPERS(fx$label)[[nm]]
    expect_no_error(f(fx$sparse))
    expect_equal(as_num(f(fx$sparse)), as_num(f(fx$dense)), tolerance = 1e-10)
  })
}))

invisible(lapply(c("idf_hdb", "iae_hdb"), function(nm) {
  test_that(paste0(nm, "() accepts dgCMatrix"), {
    skip_if_not_installed("dbscan")
    fx <- sparse_fixture()
    expect_no_error(get(nm, envir = asNamespace("smartid"))(fx$sparse))
  })
}))

## ---------------------------------------------------------------------
## 3. The `features` subset contract
## ---------------------------------------------------------------------

invisible(lapply(c("idf", "idf_sd", "iae", "iae_sd"), function(nm) {
  test_that(paste0(nm, "() honours a features subset on dgCMatrix"), {
    fx <- sparse_fixture()
    fs <- rownames(fx$dense)[1:10]
    f  <- function(x) get(nm, envir = asNamespace("smartid"))(x, features = fs)
    expect_no_error(f(fx$sparse))
    expect_equal(as_num(f(fx$sparse)), as_num(f(fx$dense)), tolerance = 1e-10)
  })
}))

## `idf_sd()` subsets `tfs` before taking row SDs; `iae_sd()` does not,
## so `sd_row` is `nrow(expr)` long while `s_row` is `length(features)`
## long and the closing `log1p()` silently recycles the shorter vector.
##
## Independent oracle: TF needs every gene (column sums run over the full
## matrix), so it is computed once and only then restricted to
## `features` -- exactly what `idf_sd()` does. `sd()` and `rowSums()` are
## used here rather than the package's own helpers so the reference
## shares no code with the implementation under test.
iae_sd_oracle <- function(expr, features) {
  tfs    <- smartid:::tf(expr, log = FALSE)
  sd_row <- apply(tfs[features, , drop = FALSE], 1, sd)
  s_row  <- rowSums(expr[features, , drop = FALSE])
  log1p(sd_row * ncol(expr) / (s_row + 1))
}

## Subset shapes matter: with the leading contiguous rows the recycled
## vector happens to line up, so the first `length(features)` values are
## accidentally correct and only the surplus entries are visible. Any
## other subset also corrupts the values themselves.
FEATURE_SUBSETS <- list(
  list(label = "leading contiguous rows",  idx = 1:10),
  list(label = "trailing contiguous rows", idx = 21:30),
  list(label = "scattered rows",           idx = c(5, 12, 20, 3, 27,
                                                   8, 15, 1, 30, 22))
)

invisible(lapply(FEATURE_SUBSETS, function(case) {
  test_that(paste0("iae_sd returns one value per feature: ", case$label), {
    fx <- sparse_fixture()
    fs <- rownames(fx$dense)[case$idx]
    ## `idf_sd()` is the correct sibling and pins the intended contract.
    expect_length(smartid:::idf_sd(fx$dense, features = fs), length(fs))
    expect_length(smartid:::iae_sd(fx$dense, features = fs), length(fs))
  })

  test_that(paste0("iae_sd matches an independent oracle: ", case$label), {
    fx <- sparse_fixture()
    fs <- rownames(fx$dense)[case$idx]
    expect_equal(unname(smartid:::iae_sd(fx$dense, features = fs)),
                 unname(iae_sd_oracle(fx$dense, fs)),
                 tolerance = 1e-10)
  })
}))

## ---------------------------------------------------------------------
## 4. gs_score family
## ---------------------------------------------------------------------

test_that("gs_score_init accepts dgCMatrix and matches dense", {
  fx <- sparse_fixture()
  fs <- rownames(fx$dense)[1:10]
  expect_no_error(smartid:::gs_score_init(fx$sparse, features = fs))
  expect_equal(
    unname(smartid:::gs_score_init(fx$sparse, features = fs)),
    unname(smartid:::gs_score_init(fx$dense,  features = fs)),
    tolerance = 1e-10
  )
})

test_that("gs_score accepts dgCMatrix for a feature vector", {
  fx <- sparse_fixture()
  fs <- rownames(fx$dense)[1:10]
  expect_no_error(gs_score(fx$sparse, features = fs))
  expect_equal(unname(gs_score(fx$sparse, features = fs)),
               unname(gs_score(fx$dense,  features = fs)),
               tolerance = 1e-10)
})

test_that("gs_score accepts dgCMatrix for a named list of signatures", {
  fx  <- sparse_fixture()
  lst <- list(sig1 = rownames(fx$dense)[1:10],
              sig2 = rownames(fx$dense)[11:20])
  expect_no_error(gs_score(fx$sparse, features = lst))
  expect_equal(as_num(gs_score(fx$sparse, features = lst)),
               as_num(gs_score(fx$dense,  features = lst)),
               tolerance = 1e-10)
})

test_that("gs_score returns a plain numeric vector for dgCMatrix input", {
  fx  <- sparse_fixture()
  out <- gs_score(fx$sparse, features = rownames(fx$dense)[1:10])
  expect_true(is.numeric(out))
  expect_length(out, ncol(fx$sparse))
})

test_that("ova_score_boxplot builds a ggplot from dgCMatrix input", {
  skip_if_not_installed("ggpubr")
  fx <- sparse_fixture()
  expect_no_error(
    ova_score_boxplot(fx$sparse, features = rownames(fx$dense)[1:10],
                      ref.group = "A", label = fx$label)
  )
})

## ---------------------------------------------------------------------
## 5. scale_mgm
## ---------------------------------------------------------------------

invisible(lapply(c(FALSE, TRUE), function(pooled) {
  test_that(paste0("scale_mgm matches dense on dgCMatrix (pooled.sd = ",
                   pooled, ")"), {
    fx <- sparse_fixture()
    expect_no_error(scale_mgm(fx$sparse, label = fx$label, pooled.sd = pooled))
    expect_equal(
      as_num(scale_mgm(fx$sparse, label = fx$label, pooled.sd = pooled)),
      as_num(scale_mgm(fx$dense,  label = fx$label, pooled.sd = pooled)),
      tolerance = 1e-10
    )
  })
}))

test_that("row_scale_zmean matches dense for dgCMatrix input", {
  fx <- sparse_fixture()
  expect_no_error(smartid:::row_scale_zmean(fx$sparse))
  expect_equal(as_num(smartid:::row_scale_zmean(fx$sparse)),
               as_num(smartid:::row_scale_zmean(fx$dense)),
               tolerance = 1e-10)
})

## Row-centring destroys sparsity by construction, so the scaled result
## is necessarily dense. It must still be a class the downstream row
## statistics have methods for: `rowSds`, `rowVars`, `rowMedians`,
## `rowMads` and `rowMaxs` all lack a `dgeMatrix` method, so leaking a
## `dgeMatrix` out of the scaling step breaks `top_markers()`.
invisible(lapply(c("rowMeans2", "rowSums2", "rowSds", "rowVars",
                   "rowMedians", "rowMads", "rowMaxs"), function(op) {
  test_that(paste0("scaled dgCMatrix output supports ", op, "()"), {
    fx     <- sparse_fixture()
    scaled <- smartid:::apply_row_scaling(
      fx$sparse, fx$label, scale = TRUE, use.mgm = TRUE, pooled.sd = FALSE
    )
    fn <- get(op, envir = asNamespace("sparseMatrixStats"))
    expect_no_error(fn(scaled, na.rm = TRUE))
  })
}))

## ---------------------------------------------------------------------
## 6. top_markers
## ---------------------------------------------------------------------

## Order-insensitive comparison on (group, gene) keys, mirroring the
## convention in test-numerical-equivalence.R.
expect_top_markers_equal <- function(got, want, tolerance = 1e-8) {
  key_got  <- paste(got$.dot,  got$Genes,  sep = "|")
  key_want <- paste(want$.dot, want$Genes, sep = "|")
  expect_setequal(key_got, key_want)
  idx <- match(key_want, key_got)
  expect_equal(got$Scores[idx], want$Scores, tolerance = tolerance)
}

TOP_MARKERS_CASES <- list(
  list(label = "glm gaussian",      args = list()),
  list(label = "abs median",        args = list(use.glm = FALSE,
                                                method = "median")),
  list(label = "abs mad",           args = list(use.glm = FALSE,
                                                method = "mad")),
  list(label = "abs mean",          args = list(use.glm = FALSE,
                                                method = "mean")),
  list(label = "no scaling",        args = list(scale = FALSE)),
  list(label = "scale without mgm", args = list(use.mgm = FALSE)),
  list(label = "pooled sd",         args = list(pooled.sd = TRUE))
)

invisible(lapply(TOP_MARKERS_CASES, function(case) {
  test_that(paste0("top_markers matches dense on dgCMatrix: ", case$label), {
    fx   <- sparse_fixture()
    base_args <- list(label = fx$label, n = 5)
    got  <- do.call(top_markers,
                    c(list(data = fx$sparse), base_args, case$args))
    want <- do.call(top_markers,
                    c(list(data = fx$dense),  base_args, case$args))
    expect_top_markers_equal(got, want)
  })
}))

## ---------------------------------------------------------------------
## 6b. Deferred row scaling
## ---------------------------------------------------------------------

## `top_markers_abs()` and `top_markers_glm()` scale the reduced
## G x K / K x G statistic rather than the G x N input, which is what
## keeps a dgCMatrix sparse. The reference below is the materialised
## route it replaced: scale the full matrix first, then reduce. Sparse
## vs dense parity alone cannot catch a mistake that moves both routes
## together, so the identity is pinned directly.

invisible(lapply(c("mean", "median", "mad"), function(method) {
  test_that(paste0("deferred scaling matches materialised: ", method), {
    fx     <- sparse_fixture()
    params <- smartid:::row_scaling_params(fx$sparse, fx$label)
    scaled <- smartid:::apply_row_scaling(fx$sparse, fx$label,
                                          TRUE, TRUE, FALSE)
    expect_equal(
      smartid:::apply_deferred_scaling(
        smartid:::aggregate_rows_by_group(fx$sparse, fx$label, method),
        params, centre = method != "mad"
      ),
      smartid:::aggregate_rows_by_group(scaled, fx$label, method),
      tolerance = 1e-10
    )
  })
}))

test_that("deferred scaling matches materialised: glm 1-vs-max contrast", {
  fx     <- sparse_fixture()
  lab    <- factor(fx$label)
  params <- smartid:::row_scaling_params(fx$sparse, lab)
  scaled <- smartid:::apply_row_scaling(fx$sparse, lab, TRUE, TRUE, FALSE)

  got <- smartid:::apply_deferred_scaling(
    smartid:::betas_to_logfc_1v_max(
      smartid:::fit_label_betas_closed_form(fx$sparse, lab, NULL)
    ),
    params, margin = 2L, centre = FALSE
  )
  want <- smartid:::betas_to_logfc_1v_max(
    smartid:::fit_label_betas_closed_form(scaled, lab, NULL)
  )
  expect_equal(unname(got), unname(want), tolerance = 1e-10)
})

## Covariate designs crossed with the label keep the closed-form fast
## path; one nested inside the label is rank-deficient, so the closed
## form declines and the materialised fallback takes over. Both routes
## must land on the same scores.
BATCH_CASES <- list(
  list(label = "no batch",      batch = function(n) NULL),
  list(label = "two levels",    batch = function(n)
    rep(c("b1", "b2"), length.out = n)),
  list(label = "four levels",   batch = function(n)
    rep(c("b1", "b2", "b3", "b4"), length.out = n)),
  list(label = "batch x donor", batch = function(n) as.character(
    interaction(rep(c("b1", "b2"), length.out = n),
                rep(c("d1", "d2", "d3"), each = 4, length.out = n),
                drop = TRUE))),
  list(label = "confounded with label", batch = function(n)
    paste0("blk_", rep(c("A", "B", "C"), length.out = n)))
)

invisible(lapply(BATCH_CASES, function(case) {
  test_that(paste0("top_markers glm matches dense with batch: ",
                   case$label), {
    fx  <- sparse_fixture(G = 60L, N = 36L)
    bat <- case$batch(ncol(fx$sparse))
    got  <- top_markers(fx$sparse, label = fx$label, n = 5, batch = bat)
    want <- top_markers(fx$dense,  label = fx$label, n = 5, batch = bat)
    expect_top_markers_equal(got, want)
  })
}))

## Guarantees the fallback branch above is genuinely exercised: if this
## design ever became full rank, the confounded case would silently stop
## covering the materialised route.
test_that("a label-confounded design makes the closed form decline", {
  fx  <- sparse_fixture(G = 60L, N = 36L)
  bat <- factor(paste0("blk_", fx$label))
  expect_null(
    smartid:::fit_label_betas_closed_form(fx$sparse, factor(fx$label), bat)
  )
})

## ---------------------------------------------------------------------
## 7. SummarizedExperiment carrying a sparse assay
## ---------------------------------------------------------------------

test_that("cal_score keeps an SE assay sparse and top_markers consumes it", {
  fx <- sparse_fixture()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = fx$sparse),
    colData = data.frame(Group = fx$label)
  )
  expect_no_error(
    se <- cal_score(se, par.idf = list(label = "Group"),
                    par.iae = list(label = "Group"))
  )
  expect_s4_class(SummarizedExperiment::assay(se, "score"), "CsparseMatrix")
  expect_no_error(top_markers(se, label = "Group", n = 5, slot = "score"))
})

## ---------------------------------------------------------------------
## 8. AnyMatrix class union coverage
## ---------------------------------------------------------------------

## `setClassUnion("AnyMatrix", c("matrix", "dgCMatrix"))` rejects every
## other matrix representation Bioconductor hands out, including the
## `SVT_SparseMatrix` that current SingleCellExperiment readers produce.
invisible(lapply(
  list(
    list(name = "dgRMatrix", to = "RsparseMatrix", from = "sparse"),
    list(name = "dgTMatrix", to = "TsparseMatrix", from = "sparse"),
    list(name = "dgeMatrix", to = "unpackedMatrix", from = "dense")
  ),
  function(rep) {
    test_that(paste0("cal_score accepts ", rep$name, " input"), {
      fx <- sparse_fixture()
      x  <- as(fx[[rep$from]], rep$to)
      expect_no_error(cal_score(x, idf = "null", iae = "null"))
    })
  }
))

## `SVT_SparseMatrix` is not a `Matrix` subclass, so widening `AnyMatrix`
## to the `Matrix` virtual class does not reach it. Supporting the format
## current SingleCellExperiment readers emit would mean adding a
## SparseArray dependency, which is a separate decision.
test_that("cal_score accepts an SVT_SparseMatrix", {
  skip("TODO: needs a SparseArray dependency; AnyMatrix covers Matrix only")
  skip_if_not_installed("SparseArray")
  fx  <- sparse_fixture()
  svt <- SparseArray::SVT_SparseArray(fx$dense)
  expect_no_error(cal_score(svt, idf = "null", iae = "null"))
})
