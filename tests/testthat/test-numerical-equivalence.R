## Numerical regression tests for the memory-optimization refactor.
##
## These tests compare current output against pre-refactor snapshots
## produced by `tests/testthat/testdata/capture_legacy_snapshots.R` on the
## same HEAD that first introduced this file. Any refactor that claims to
## preserve scoring behaviour must keep all `all.equal()` calls green.
##
## Tolerance is deliberately tight (1e-10) because the refactor only
## rearranges algebra and subsetting order; no stochastic steps are
## involved upstream of the matrices compared here.

snap_path <- test_path("testdata", "legacy_scores.rds")

skip_if_no_snapshot <- function() {
  if (!file.exists(snap_path)) {
    skip(paste0("Legacy snapshot missing: ", snap_path))
  }
}

test_that("cal_score on dense matrix reproduces legacy score", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  res <- cal_score(
    data    = snap$inputs$counts_dense,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = as.character(snap$inputs$label)),
    par.iae = list(label = as.character(snap$inputs$label)),
    return.intermediate = TRUE
  )
  expect_equal(
    unname(as.matrix(res$score)),
    unname(as.matrix(snap$cal_score$matrix_out$score)),
    tolerance = 1e-10
  )
})

test_that("cal_score default no longer stores intermediates in metadata", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(
      Group = snap$inputs$label,
      Batch = snap$inputs$batch
    )
  )
  out <- cal_score(
    data    = se,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group")
  )
  expect_null(out@metadata$tf)
  expect_null(out@metadata$idf)
  expect_null(out@metadata$iae)
  expect_equal(
    unname(as.matrix(SummarizedExperiment::assay(out, "score"))),
    unname(snap$cal_score$se_assay),
    tolerance = 1e-10
  )
})

test_that("cal_score return.intermediate=TRUE restores legacy metadata", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(
      Group = snap$inputs$label,
      Batch = snap$inputs$batch
    )
  )
  out <- cal_score(
    data    = se,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group"),
    return.intermediate = TRUE
  )
  md <- out@metadata
  expect_false(is.null(md$tf))
  expect_false(is.null(md$idf))
  expect_false(is.null(md$iae))
  ## tf stays G x N; direct comparison.
  expect_equal(
    unname(as.matrix(md$tf)),
    unname(snap$cal_score$se_tf),
    tolerance = 1e-10
  )
  ## Phase B: idf/iae for labelled prob/rf methods now return compact
  ## G x K matrices. Expand via the Group label to recover the legacy
  ## G x N representation and compare element-wise.
  label_ch <- as.character(snap$inputs$label)
  expect_equal(
    unname(as.matrix(md$idf)[, label_ch, drop = FALSE]),
    unname(snap$cal_score$se_idf),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.matrix(md$iae)[, label_ch, drop = FALSE]),
    unname(snap$cal_score$se_iae),
    tolerance = 1e-10
  )
})

test_that("top_markers GLM path reproduces legacy scores", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(
      Group = snap$inputs$label,
      Batch = snap$inputs$batch
    )
  )
  se <- cal_score(
    data    = se,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group")
  )
  tm <- top_markers(
    data    = se,
    label   = "Group",
    n       = 10L,
    use.glm = TRUE,
    slot    = "score"
  )
  tm_ref <- snap$top_markers$glm
  ## order-sensitive comparison on the same keys
  key_new <- paste(tm$.dot,       tm$Genes,       sep = "|")
  key_ref <- paste(tm_ref$.dot,   tm_ref$Genes,   sep = "|")
  expect_setequal(key_new, key_ref)
  ## align and compare scores
  idx <- match(key_ref, key_new)
  expect_equal(tm$Scores[idx], tm_ref$Scores, tolerance = 1e-8)
})

test_that("top_markers GLM with batch reproduces legacy scores", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(
      Group = snap$inputs$label,
      Batch = snap$inputs$batch
    )
  )
  se <- cal_score(
    data    = se,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group")
  )
  tm <- top_markers(
    data    = se,
    label   = "Group",
    n       = 10L,
    use.glm = TRUE,
    batch   = "Batch",
    slot    = "score"
  )
  tm_ref <- snap$top_markers$glm_batch
  key_new <- paste(tm$.dot, tm$Genes, sep = "|")
  key_ref <- paste(tm_ref$.dot, tm_ref$Genes, sep = "|")
  expect_setequal(key_new, key_ref)
  idx <- match(key_ref, key_new)
  expect_equal(tm$Scores[idx], tm_ref$Scores, tolerance = 1e-8)
})

test_that("top_markers abs mean path reproduces legacy scores", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(
      Group = snap$inputs$label,
      Batch = snap$inputs$batch
    )
  )
  se <- cal_score(
    data    = se,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group")
  )
  tm <- top_markers(
    data    = se,
    label   = "Group",
    n       = 10L,
    use.glm = FALSE,
    method  = "mean",
    slot    = "score"
  )
  tm_ref <- snap$top_markers$abs_mean
  key_new <- paste(tm$.dot, tm$Genes, sep = "|")
  key_ref <- paste(tm_ref$.dot, tm_ref$Genes, sep = "|")
  expect_setequal(key_new, key_ref)
  idx <- match(key_ref, key_new)
  expect_equal(tm$Scores[idx], tm_ref$Scores, tolerance = 1e-8)
})

test_that("top_markers abs median path reproduces legacy scores", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(
      Group = snap$inputs$label,
      Batch = snap$inputs$batch
    )
  )
  se <- cal_score(
    data    = se,
    tf      = "logtf",
    idf     = "prob",
    iae     = "prob",
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group")
  )
  tm <- top_markers(
    data    = se,
    label   = "Group",
    n       = 10L,
    use.glm = FALSE,
    method  = "median",
    slot    = "score"
  )
  tm_ref <- snap$top_markers$abs_median
  key_new <- paste(tm$.dot, tm$Genes, sep = "|")
  key_ref <- paste(tm_ref$.dot, tm_ref$Genes, sep = "|")
  expect_setequal(key_new, key_ref)
  idx <- match(key_ref, key_new)
  expect_equal(tm$Scores[idx], tm_ref$Scores, tolerance = 1e-8)
})

## ---------------------------------------------------------------------
## Edge cases specific to the memory-optimization refactor
## ---------------------------------------------------------------------

test_that("dgCMatrix input yields identical scores to dense input", {
  skip_if_no_snapshot()
  snap   <- readRDS(snap_path)
  dense  <- snap$inputs$counts_dense
  sparse <- as(dense, "CsparseMatrix")
  lab    <- as.character(snap$inputs$label)

  r_dense <- cal_score(
    dense,  tf = "logtf", idf = "prob", iae = "prob",
    par.idf = list(label = lab), par.iae = list(label = lab),
    return.intermediate = TRUE
  )
  r_sparse <- cal_score(
    sparse, tf = "logtf", idf = "prob", iae = "prob",
    par.idf = list(label = lab), par.iae = list(label = lab),
    return.intermediate = TRUE
  )
  expect_s4_class(r_sparse$score, "dgCMatrix")
  expect_true(is.matrix(r_dense$score))
  expect_equal(
    as.matrix(r_sparse$score), as.matrix(r_dense$score),
    tolerance = 1e-10
  )
  expect_equal(
    as.matrix(r_sparse$tf), as.matrix(r_dense$tf),
    tolerance = 1e-10
  )
})

test_that("top_markers gaussian closed-form matches glm apply-loop", {
  skip_if_no_snapshot()
  snap <- readRDS(snap_path)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = snap$inputs$counts_dense),
    colData = data.frame(Group = snap$inputs$label)
  )
  se <- cal_score(se,
    par.idf = list(label = "Group"),
    par.iae = list(label = "Group")
  )
  scored <- SummarizedExperiment::assay(se, "score")
  lab    <- SummarizedExperiment::colData(se)$Group

  ## closed form (default gaussian)
  betas_cf <- smartid:::fit_label_betas_closed_form(scored, factor(lab), NULL)
  ## glm apply-loop reference
  betas_gl <- smartid:::fit_label_betas_glm_loop(
    scored, factor(lab), NULL, stats::gaussian()
  )
  expect_equal(unname(betas_cf), unname(betas_gl), tolerance = 1e-8)
})

test_that("rows with zero SD collapse to zero after row_scale_zmean", {
  m <- matrix(rnorm(60), nrow = 6)
  m[3, ] <- 5                    # constant row -> SD 0
  m[5, ] <- NA_real_             # all-NA row -> SD NA
  rownames(m) <- paste0("g", seq_len(6))
  scaled <- smartid:::row_scale_zmean(m)
  expect_true(all(scaled[3, ] == 0))
  expect_true(all(scaled[5, ] == 0))
})

test_that("pmax0_offset preserves sparsity and masks negatives", {
  x <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(-1, 2, 5), dims = c(4, 4)
  )
  ## thres == 0 short-circuits to identity
  expect_identical(smartid:::pmax0_offset(x, 0), x)
  ## thres == 3 zeroes out 2 (< 3) and clips 5 to 2
  out <- smartid:::pmax0_offset(x, 3)
  expect_s4_class(out, "CsparseMatrix")
  expect_equal(as.numeric(out[2, 2]), 0)
  expect_equal(as.numeric(out[3, 3]), 2)
  expect_equal(as.numeric(out[1, 1]), 0)
})

test_that("rowwise_notin_max vectorised form matches naive apply", {
  set.seed(7)
  mat <- matrix(runif(40, 0, 10), nrow = 8)
  colnames(mat) <- paste0("k", seq_len(ncol(mat)))
  fast <- smartid:::rowwise_notin_max(mat)
  slow <- vapply(
    colnames(mat), function(type)
      apply(mat, 1, function(x) max(x[names(x) != type])),
    numeric(nrow(mat))
  )
  expect_equal(unname(fast), unname(slow), tolerance = 1e-12)
})
