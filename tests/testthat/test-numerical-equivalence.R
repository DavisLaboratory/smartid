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
    colData = S4Vectors::DataFrame(
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
  expect_null(S4Vectors::metadata(out)$tf)
  expect_null(S4Vectors::metadata(out)$idf)
  expect_null(S4Vectors::metadata(out)$iae)
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
    colData = S4Vectors::DataFrame(
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
  md <- S4Vectors::metadata(out)
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
    colData = S4Vectors::DataFrame(
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
    colData = S4Vectors::DataFrame(
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
    colData = S4Vectors::DataFrame(
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
    colData = S4Vectors::DataFrame(
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
