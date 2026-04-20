# Capture legacy outputs of cal_score() and top_markers() from the current
# pre-refactor HEAD. Run ONCE from the package root with the source branch
# loaded via devtools::load_all(). The resulting .rds file is committed and
# consumed by tests/testthat/test-numerical-equivalence.R to verify the
# memory-optimization refactor preserves numerical equivalence.
#
# Usage (from package root):
#   Rscript tests/testthat/testdata/capture_legacy_snapshots.R
#
# Prerequisites: devtools, Matrix, SummarizedExperiment, S4Vectors installed.

suppressPackageStartupMessages({
  library(devtools)
  library(Matrix)
  library(SummarizedExperiment)
  library(S4Vectors)
})

## Load source package at current HEAD (pre-refactor)
devtools::load_all(".", quiet = TRUE)

## Deterministic test inputs ----------------------------------------------
set.seed(2026)

G <- 40L   # features
N <- 60L   # cells
K <- 3L    # groups

counts_dense <- matrix(rpois(G * N, lambda = 1.5), nrow = G)
rownames(counts_dense) <- paste0("gene", seq_len(G))
colnames(counts_dense) <- paste0("cell", seq_len(N))

counts_sparse <- as(counts_dense, "CsparseMatrix")  # dgCMatrix

label <- factor(rep(LETTERS[seq_len(K)], length.out = N))
batch <- factor(rep(c("b1", "b2"), length.out = N))

se_dense <- SummarizedExperiment(
  assays  = list(counts = counts_dense),
  colData = DataFrame(Group = label, Batch = batch)
)

## Reference outputs ------------------------------------------------------

## (1) cal_score on dense matrix, default idf/iae = "prob", logtf
score_mat_dense <- cal_score(
  data    = counts_dense,
  tf      = "logtf",
  idf     = "prob",
  iae     = "prob",
  par.idf = list(label = as.character(label)),
  par.iae = list(label = as.character(label))
)

## (2) cal_score on SummarizedExperiment (currently stores intermediates
##     in metadata; post-refactor we will compare via return.intermediate)
se_score <- cal_score(
  data    = se_dense,
  tf      = "logtf",
  idf     = "prob",
  iae     = "prob",
  par.idf = list(label = "Group"),
  par.iae = list(label = "Group")
)

## (3) top_markers with use.glm = TRUE (gaussian default)
tm_glm <- top_markers(
  data    = se_score,
  label   = "Group",
  n       = 10L,
  use.glm = TRUE,
  slot    = "score"
)

## (4) top_markers with use.glm = TRUE + batch
tm_glm_batch <- top_markers(
  data    = se_score,
  label   = "Group",
  n       = 10L,
  use.glm = TRUE,
  batch   = "Batch",
  slot    = "score"
)

## (5) top_markers with use.glm = FALSE, method = "mean"
tm_abs_mean <- top_markers(
  data    = se_score,
  label   = "Group",
  n       = 10L,
  use.glm = FALSE,
  method  = "mean",
  slot    = "score"
)

## (6) top_markers with use.glm = FALSE, method = "median"
tm_abs_median <- top_markers(
  data    = se_score,
  label   = "Group",
  n       = 10L,
  use.glm = FALSE,
  method  = "median",
  slot    = "score"
)

## Assemble & save --------------------------------------------------------
snap <- list(
  inputs = list(
    counts_dense  = counts_dense,
    counts_sparse = counts_sparse,
    label         = label,
    batch         = batch
  ),
  cal_score = list(
    matrix_out = score_mat_dense,        # list(score = ..., tf, idf, iae)
    se_assay   = as.matrix(assay(se_score, "score")),
    se_tf      = as.matrix(metadata(se_score)$tf),
    se_idf     = as.matrix(metadata(se_score)$idf),
    se_iae     = as.matrix(metadata(se_score)$iae)
  ),
  top_markers = list(
    glm        = tm_glm,
    glm_batch  = tm_glm_batch,
    abs_mean   = tm_abs_mean,
    abs_median = tm_abs_median
  ),
  meta = list(
    captured_at = Sys.time(),
    R_version   = R.version.string,
    smartid_ver = as.character(packageVersion("smartid"))
  )
)

out <- "tests/testthat/testdata/legacy_scores.rds"
saveRDS(snap, file = out, version = 2L)

message("Legacy snapshot written to: ", out)
message("  cal_score$score dim:  ", paste(dim(snap$cal_score$matrix_out$score), collapse = " x "))
message("  SE assay 'score' dim: ", paste(dim(snap$cal_score$se_assay), collapse = " x "))
message("  top_markers$glm rows: ", nrow(snap$top_markers$glm))
