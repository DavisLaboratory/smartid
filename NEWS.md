# smartid 0.99.0

* Added a `NEWS.md` file to track changes to the package.

# smartid 0.99.1

* Ready for submission to Bioconductor.

# smartid 0.99.2

* Added test for `gs_score()` function.

# smartid 0.99.3

* Bump R version dependency to >= 4.4 and add details for TF, IDF, IAE functions.

# smartid 0.99.4

* Add details for TF, IDF, IAE functions.

# smartid 0.99.5

* Update `scale_mgm()` function adding pooled SD option, add details for scale function.

# smartid 1.1.1

* Update `cal_score()` function to convert input sparse matrix into dense matrix.

# smartid 1.1.2

* Update marker selection functions to fix wrong names of marker list.

# smartid 1.3.1

* Update top_markers function to allow batch correction for glm method.

# smartid 1.3.2

* Update batch param in top_markers function.

# smartid 1.7.1

* Add bioRxiv citation.

# smartid 1.7.2

* Fix dplyr defunct.

# smartid 1.7.3

* Major memory and performance optimization for `cal_score()` and
  `top_markers()`. On a 20,000 gene x 100,000 cell sparse input peak
  memory drops from roughly 100 GB to a few GB, and `top_markers()`
  with the default `gaussian()` family runs in seconds instead of
  hours.
* **Breaking change**: `cal_score()` no longer stores the intermediate
  `tf`, `idf` and `iae` matrices in `metadata()` by default. Callers
  that relied on `metadata(se)$tf / $idf / $iae` (introduced in v1.1.1)
  must now pass `return.intermediate = TRUE`. When the flag is `TRUE`,
  the stored `idf` / `iae` for labelled methods (`prob`, `rf`) are now
  compact `G x K` matrices (columns = unique labels); expand with
  `md$idf[, as.character(label)]` to recover the legacy per-cell form.
* Internal refactor: labelled `idf_prob`, `idf_rf`, `iae_prob` and
  `iae_rf` helpers now return a compact `G x K` matrix. `cal_score()`
  composes the final score through per-group column-block
  multiplication, avoiding the materialisation of full `G x N`
  intermediates.
* `cal_score()` no longer forces dense conversion of `dgCMatrix`
  inputs; the score assay stays sparse throughout the pipeline when
  the input is sparse.
* `tf()`, `idf_hdb()`, `iae_hdb()` and all IAE helpers now preserve
  `dgCMatrix` sparsity by routing column scaling through
  `Matrix::Diagonal` and replacing the densifying
  `x[x < 0] <- 0` pattern with `pmax0_offset()`.
* `top_markers_glm()` has a vectorised closed-form least-squares fast
  path for the default `gaussian()` + identity link; non-gaussian
  families or rank-deficient designs automatically fall back to the
  legacy per-gene `glm()` loop with no behaviour change.
* `top_markers_abs()` aggregates directly on the scored matrix via
  `sparseMatrixStats::rowMeans2 / rowMedians / rowMads`, removing the
  intermediate wide `data.frame` that previously reached tens of GB.
* `scale_mgm()` caches per-group column indices and collapses the
  two-step `(expr - mgm) / (sds + 1e-8)` into a single broadcast.
* The `multi = TRUE` branch of the labelled IDF/IAE helpers switched
  from an O(G * K^2) `apply()` to an O(G * K) top-1 + top-2 trick via
  the new `rowwise_notin_max()` helper.
* New `inst/bench/benchmark_smartid.R` micro-benchmark script; new
  `tests/testthat/test-numerical-equivalence.R` pins `cal_score()` and
  `top_markers()` outputs to a frozen pre-refactor snapshot at 1e-10
  tolerance.
* No new dependencies; the refactor relies entirely on `Matrix`,
  `sparseMatrixStats` and base R.
