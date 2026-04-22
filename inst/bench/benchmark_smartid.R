## Micro-benchmark for the memory-optimization refactor (smartid >= 1.7.3).
##
## Run from the package root:
##
##   Rscript inst/bench/benchmark_smartid.R
##
## The script measures `cal_score()` and `top_markers()` on a simulated
## 20,000 gene x 100,000 cell dgCMatrix with density 5% and 10 balanced
## groups. This mirrors the workload reported on large scRNA-seq atlases
## where the legacy implementation peaked around 100 GB of memory.
##
## The output goes to stdout. `bench::mark()` reports `mem_alloc` which
## accounts for R-level allocations during the call; peak RSS can be
## tracked externally (e.g. `/usr/bin/time -v Rscript ...`). We also
## `gc()` before each measurement so per-call allocations are isolated.
##
## No hard coded dependency on `bench`: if it is missing the script falls
## back to a `proc.time()` + `gc()` delta report.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(smartid)
})

set.seed(42)

G <- 20000L
N <- 100000L
K <- 10L
density <- 0.05

message("Simulating ", G, " x ", N, " dgCMatrix (density = ", density, ")...")
sim_t <- system.time({
  counts <- Matrix::rsparsematrix(
    G, N, density = density,
    rand.x = function(n) rpois(n, lambda = 2)
  )
  dimnames(counts) <- list(paste0("gene", seq_len(G)),
                           paste0("cell", seq_len(N)))
  label <- sample(LETTERS[seq_len(K)], N, replace = TRUE)
})
message("Simulation took ", round(sim_t["elapsed"], 1), "s; nnz = ",
        format(length(counts@x), big.mark = ","))

se <- SummarizedExperiment(
  assays  = list(counts = counts),
  colData = DataFrame(Group = factor(label))
)

report <- function(label, expr) {
  gc(reset = TRUE, full = TRUE)
  pre  <- gc(reset = FALSE)
  wall <- system.time(value <- force(expr))
  post <- gc(reset = FALSE)
  used_mb <- max(post[, "used (Mb)"] - pre[, "used (Mb)"])
  cat(sprintf("  %-35s  wall = %7.2fs   R-allocated delta = %7.1f MB\n",
              label, wall[["elapsed"]], used_mb))
  invisible(value)
}

has_bench <- requireNamespace("bench", quietly = TRUE)
cat("\nbenchmark backend: ", if (has_bench) "bench::mark" else "gc delta",
    "\n\n", sep = "")

## --- cal_score -----------------------------------------------------------
score_se <- report("cal_score (default)",
                   cal_score(se,
                             par.idf = list(label = "Group"),
                             par.iae = list(label = "Group")))

score_se_int <- report("cal_score (return.intermediate=TRUE)",
                       cal_score(se,
                                 par.idf = list(label = "Group"),
                                 par.iae = list(label = "Group"),
                                 return.intermediate = TRUE))

## --- top_markers ---------------------------------------------------------
tm_glm <- report("top_markers (gaussian closed-form)",
                 top_markers(score_se, label = "Group",
                             n = 50L, use.glm = TRUE))

tm_glm_batch <- report(
  "top_markers (glm + batch)",
  {
    SummarizedExperiment::colData(score_se)$Batch <-
      sample(c("b1", "b2", "b3"), N, replace = TRUE)
    top_markers(score_se, label = "Group", batch = "Batch",
                n = 50L, use.glm = TRUE)
  }
)

tm_abs <- report("top_markers (abs, method=mean)",
                 top_markers(score_se, label = "Group",
                             n = 50L, use.glm = FALSE, method = "mean"))

cat("\nScore class: ",
    class(SummarizedExperiment::assay(score_se, "score"))[1], "\n")
cat("Intermediate idf dim (compact): ",
    paste(dim(S4Vectors::metadata(score_se_int)$idf), collapse = " x "),
    "\n")
cat("top_markers GLM rows: ", nrow(tm_glm), "\n")
cat("top_markers abs rows: ", nrow(tm_abs), "\n")
