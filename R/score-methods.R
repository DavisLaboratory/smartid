#' @include score.R smartid-package.R
NULL

#' @rdname cal_score
setMethod(
  "cal_score", signature(
    data = "AnyMatrix"
  ),
  function(data,
           tf = c("logtf", "tf"),
           idf = "prob",
           iae = "prob",
           par.idf = NULL,
           par.iae = NULL) {
    score <- cal_score_init(
      expr = as.matrix(data),
      tf = tf,
      idf = idf,
      iae = iae,
      par.idf = par.idf,
      par.iae = par.iae
    )

    return(score)
  }
)

#' @rdname cal_score
setMethod(
  "cal_score", signature(
    data = "SummarizedExperiment"
  ),
  function(data,
           tf = c("logtf", "tf"),
           idf = "prob",
           iae = "prob",
           slot = "counts",
           new.slot = "score",
           par.idf = NULL,
           par.iae = NULL) {
    ## get expr
    expr <- SummarizedExperiment::assay(data, i = slot)
    ## get label
    if (!is.null(par.idf) & !is.null(par.idf$label)) {
      par.idf$label <- SummarizedExperiment::colData(data)[[par.idf$label]]
    }
    if (!is.null(par.iae) & !is.null(par.iae$label)) {
      par.iae$label <- SummarizedExperiment::colData(data)[[par.iae$label]]
    }

    res <- cal_score(
      data = expr,
      tf = tf,
      idf = idf,
      iae = iae,
      par.idf = par.idf,
      par.iae = par.iae
    )

    SummarizedExperiment::assay(data, i = new.slot) <- res$score
    slot(data, "metadata")$tf <- res$tf
    slot(data, "metadata")$idf <- res$idf
    slot(data, "metadata")$iae <- res$iae

    return(data)
  }
)
