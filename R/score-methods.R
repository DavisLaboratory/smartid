#' @include score.R
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
      expr = data,
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
      par.idf$label <- data@colData[[par.idf$label]]
    }
    if (!is.null(par.iae) & !is.null(par.iae$label)) {
      par.iae$label <- data@colData[[par.iae$label]]
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
    data@metadata$tf <- res$tf
    data@metadata$idf <- res$idf
    data@metadata$iae <- res$iae

    return(data)
  }
)
