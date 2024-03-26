#' @include score.R
NULL

#' @rdname gs_score
setMethod(
  "gs_score", signature(
    data = "AnyMatrix"
  ),
  function(data,
           features = NULL) {
    ## compute overall score
    score <- gs_score_init(score = data, features = features)

    return(score)
  }
)

#' @rdname gs_score
setMethod(
  "gs_score", signature(
    data = "AnyMatrix",
    features = "list"
  ),
  function(data,
           features = NULL,
           suffix = "score") {
    ## check
    stopifnot(
      "suffix must be a character" = is.character(suffix) & length(suffix) == 1
    )

    ## compute overall score
    score <- vapply(
      names(features), \(i)
      gs_score(data = data, features = features[[i]]),
      rep(1, ncol(data))
    ) |>
      data.frame()
    ## set colnames
    colnames(score) <- paste(names(features), suffix, sep = ".")

    return(score)
  }
)

#' @rdname gs_score
setMethod(
  "gs_score", signature(
    data = "SummarizedExperiment"
  ),
  function(data,
           features = NULL,
           slot = "score",
           suffix = "score") {
    ## check
    stopifnot(
      "slot must be one character" = is.character(slot) & length(slot) == 1,
      "suffix must be a character" = is.character(suffix) & length(suffix) == 1
    )

    ## get expr
    expr <- SummarizedExperiment::assay(data, i = slot)

    ## compute score
    score <- gs_score(data = expr, features = features, suffix = suffix)

    SummarizedExperiment::colData(data) <- cbind(
      SummarizedExperiment::colData(data),
      score
    )

    return(data)
  }
)
