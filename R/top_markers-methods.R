#' @include top_markers.R
NULL

#' @rdname top_markers
setMethod(
  "top_markers", signature(
    data = "AnyMatrix"
  ),
  function(data,
           label,
           n = 10,
           use.glm = TRUE,
           scale = TRUE,
           use.mgm = TRUE,
           softmax = TRUE,
           ...) {
    ## check
    stopifnot(
      "n must be an integer" = round(n) == n,
      "label must be atomic and the same length as data column number" =
        is.atomic(label) & length(label) == ncol(data),
      "use.glm must be logical" = is.logical(use.glm),
      "scale must be logical" = is.logical(scale),
      "use.mgm must be logical" = is.logical(use.mgm),
      "softmax must be logical" = is.logical(softmax)
    )

    top_m <- top_markers_init(
      data = data,
      label = label,
      n = n,
      use.glm = use.glm,
      scale = scale,
      use.mgm = use.mgm,
      softmax = softmax,
      ...
    )

    return(top_m)
  }
)

#' @rdname top_markers
setMethod(
  "top_markers", signature(
    data = "SummarizedExperiment"
  ),
  function(data,
           label,
           n = 10,
           use.glm = TRUE,
           scale = TRUE,
           use.mgm = TRUE,
           softmax = TRUE,
           slot = "score",
           ...) {
    ## check
    stopifnot(
      "label must be a single character for se object" =
        is.character(label) & length(label) == 1
    )

    ## get expr
    expr <- SummarizedExperiment::assay(data, i = slot)
    ## get label
    label <- colData(data)[[label]]

    top_m <- top_markers(
      data = expr,
      label = label,
      n = n,
      use.glm = use.glm,
      scale = scale,
      use.mgm = use.mgm,
      softmax = softmax,
      ...
    )

    return(top_m)
  }
)
