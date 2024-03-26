## =========================================================================##
## =========================================================================##
##                       Generic methods definitions                        ##
## =========================================================================##
## =========================================================================##


## ===========================================================================
## Generic for calculating score
## ---------------------------------------------------------------------------
#' @title calculate combined score
#'
#' @description compute TF (term/feature frequency), IDF (inverse document/cell
#'     frequency), IAE (inverse average expression of features) and combine the
#'     the final score
#'
#' @inheritParams cal_score_init
#' @param data an expression object, can be matrix or SummarizedExperiment
#' @param slot a character, specify which slot to use when data is se object,
#'     optional, default 'counts'
#' @param new.slot a character, specify the name of slot to save score in se object,
#'     optional, default 'score'
#'
#' @return A list of matrices or se object containing combined score
#'
#' @include score.R
#' @export
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' cal_score(
#'   data,
#'   par.idf = list(label = sample(c("A", "B"), 10, replace = TRUE)),
#'   par.iae = list(label = sample(c("A", "B"), 10, replace = TRUE))
#' )
setGeneric(
  "cal_score",
  function(data,
           tf = c("logtf", "tf"),
           idf = "prob",
           iae = "prob",
           slot = "counts",
           new.slot = "score",
           par.idf = NULL,
           par.iae = NULL) {
    standardGeneric("cal_score")
  }
)

## ===========================================================================
## Generic for scaling score and get top markers
## ---------------------------------------------------------------------------
#' @title scale score and return top markers
#'
#' @description scale and transform score and output top markers for groups
#'
#' @inheritParams top_markers_init
#' @param data an expression object, can be matrix or SummarizedExperiment
#' @param slot a character, specify which slot to use when data is se object,
#'     optional, default 'score'
#'
#' @return A tibble with top n feature names, group labels and ordered scores
#'
#' @include score.R
#' @export
#'
#' @examples
#' data <- matrix(rgamma(100, 2), 10, dimnames = list(1:10))
#' top_markers(data, label = rep(c("A", "B"), 5))
setGeneric(
  "top_markers",
  function(data,
           label,
           n = 10,
           use.glm = TRUE,
           scale = TRUE,
           use.mgm = TRUE,
           softmax = TRUE,
           slot = "score",
           ...) {
    standardGeneric("top_markers")
  }
)

## ===========================================================================
## Generic for compute overall score on given marker list
## ---------------------------------------------------------------------------
#' @title compute overall score based on the given marker list
#'
#' @description compute overall score based on the given marker list
#'
#' @inheritParams gs_score_init
#' @param data an expression object, can be matrix or SummarizedExperiment
#' @param features vector or named list, feature names to compute score
#' @param slot a character, specify which slot to use when data is se object,
#'     optional, default 'score'
#' @param suffix a character, specify the name suffix to save score when
#'     features is a named list
#'
#' @return A vector of overall score for each sample
#'
#' @include score.R
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(seq_len(10)))
#' gs_score(data, features = seq_len(3))
setGeneric(
  "gs_score",
  function(data,
           features = NULL,
           slot = "score",
           suffix = "score") {
    standardGeneric("gs_score")
  }
)
