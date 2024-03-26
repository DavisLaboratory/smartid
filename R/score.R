#' @include tf_idf_iae_wrappers.R
NULL

#' Calculate scores of each cell on given features
#'
#' @param score matrix, features in row and samples in column
#' @param features vector, feature names to compute score
#'
#' @return a vector of score
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' gs_score_init(data, 1:5)
gs_score_init <- function(score, features = NULL) {
  if (is.null(features)) features <- rownames(score)

  ## check features
  if (!all(features %in% rownames(score))) {
    warning(sprintf(
      "Feature %s is not in score!\n",
      setdiff(features, rownames(score))
    ))
  }
  features <- intersect(rownames(score), features)
  stopifnot("less than 2 features are in score rows!" = length(features) > 1)

  ## calculate mean score of features
  m_score <- colMeans(score[features, , drop = FALSE], na.rm = TRUE)
  return(m_score)
}

#' @title Get names of available IDF and IAE methods
#'
#' @description Returns a named vector of IDF/IAE methods
#' @return names of methods implemented
#' @export
#'
#' @examples
#' idf_iae_methods()
idf_iae_methods <- function() {
  return(sort(c(
    "label probability" = "prob", "label relative frequency" = "rf",
    "label IGM" = "igm", "null" = "null",
    "unlabel max" = "m", "unlabel SD" = "sd",
    "unlabel HDBSCAN" = "hdb", "unlabel standard" = "standard"
  )))
}

#' Calculate score for each feature in each cell
#'
#' @param expr a count matrix, features in row and cells in column
#' @param tf a character, specify the TF method to use, can be "tf" or "logtf"
#' @param idf a character, specify the IDF method to use. Available methods can
#'     be accessed using [idf_iae_methods()]
#' @param iae a character, specify the IAE method to use. Available methods can
#'     be accessed using [idf_iae_methods()]
#' @param par.idf other parameters for specified IDF methods
#' @param par.iae other parameters for specified IAE methods
#'
#' @return a list of combined score, tf, idf and iae
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' label <- sample(c("A", "B"), 10, replace = TRUE)
#' smartid:::cal_score_init(data,
#'   par.idf = list(label = label),
#'   par.iae = list(label = label)
#' )
cal_score_init <- function(expr, tf = c("logtf", "tf"),
                           idf = "prob", iae = "prob",
                           par.idf = NULL, par.iae = NULL) {
  ## check
  tf <- match.arg(tf)
  idf <- match.arg(idf, choices = idf_iae_methods())
  iae <- match.arg(iae, choices = idf_iae_methods())
  stopifnot(
    "par.idf must be a named list or NULL" = is.null(par.idf) | is.list(par.idf),
    "par.iae must be a named list or NULL" = is.null(par.iae) | is.list(par.iae)
  )

  ## compute tf
  tf <- tf(expr, log = (tf == "logtf"))

  ## compute idf
  if (idf == "null") {
    idf <- 1
  } else {
    idf <- ifelse(idf == "standard", "idf", paste0("idf_", idf))
    idf <- do.call(idf, c(list(expr = expr), par.idf))
  }

  ## compute iae
  if (iae == "null") {
    iae <- 1
  } else {
    iae <- ifelse(iae == "standard", "iae", paste0("iae_", iae))
    iae <- do.call(iae, c(list(expr = expr), par.iae))
  }

  ## combined score
  score <- tf * idf * iae

  return(list(score = score, tf = tf, idf = idf, iae = iae))
}
