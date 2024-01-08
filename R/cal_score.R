#' Calculate scores on given features
#'
#' @param score matrix, features in row and samples in column
#' @param features vector, feature names to compute score
#'
#' @return a vector of score
#' @export
#'
#' @examples
#' data <- matrix(rnorm(100), 10, dimnames = list(1:10))
#' cal_score_init(data, 1:5)
cal_score_init <- function(score, features = NULL) {
  if(is.null(features)) features <- rownames(score)

  ## check features
  if(!all(features %in% rownames(score)))
    warning(sprintf("Feature %s is not in score!\n",
                    setdiff(features, rownames(score))))
  features <- intersect(rownames(score), features)
  stopifnot("less than 2 features are in score rows!" = length(features) > 1)

  ## calculate mean score of features
  m_score <- colMeans(score[features, ], na.rm = TRUE)
  return(m_score)
}
