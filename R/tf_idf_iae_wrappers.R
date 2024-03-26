#################################################
#-----------------TF variants------------------#
#################################################

### $$\mathbf{TF_{i,j}}=\frac{N_{i,j}}{\sum_j{N_{i,j}}}$$

## term frequency
#' compute term/feature frequency within each cell
#'
#' @details
#' \deqn{\mathbf{TF_{i,j}}=\frac{N_{i,j}}{\sum_j{N_{i,j}}}}
#'
#' @param expr a count matrix, features in row and cells in column
#' @param log logical, if to do log-transformation
#'
#' @return a matrix of term/gene frequency
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::tf(data)
tf <- function(expr, log = FALSE) {
  t.f <- sweep(expr, 2, colSums(expr, na.rm = TRUE) + 0.01, FUN = "/")
  if(log) {
    t.f <- log1p(t.f)
  }

  return(t.f)
}

#################################################
#-----------------IDF variants------------------#
#################################################

##-----------------unlabeled-------------------##

## inverse document frequency
### $$\mathbf{IDF_i} = log(1+\frac{n}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' standard inverse cell frequency
#'
#' @details
#' \deqn{\mathbf{IDF_i} = log(1+\frac{n}{n_i+1})}
#'
#' @inheritParams idf_rf
#'
#' @return a vector of inverse cell frequency score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf(data)
idf <- function(expr, features = NULL, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr)  ## number of total obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  n_sub <- rowSums(expr[features, ] > thres) ## num of obs contain feature i > thres

  idf <- log1p(n_obs/(n_sub + 1))
  return(idf)
}

## inverse document frequency max
### $$\mathbf{IDF_{i,j}} = log(\frac{max_{\{i^{'}\in j\}}(n_{i^{'}})}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse cell frequency: max
#'
#' @details
#' \deqn{\mathbf{IDF_{i,j}} = log(\frac{max_{\{i^{'}\in j\}}(n_{i^{'}})}{n_i+1})}
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse cell frequency score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_m(data)
idf_m <- function(expr, features = NULL, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr)  ## number of total obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr, probs = 0.25, na.rm = TRUE)
  n_sub <- rowSums(expr > thres) ## num of obs contain feature i > thres

  n_max <- ifelse(expr > thres, n_sub, 0) |> sparseMatrixStats::colMaxs()

  idf <- matrix(1/(1 + n_sub), ncol = 1) %*% matrix(n_max, nrow = 1)
  dimnames(idf) <- dimnames(expr)

  idf <- log(idf[features, ])
  return(idf)
}

## inverse document frequency sd
### $$\mathbf{IDF_i} = log(1+sd(N_{i})*\frac{n}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse cell frequency using standard deviation (SD)
#'
#' @details
#' \deqn{\mathbf{IDF_i} = log(1+sd(N_{i})*\frac{n}{n_i+1})}
#'
#' @inheritParams idf_rf
#'
#' @return a vector of inverse cell frequency score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_sd(data)
idf_sd <- function(expr, features = NULL, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr)  ## number of total obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  n_sub <- rowSums(expr[features, ] > thres) ## num of obs contain feature i > thres
  sd_row <- sparseMatrixStats::rowSds(expr[features,], na.rm = TRUE)

  idf <- log1p(sd_row*n_obs/(n_sub + 1))
  return(idf)
}

## inverse document frequency using hdbscan cluster as label
#' inverse document frequency using hdbscan cluster as label
#'
#' @details
#' Details as [smartid:::idf_prob()].
#'
#' @inheritParams idf_rf
#' @param minPts integer, minimum size of clusters, default 2.
#'               Details in [dbscan::hdbscan()].
#' @param ... parameters for [dbscan::hdbscan()]
#'
#' @return a matrix of inverse cell frequency score
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_hdb(data)
idf_hdb <- function(expr, features = NULL, multi = TRUE,
                    thres = 0, minPts = 2, ...) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  ## initially compute naive tf-idf
  # tf <- (edgeR::cpm(expr)/1e6)[features, ]
  tf <- sweep(expr, 2, colSums(expr, na.rm = TRUE), "/")[features, ]
  tfidf <- tf*idf(expr, features = features, thres = thres)

  ## cluster obs based on given features
  cluster <- dbscan::hdbscan(t(tfidf), minPts = minPts, ...)$cluster
  ## factor cluster
  cluster <- factor(cluster)

  idf <- idf_prob(expr = expr, features = features,
                  label = cluster, multi = multi,
                  thres = thres)
  return(idf)
}

##------------------labeled--------------------##

## labeled inverse document frequency: relative frequency
### $$\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}})$$

#' labeled inverse cell frequency: relative frequency
#'
#' @details
#' \deqn{\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}})}
#'
#' @param expr a matrix, features in row and cells in column
#' @param features vector, feature names or indexes to compute
#' @param label vector, group label of each cell
#' @param multi logical, if to compute based on binary (FALSE) or multi-class (TRUE)
#' @param thres numeric, cell only counts when expr > threshold, default 0
#'
#' @return a matrix of inverse cell frequency score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_rf(data, label = sample(c("A", "B"), 10, replace = TRUE))
idf_rf <- function(expr, features = NULL, label,
                   multi = TRUE, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  df_n <- expr[features, ] > thres ## if contain feature i > thres

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- sapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(df_n[, label == type], na.rm = TRUE)
  })  ## mean counts for each gene in the group
  if(multi == TRUE) {
    mean_row_notin <- sapply(colnames(mean_row_in), function(type) {
      apply(mean_row_in, 1, function(x) max(x[names(x) != type]))
    })  ## doc freq for each gene not in group for multi-class: max(mean(N in Gi))
  }else {
    mean_row_notin <- sapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(df_n[, label != type], na.rm = TRUE)
    })  ## doc freq for each gene not in group for bi-class
  }

  idf <- log1p((mean_row_in/(mean_row_notin+1e-8))[,label])  ## IDF scores

  return(idf)
}

## labeled inverse document frequency: probability based
### $$\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}}\frac{n_{i,j\in D}}{n_{j\in D}})$$
### modified from
### $$\mathbf{IDF_{i,j}} = log(1+\frac{A}{max(B)+1}*\frac{A}{C+1})$$
### A denotes the number of cells belonging to category D where the gene i occurs at least once; B denotes the number of cells not belonging to category D where the gene i occurs at least once; C denotes the number of cells belonging to category D where the gene i does not occur; D denotes the number of cells not belonging to category D where the gene i does not occur.

#' labeled inverse cell frequency: probability based
#'
#' @details
#' \deqn{\mathbf{IDF_{i,j}} = log(1+\frac{\frac{n_{i,j\in D}}{n_{j\in D}}}{max(\frac{n_{i,j\in \hat D}}{n_{j\in \hat D}})+ e^{-8}}\frac{n_{i,j\in D}}{n_{j\in D}})}
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse cell frequency score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_prob(data, label = sample(c("A", "B"), 10, replace = TRUE))
idf_prob <- function(expr, features = NULL, label,
                     multi = TRUE, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  df_n <- expr[features, ] > thres ## if contain feature i > thres
  df_n_inv <- !df_n ## if not contain feature i > thres

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- sapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(df_n[, label == type], na.rm = TRUE)
  })  ## mean counts for each gene in the group
  if(multi == TRUE) {
    mean_row_notin <- sapply(colnames(mean_row_in), function(type) {
      apply(mean_row_in, 1, function(x) max(x[names(x) != type]))
    })  ## doc freq for each gene not in group for multi-class: max(mean(N in Gi))
  }else {
    mean_row_notin <- sapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(df_n[, label != type], na.rm = TRUE)
    })  ## doc freq for each gene not in group for bi-class
  }

  idf <- log1p((mean_row_in^2/(mean_row_notin+1e-8))[,label])  ## IDF scores

  return(idf)
}

## labeled inverse document frequency IGM
### $$\mathbf{IGM_i} = log(1+\lambda\frac{max(n_{i,j\in D})_{k}}{\sum_{k}^{K}((n_{i,j\in D})_{k}*r_{k})+e^{-8}})$$
### $$\mathbf{k}: type\ in\ total\ group\ K$$
### $$\mathbf{r_{k}}: rank\ of\ n_{i,j\in D})_{k}\ in\ total\ group\ K$$

#' labeled inverse cell frequency: IGM
#'
#' @details
#' \deqn{\mathbf{IGM_i} = log(1+\lambda\frac{max(n_{i,j\in D})_{k}}{\sum_{k}^{K}((n_{i,j\in D})_{k}*r_{k})+e^{-8}})}
#'
#' @inheritParams idf_rf
#' @param lambda numeric, hyperparameter for IGM
#'
#' @return a vector of inverse gravity moment score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::idf_igm(data, label = sample(c("A", "B"), 10, replace = TRUE))
idf_igm <- function(expr, features = NULL, label, lambda = 7, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  df_n <- expr[features, ] > thres ## if contain feature i > thres

  mean_row_in <- sapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(df_n[, label == type], na.rm = TRUE)
  })  ## mean counts for each gene in the group
  ## calculate igm for each gene
  igm <- apply(mean_row_in, 1, function(x) {
    max(x)/(sum(x*rank(-x), na.rm = TRUE)+1e-8)
  })

  idf <- log1p(lambda*igm)  ## TF-IDF scores
  return(idf)
}

#################################################
#-----------------IAE variants------------------#
#################################################

##-----------------unlabeled-------------------##

## inverse average expression
### $$\mathbf{IAE_i} = log(1+\frac{n}{\sum_j^n\hat N_{i,j}+1})$$
### $$n: total\ counts\ of\ documents;\ \hat N_{i,j}: max(0, N_{i,j} - threshold)$$

#' standard inverse average expression
#'
#' @details
#' \deqn{\mathbf{IAE_i} = log(1+\frac{n}{\hat N_{i,j}+1})}
#'
#' @inheritParams idf_rf
#'
#' @return a vector of inverse average expression score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae(data)
iae <- function(expr, features = NULL, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr)  ## number of total obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  expr_offset <- expr[features, ] - thres ## subtract offset
  expr_offset[expr_offset < 0] <- 0
  s_row <- rowSums(expr_offset) ## total counts of feature i across all cells

  iae <- log1p(n_obs/(s_row + 1))
  return(iae)
}

## inverse average expression max
### $$\mathbf{IAE_{i,d}} = log(1+\frac{max_{\{i^{'}\in d\}}(n_{i^{'}})}{n_i+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse average expression: max
#'
#' @details
#' \deqn{\mathbf{IAE_{i,j}} = log(1+\frac{max_{\{i^{'}\in j\}}(n_{i^{'}})}{\sum_{j = 1}^{n} sign(N_{i,j} > threshold)+1})}
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse average expression score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_m(data)
iae_m <- function(expr, features = NULL, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr)  ## number of total obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr, probs = 0.25, na.rm = TRUE)
  expr_offset <- expr - thres ## subtract offset
  expr_offset[expr_offset < 0] <- 0
  s_row <- rowSums(expr_offset) ## total counts of feature i across all cells

  s_max <- ifelse(expr_offset > 0, s_row, 0) |> sparseMatrixStats::colMaxs()

  iae <- matrix(1/(1 + s_row), ncol = 1) %*% matrix(s_max, nrow = 1)
  dimnames(iae) <- dimnames(expr)
  iae <- log1p(iae[features, ])
  return(iae)
}

## inverse average expression sd
### $$\mathbf{IAE} = log(1+sd(N_{i})*\frac{n}{\sum_{j=1}^{n}N_{i,j}+1})$$
### $$n: total\ counts\ of\ documents;\ n_i: \sum_{j = 1}^{n} sign(N_{i,j} > threshold)$$

#' inverse average expression using standard deviation (SD)
#'
#' @details
#' \deqn{\mathbf{IAE} = log(1+sd(N_{i})*\frac{n}{\sum_{j=1}^{n}N_{i,j}+1})}
#'
#' @inheritParams idf_rf
#'
#' @return a vector of inverse average expression score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_sd(data)
iae_sd <- function(expr, features = NULL, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))
  n_obs <- ncol(expr)  ## number of obs

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  expr_offset <- expr[features, ] - thres ## subtract offset
  expr_offset[expr_offset < 0] <- 0
  s_row <- sparseMatrixStats::rowSums2(expr_offset, na.rm = TRUE)  ## summed counts for each gene
  sd_row <- sparseMatrixStats::rowSds(expr_offset, na.rm = TRUE)

  iae <- log1p(sd_row*n_obs/(s_row + 1))  ## IDF scores
  return(iae)
}

## inverse average expression using hdbscan cluster as label
#' inverse average expression using hdbscan cluster as label
#'
#' @details
#' Details as [smartid:::iae_prob()].
#'
#' @inheritParams idf_rf
#' @param minPts integer, minimum size of clusters, default 2.
#'               Details in [dbscan::hdbscan()].
#' @param ... parameters for [dbscan::hdbscan()]
#'
#' @return a matrix of inverse average expression score
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_hdb(data)
iae_hdb <- function(expr, features = NULL, multi = TRUE,
                    thres = 0, minPts = 2, ...) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  ## initially compute naive tf-idf
  # tf <- (edgeR::cpm(expr)/1e6)[features, ]
  tf <- sweep(expr, 2, colSums(expr, na.rm = TRUE), "/")[features, ]
  tfidf <- tf*iae(expr, features = features, thres = thres)

  ## cluster obs based on given features
  cluster <- dbscan::hdbscan(t(tfidf), minPts = minPts, ...)$cluster
  ## factor cluster
  cluster <- factor(cluster)

  # # thres <- 0
  # # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  # expr_offset <- expr[features, ] - thres ## subtract offset
  # expr_offset[expr_offset < 0] <- 0
  #
  # mean_row_in <- sapply(levels(cluster), function(type) {
  #   rowMeans(expr_offset[, cluster == type], na.rm = TRUE)
  # }) |> setNames(levels(cluster))  ## mean counts for each gene in the group
  # mean_row_notin <- sapply(levels(cluster), function(type) {
  #   apply(mean_row_in, 1, function(x) max(x[names(x) != type]))
  # }) |> setNames(levels(cluster))  ## mean counts for each gene not in group
  #
  # iae <- log1p((mean_row_in/(mean_row_notin+0.01))[, cluster])  ## IDF scores

  iae <- iae_prob(expr = expr, features = features,
                  label = cluster, multi = multi,
                  thres = thres)
  return(iae)
}

##------------------labeled--------------------##

## labeled inverse average expression: relative frequency
### $$\mathbf{IAE} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}})$$
###   or
### $$\mathbf{IAE} = log(1+\frac{mean(N_{i,j\in D})}{mean(N_{i,j\notin D})+ e^{-8}})$$

#' labeled inverse average expression: relative frequency
#'
#' @details
#' \deqn{\mathbf{IAE} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}})}
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse average expression score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_rf(data, label = sample(c("A", "B"), 10, replace = TRUE))
iae_rf <- function(expr, features = NULL, label,
                   multi = TRUE, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  expr_offset <- expr[features, ] - thres ## subtract offset
  expr_offset[expr_offset < 0] <- 0

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- sapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(expr_offset[, label == type], na.rm = TRUE)
  })  ## mean counts for each gene in the group
  if(multi == TRUE) {
    mean_row_notin <- sapply(colnames(mean_row_in), function(type) {
      apply(mean_row_in, 1, function(x) max(x[names(x) != type]))
    })  ## mean counts for each gene not in group for multi-class: max(mean(N in Gi))
  }else {
    mean_row_notin <- sapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(expr_offset[, label != type], na.rm = TRUE)
    })  ## mean counts for each gene not in group for bi-class
  }

  iae <- log1p((mean_row_in/(mean_row_notin+1e-8))[,label])  ## IDF scores
  return(iae)
}

## labeled inverse average expression: probability based
### $$\mathbf{IAE_{i,j}} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}}*mean(N_{i,j\in D}))$$
### modified from
### $$\mathbf{IDF_{i,j}} = log(1+\frac{A}{max(B)+1}*\frac{A}{C+1})$$
### A denotes the number of cells belonging to category D where the gene i occurs at least once; B denotes the number of cells not belonging to category D where the gene i occurs at least once; C denotes the number of cells belonging to category D where the gene i does not occur; D denotes the number of cells not belonging to category D where the gene i does not occur.

#' labeled inverse average expression: probability based
#'
#' @details
#' \deqn{\mathbf{IAE_{i,j}} = log(1+\frac{mean(N_{i,j\in D})}{max(mean(N_{i,j\in \hat D}))+ e^{-8}}*mean(N_{i,j\in D}))}
#'
#' @inheritParams idf_rf
#'
#' @return a matrix of inverse average expression score
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_prob(data, label = sample(c("A", "B"), 10, replace = TRUE))
iae_prob <- function(expr, features = NULL, label,
                     multi = TRUE, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  expr_offset <- expr[features, ] - thres ## subtract offset
  expr_offset[expr_offset < 0] <- 0

  ## convert label into character in case problem for factor
  label <- as.character(label)
  mean_row_in <- sapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(expr_offset[, label == type], na.rm = TRUE)
  })  ## mean counts for each gene in the group
  if(multi == TRUE) {
    mean_row_notin <- sapply(colnames(mean_row_in), function(type) {
      apply(mean_row_in, 1, function(x) max(x[names(x) != type]))
    })  ## mean counts for each gene not in group for multi-class: max(mean(N in Gi))
  }else {
    mean_row_notin <- sapply(unique(label), function(type) {
      sparseMatrixStats::rowMeans2(expr_offset[, label != type], na.rm = TRUE)
    })  ## mean counts for each gene not in group for bi-class
  }

  iae <- log1p((mean_row_in^2/(mean_row_notin+1e-8))[,label])  ## IDF scores
  return(iae)
}

## labeled inverse average expression IGM
### $$\mathbf{IGM_i} = log(1+\lambda\frac{max(mean(N_{i,j\in D})_{k})}{\sum_{k}^{K}(mean(N_{i,j\in D})_{k}*r_{k})+e^{-8}})$$
### $$\mathbf{k}: type\ in\ total\ group\ K$$
### $$\mathbf{r_{k}}: rank\ of\ mean(N_{i,j\in D})_{k}\ in\ total\ group\ K$$

#' labeled inverse average expression: IGM
#'
#' @details
#' \deqn{\mathbf{IGM_i} = log(1+\lambda\frac{max(mean(N_{i,j\in D})_{k})}{\sum_{k}^{K}(mean(N_{i,j\in D})_{k}*r_{k})+e^{-8}})}
#'
#' @inheritParams idf_rf
#' @param lambda numeric, hyperparameter for IGM
#'
#' @return a vector of inverse gravity moment score for each feature
#'
#' @examples
#' data <- matrix(rpois(100, 2), 10, dimnames = list(1:10))
#' smartid:::iae_igm(data, label = sample(c("A", "B"), 10, replace = TRUE))
iae_igm <- function(expr, features = NULL, label, lambda = 7, thres = 0) {
  if(is.null(features)) features <- seq_len(nrow(expr))

  # thres <- 0
  # thres <- sparseMatrixStats::rowQuantiles(expr[features, ], probs = 0.25, na.rm = TRUE)
  expr_offset <- expr[features, ] - thres ## subtract offset
  expr_offset[expr_offset < 0] <- 0

  mean_row_in <- sapply(unique(label), function(type) {
    sparseMatrixStats::rowMeans2(expr_offset[, label == type], na.rm = TRUE)
  })  ## mean counts for each gene in the group
  ## calculate igm for each gene
  igm <- apply(mean_row_in, 1, function(x) {
    max(x)/(sum(x*rank(-x), na.rm = TRUE)+1e-8)
  })

  iae <- log1p(lambda*igm)  ## TF-IDF scores
  return(iae)
}

