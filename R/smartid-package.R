#' @import ggplot2
#' @import methods
#' @import stats
#' @import graphics
#' @importClassesFrom Matrix dgCMatrix
NULL

#' Scoring and Marker Selection method based on modified TF-IDF
#'
#' `smartid` This package enables automated selection of group specific signature,
#'   especially for rare population. The package is developed for generating
#'   specifc lists of signature genes based on TF-IDF modified methods. It can
#'   also be used as a new gene-set scoring method or data transformation method.
#'   Multiple visualization functions are implemented in this package.
#'
#' @author Jinjin Chen \email{chen.j@@wehi.edu.au}
#' @name smartid_Package
#' @docType package
#' @aliases smartid smartid_Package
#' @keywords internal
#' @return Marker list and scores
#'
NULL

#' scRNA-seq test data of 4 groups simulated by `splatter`.
#'
#' A SingleCellExperiment object containing 4 groups with each group up-regulated
#' DEGs saved in metadata.
#'
#' @format A SingleCellExperiment object of 100genes * 400 cells.
#' @usage data(sim_sce_test)
#' @return SingleCellExperiment
#' @source [splatter::splatSimulate()]
"sim_sce_test"
