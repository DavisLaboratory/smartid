% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scale_mgm.R
\name{scale_mgm}
\alias{scale_mgm}
\title{scale by mean of group mean in case extreme unbalanced data}
\usage{
scale_mgm(expr, label)
}
\arguments{
\item{expr}{matrix}

\item{label}{a vector of group label}
}
\value{
scaled matrix
}
\description{
scale by mean of group mean in case extreme unbalanced data
}
\examples{
scale_mgm(matrix(rnorm(100), 10), label = rep(letters[1:2], 5))
}
