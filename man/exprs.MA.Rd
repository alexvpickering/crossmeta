% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_agil.R
\name{exprs.MA}
\alias{exprs.MA}
\title{Extract Log-Expression Matrix from MAList}
\usage{
exprs.MA(MA)
}
\arguments{
\item{MA}{an \code{MAList} object.}
}
\value{
A numeric matrix with twice the columns of the input.
}
\description{
Converts M and A-values to log-expression values. The output matrix will have two columns for each array, in the order all red then all green.
Adapted from \link[limma]{plotDensities.MAList} instead of \link[limma]{exprs.MA} so that order is same as \link{phenoData.ch2}.
}
