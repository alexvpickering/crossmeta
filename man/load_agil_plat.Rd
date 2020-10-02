% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_agil.R
\name{load_agil_plat}
\alias{load_agil_plat}
\title{Load Agilent raw data}
\usage{
load_agil_plat(eset, ch2, gse_dir, gse_name, ensql)
}
\arguments{
\item{eset}{ExpressionSet from \link{getGEO}}

\item{ch2}{Boolean indicating in two-channel array}

\item{gse_dir}{Direction with Agilent raw data}

\item{gse_name}{Accession name for \code{eset}.}

\item{ensql}{For development. Path to sqlite file with ENTREZID and SYMBOL columns created in data-raw/entrezdt.}
}
\value{
ExpressionSet
}
\description{
Load Agilent raw data
}