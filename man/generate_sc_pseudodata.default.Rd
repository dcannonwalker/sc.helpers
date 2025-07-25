% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_sc_pseudodata.R
\name{generate_sc_pseudodata.default}
\alias{generate_sc_pseudodata.default}
\title{Generate pseudodata from a matrix of single cell counts}
\usage{
\method{generate_sc_pseudodata}{default}(
  y,
  group = NULL,
  lib.size = NULL,
  dispersion,
  prior.df = 10,
  tol = 1e-06,
  niter = 4,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{y}{A matrix of single cell counts}
}
\description{
Generate pseudodata from a matrix of single cell counts
}
