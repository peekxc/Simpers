% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simpers.R
\name{simpers}
\alias{simpers}
\title{SimPers}
\usage{
simpers(x, sm, validate_maps = TRUE, elementary_mode = TRUE,
  dim = 2L, n_maps = NULL, thresh = 0,
  generators = elementary_mode, ...)
}
\arguments{
\item{x}{the file name for the input simplicial complex, or a simplex tree. See details.}

\item{sm}{input simplicial maps (either filename or character vector).}

\item{validate_maps}{whether to check to ensure the simplicial maps are 'valid'. Defaults to true. See details.}

\item{elementary_mode}{whether to run in elementary mode or general mode.}

\item{dim}{maximum dimension of the barcodes.}

\item{n_maps}{number of maps (general-mode only).}

\item{thresh}{barcode length threshold}

\item{...}{Other simpers parameters.}
}
\description{
Simple R-wrapper for the SimPers software for computing topological persistence under simplicial maps (see 1).
}
\details{
Other parameters are converted to their <option>=<value> string equivalents and passed as arguments
to the main options parser. \cr
\cr
If \code{validate_maps} and \code{elementary_mode} are TRUE, the simplicial maps are inspected prior to passing to Simpers.
Namely, its checked that the faces of the each simplex being inserted exist prior to insertion, and if not they are
preprended to the insertion at the same time. Elementary collapses that don't collapse to a non-existent vertex throw an exception.
All operations are changed such that their vertices are listed in increasing order.
Simplex ids cannot be changed; it is up to the user to ensure the indices of the operation are consistent across all simplicial maps.
}
\examples{
## Elementary mode example
sp_dir <- system.file("extdata", package = "Simpers")
cmd <- sprintf(" -d \%s -m \%s", file.path(sp_dir, "Elementary_Mode", "iDC.txt"), file.path(sp_dir, "Elementary_Mode", "operations.txt"))
}
\references{
1. Dey, Tamal K., Fengtao Fan, and Yusu Wang. "Computing topological persistence for simplicial maps." Proceedings of the thirtieth annual symposium on Computational geometry. ACM, 2014. \cr
2. \href{http://web.cse.ohio-state.edu/~dey.8/SimPers/Simpers.html}{SimPers homepage}
}
\author{
R-wrapper written by Matt Piekenbrock. Source code for computing the persistence copyright 2015 by Fengtao Fan and Dayu Shi, developed by the Jyamiti group at The Ohio State University
}
