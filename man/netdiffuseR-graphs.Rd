% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_data.r
\name{netdiffuseR-graphs}
\alias{netdiffuseR-graphs}
\title{Network data formats}
\value{
No return value (this manual entry only provides information).
}
\description{
List of accepted graph formats
}
\details{
The \pkg{netdiffuseR} package can handle different types of graph
objects. Two general classes are defined across the package's functions:
static graphs, and dynamic graphs.
\itemize{
 \item{In the case of \strong{static graphs}, these are represented as adjacency
 matrices of size \eqn{n\times n}{n * n} and can be either \code{\link{matrix}}
 (dense matrices) or \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
 (sparse matrix from the \pkg{\link[Matrix:Matrix]{Matrix}} package). While
 most of the package functions are defined for both classes, the default output
 graph is sparse, i.e. \code{dgCMatrix}.}
 \item{With respect to \strong{dynamic graphs}, these are represented by either
 a \code{\link{diffnet}} object, an \code{\link{array}} of size
 \eqn{n\times n \times T}{n * n * T}, or a list of size \eqn{T}
 with sparse matrices (class \code{dgCMatrix}) of size \eqn{n\times n}{n * n}.
 Just like the static graph case, while most of the functions accept both
 graph types, the default output is \code{dgCMatrix}.}
}
}
\section{diffnet objects}{

 In the case of \code{diffnet}-class objects, the following arguments can be omitted
 when calling fuictions suitable for graph objects:
 \itemize{
   \item{\code{toa}: Time of Adoption vector}
   \item{\code{adopt}: Adoption Matrix}
   \item{\code{cumadopt}: Cumulative Adoption Matrix}
   \item{\code{undirected}: Whether the graph is directed or not}
 }
}

\section{Objects' names}{

When possible, \pkg{netdiffuseR} will try to reuse graphs dimensional names,
this is, \code{\link{rownames}}, \code{\link{colnames}}, \code{\link{dimnames}}
and \code{\link{names}} (in the case of dynamic graphs as lists). Otherwise,
when no names are provided, these will be created from scratch.
}

\author{
George G. Vega Yon
}
\concept{graph formats}
