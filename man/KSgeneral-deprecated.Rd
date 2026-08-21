\name{KSgeneral-deprecated}
\alias{KSgeneral-deprecated}
\alias{ks_c_cdf_Rcpp}
\title{Deprecated Functions in Package 'KSgeneral'}
\description{
These functions are provided for compatibility with earlier versions of
\pkg{KSgeneral} and may be removed in a future release.
}
\usage{
ks_c_cdf_Rcpp(n)
}
\arguments{
  \item{n}{A positive integer giving the sample size.}
}
\details{
\code{ks_c_cdf_Rcpp()} is deprecated as of \pkg{KSgeneral} 2.1.0.
It uses the historical file-based interface, reading boundary values from
\file{Boundary_Crossing_Time.txt} in the current working directory.

New code should use \code{\link{ks_c_cdf}(n, A, B)}, which accepts the lower
and upper boundary vectors directly and does not require intermediate file I/O.
}
\value{
A numeric value giving the complementary boundary-crossing probability.
}
\seealso{
\code{\link{ks_c_cdf}}, \code{\link[base]{Deprecated}}
}
