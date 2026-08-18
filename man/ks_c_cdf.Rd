\name{ks_c_cdf}
\alias{ks_c_cdf}
\title{Exact-KS-FFT Complementary CDF from Pre-computed Boundaries}
\description{
Computes the complementary cumulative distribution function used by the
one-sample Exact-KS-FFT method from pre-computed lower and upper boundary
vectors, without using an intermediate file.
}
\usage{
ks_c_cdf(n, A, B)
}
\arguments{
  \item{n}{A positive integer giving the sample size.}
  \item{A}{A numeric vector of length \code{n} containing the lower boundary
  values \eqn{A_i}, in nondecreasing order, with values in \eqn{[0,1]}.}
  \item{B}{A numeric vector of length \code{n} containing the upper boundary
  values \eqn{B_i}, in nondecreasing order, with values in \eqn{[0,1]}.}
}
\details{
The function evaluates the Exact-KS-FFT boundary-crossing calculation directly
from the vectors \eqn{A_i} and \eqn{B_i}. These vectors define the rectangular
region for the uniform order statistics,
\deqn{A_i \le U_{(i)} \le B_i, \qquad i=1,\ldots,n.}
It returns the complementary probability
\deqn{1 - P(A_i \le U_{(i)} \le B_i,\ i=1,\ldots,n).}

The boundary values can be constructed as in Steps 1 and 2 of the Exact-KS-FFT
method of Dimitrova, Kaishev and Tan (2020). For the standard two-sided
one-sample KS problem with a continuous null distribution and fixed \eqn{q},
one may use
\deqn{A_i = \max(0, i/n-q), \qquad
      B_i = \min(1, (i-1)/n+q).}

This interface replaces the historical package-internal mechanism that wrote
the boundary vectors to \file{Boundary_Crossing_Time.txt} before calling C++.
The vectors are now passed directly to the C++ implementation.

}
\value{
A numeric value giving the complementary probability associated with the
specified boundary vectors. For boundaries corresponding to a one-sample KS
statistic at threshold \eqn{q}, this is the Exact-KS-FFT p-value
\eqn{P(D_n \ge q)}.
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020)
"Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is
Purely Discrete, Mixed or Continuous". Journal of Statistical Software,
\bold{95}(10), 1--42. doi:10.18637/jss.v095.i10.

Moscovich A., Nadler B. (2017). "Fast Calculation of Boundary Crossing
Probabilities for Poisson Processes". Statistics and Probability Letters,
\bold{123}, 177--182.
}
\examples{
## Two-sided one-sample KS boundaries for a continuous null distribution
n <- 10
q <- 0.1
A <- pmax(0, (1:n) / n - q)
B <- pmin(1, ((1:n) - 1) / n + q)

ks_c_cdf(n, A, B)
}
\keyword{distribution}
