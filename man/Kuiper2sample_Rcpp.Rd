\name{Kuiper2sample_Rcpp}
\alias{Kuiper2sample_Rcpp}
\title{
R function calling the C++ routines that compute the p-value for a (unweighted) two-sample Kuiper test, given arbitrary data samples on the real line or on the circle with possibly repeated observations (i.e. ties)
}
\description{
Function calling directly the C++ routines that compute the exact p-value  \eqn{P(V_{m,n} \ge q)}{}, at a fixed \eqn{q}{}, \eqn{q\in [0,2]}{}, given the sample sizes \code{m}, \code{n} and repeated times \code{M} for each distinct observation appeared in the pooled sample.
}
\usage{Kuiper2sample_Rcpp(m, n, M, q, tol)}
\arguments{
  \item{m}{
  the sample size of first tested sample.
}
  \item{n}{
  the sample size of second tested sample. 
}
  \item{M}{
  an \eqn{1}{}-dimensional integer-valued vector with \eqn{K}{} cells, where \eqn{K}{} denotes the number of distinct values in the ordered pooled sample of tested pair of samples(i.e. \eqn{A_1<A_2<\ldots<A_K}{}). \code{M[i]} is the number of times that \eqn{A_i}{} has repeated in the pooled sample. A valid \code{M} must have strictly positive integer values and have the sum of all cells equal to \code{m+n}.
  }
  \item{q}{
  numeric value between 0 and 2, at which the p-value \eqn{P(V_{m,n}>q)}{} is computed.
  }
  \item{tol}{the value of \eqn{\epsilon}{} for computing \eqn{P(V_{m,n}\geq q+ \epsilon)}{}, non-positive value and large value (\code{tol>1e-6}) is replaced by default \code{tol=1e-6}. In cases when \code{m} and \code{n} have large least common multiple, a smaller value is highly recommended.
  }
  
}
\details{
Given a pair of random samples, either on the real line or the circle, denoted by \eqn{(X_{1},..., X_{m})}{} and \eqn{(Y_{1},..., Y_{n})}{}, of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from unknown CDFs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. We want to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{}, against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}. The two-sample Kuiper goodness-of-fit statistic that is used to test this hypothesis is defined as:
\deqn{\varsigma_{m,n} = \sup [F_{m}(t) - G_n(t)] - \inf [F_{m}(t) - G_n(t)].}{}

	
The numeric array \code{M} specifies the number of \emph{repeated observations} in the pooled sample. Given two data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{}, let there be \eqn{K}{}, \eqn{K\leq }{}\code{m+n}, distinct observation (\eqn{K}{} is the length of the array \code{M}), \eqn{A_1<A_2<\ldots<A_K}{}, in the pooled sample \eqn{\{x_{1},..., x_{m},y_{1},..., y_{n}\}}{} and denote by \eqn{M_i}{} the number of times \eqn{A_i}{} is repeated in the pooled sample (the value of \code{M[i]} = \eqn{M_i}{}), and by \eqn{r_i}{}  and \eqn{l_i}{}  the number of times \eqn{A_i}{} is repeated in the first and second samples respectively. The calculated p-value is the conditional probability:
\deqn{P(V_{m,n}\geq q|M_1=r_1+l_1,\ldots,M_K=r_K+l_k)}{}
where \eqn{V_{m,n}}{} is the two-sample Kuiper test statistic defined as \eqn{\varsigma_{m,n}}{}, for two samples of sizes \code{m} and \code{n}, \emph{randomly drawn from the pooled sample without replacement}, as before, \code{q} is the value of two-sample Kuiper statistic computed based on the data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{}.


  Note that \eqn{V_{m,n}}{} is defined on the space \eqn{\Omega}{} of all possible \eqn{\frac{(m+n)!}{m!n!}}{}  pairs of samples,  \eqn{(X_1,\ldots,X_m)}{}  and \eqn{(Y_1,\ldots,Y_n)}{}  of sizes \eqn{m}{}  and \eqn{n}{}, with empirical cdfs \eqn{F_m(x)}{}  and \eqn{G_n(x)}{}  respectively, that are randomly drawn without replacement from the pooled sample \eqn{\{x_{1},..., x_{m},y_{1},..., y_{n}\}}{}.
	
\code{\link{Kuiper2sample_Rcpp}}  implements an algorithm that is based on extending the algorithm provided by Nikiforov (1994) and generalizing the method due to  Maag and Stephens (1968) and Hirakawa (1973). It is accurate and valid for \emph{arbitrary (possibly large) sample sizes}. This algorithm ensures a total worst-case run-time of order \eqn{O((mn)^{2})}{}. When  \code{m} and \code{n} have large greatest common divisor (an extreme case is \code{m} = \code{n}), it ensures a total worst-case run-time of order \eqn{O((m)^{2}n)}{}. 

Other known implementations for the two-sample Kuiper test mainly use the approximation method or Monte Carlo simulation (See also \code{\link{Kuiper2sample}}). The former method is invalid for data with ties and often gives p-values with large errors when sample sizes are small, the latter method is usually slow and inaccurate. Compared with other known algorithms, \code{\link{Kuiper2sample_Rcpp}}  allows data samples to come from \emph{continuous, discrete or mixed distribution} (i.e. ties may appear), and is more accurate and generally applicable for \emph{large sample sizes}.

}
\value{
Numeric value corresponding to \eqn{P(V_{m,n}\geq q|M_1=r_1+l_1,\ldots,M_K=r_K+l_k)}{}, given sample sizes \code{m}, \code{n} and \code{M}. If the value of \code{m}, \code{n} are non-positive, then the function returns \code{-1}, the non-permitted value of \code{M} returns \code{-2}, numerically unstable calculation returns \code{-3}.
}

\references{

Maag, U. R., Stephens, M. A. (1968). The \eqn{V_{NM}} Two-Sample Test. The Annals of Mathematical Statistics, \bold{39}(3), 923-935.

Hirakawa, K. (1973). The two-sample Kuiper test. TRU Mathematics, \bold{9}, 99-118.

Nikiforov, A. M. (1994). "Algorithm AS 288: Exact Smirnov Two-Sample Tests for Arbitrary Distributions." Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{43}(1), 265–270.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests". \emph{submitted}
}
\examples{
## Computing the unweighted two-sample Kolmogorov-Smirnov test
## Example see in Nikiforov (1994)

m <- 120
n <- 150
q <- 0.183333333
M <- c(80,70,40,80)
tol <- 1e-6
Kuiper2sample_Rcpp(m, n, M, q, tol)
}

