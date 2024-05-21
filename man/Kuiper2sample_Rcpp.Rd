\name{Kuiper2sample_Rcpp}
\alias{Kuiper2sample_Rcpp}
\title{
R function calling the C++ routines that compute the p-value for a (unweighted) two-sample Kuiper test, given arbitrary data samples on the real line or on the circle with possibly repeated observations (i.e. ties)
}
\description{
Function calling directly the C++ routines that compute the exact p-value  \eqn{P(V_{m,n} \ge q)}{} for the two-sample Kuiper test, at a fixed \eqn{q}{}, \eqn{q\in [0,2]}{}, given the sample sizes \code{m}, \code{n} and the vector \code{M} containing the number of times each distinct observation is repeated in the pooled sample.
}
\usage{Kuiper2sample_Rcpp(m, n, M, q)}
\arguments{
  \item{m}{
  the sample size of first tested sample.
}
  \item{n}{
  the sample size of second tested sample. 
}
  \item{M}{
an integer-valued vector with \eqn{k}{} cells, where \eqn{k}{} denotes the number of distinct values in the ordered pooled sample of tested pair of samples(i.e. \eqn{a_1<a_2<\ldots<a_k}{}). \code{M[i]} is the number of times that \eqn{a_i}{} is repeated in the pooled sample. A valid \code{M} must have strictly positive integer values and have the sum of all cells equals to \code{m+n}.
  }
  \item{q}{
  numeric value between 0 and 2, at which the p-value \eqn{P(V_{m,n}\ge q)}{} is computed.
  }
}
\details{
Given a pair of random samples, either on the real line or the circle, denoted by \eqn{\bm{X}_m=(X_{1},..., X_{m})}{} and \eqn{\bm{Y}_n=(Y_{1},..., Y_{n})}{}, of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from some unknown cdfs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. The task is to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{}, against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}. The two-sample Kuiper goodness-of-fit statistic that is used to test this hypothesis is defined as:
\deqn{\varsigma_{m,n} = \sup [F_{m}(t) - G_n(t)] - \inf [F_{m}(t) - G_n(t)].}{}

	
The numeric array \code{M} specifies the number of \emph{repeated observations} in the pooled sample. For a particular realization of the pooled sample \eqn{\bm{Z}_{m,n}=(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, let there be \eqn{k}{} distinct values, \eqn{a_1<a_2<...<a_k}{}, in the ordered, pooled sample \eqn{(z_1\leq z_2\leq \ldots \leq z_{m+n})}{}, where \eqn{k\leq m+n}{}, and where \eqn{m_i}{} = \code{M[i]} is the number of times \eqn{a_i}{}, \eqn{i=1,\ldots,k}{} appears in the pooled sample. The p-value is then defined as the probability
\deqn{P\left(V_{m,n}\geq q\right),}{}
where \eqn{V_{m,n}}{} is the two-sample Kuiper test statistic defined as \eqn{\varsigma_{m,n}}{}, for two samples \eqn{\bm{X}'_m}{} and  \eqn{\bm{Y}'_n}{} of sizes \eqn{m}{} and \eqn{n}{}, \emph{randomly drawn from the pooled sample without replacement}, i.e. \eqn{V_{m,n}}{} is defined on the space \eqn{\Omega}{} (see further details in \code{\link{Kuiper2sample}}), and \eqn{q}{} is the observed value of \eqn{\varsigma_{m,n}}{} with respect to the two observed samples \eqn{\bm{X}_m}{} and \eqn{\bm{Y}_n}{}.


\code{\link{Kuiper2sample_Rcpp}} implements an algorithm from Dimitrova, Jia, Kaishev (2024),  that is based on extending the algorithm provided by Nikiforov (1994) and generalizing the method due to  Maag and Stephens (1968) and Hirakawa (1973). A version of the Nikiforov's recurrence proposed recently by Viehmann (2021) is further incorporated, which computes directly the p-value, with up to 4 digits extra accuracy, but at up to 3 times higher computational cost than \code{\link{Kuiper2sample_c_Rcpp}}. It is accurate and valid for \emph{arbitrary (possibly large) sample sizes}. This algorithm ensures a total worst-case run-time of order \eqn{O((mn)^{2})}{}. When  \code{m} and \code{n} have large greatest common divisor (an extreme case is \code{m} = \code{n}), it ensures a total worst-case run-time of order \eqn{O((m)^{2}n)}{}. 

Other known implementations for the two-sample Kuiper test mainly use the approximation method or Monte Carlo simulation (See also \code{\link{Kuiper2sample}}). The former method is invalid for data with ties and often gives p-values with large errors when sample sizes are small, the latter method is usually slow and inaccurate. Compared with other known algorithms, \code{\link{Kuiper2sample_Rcpp}}  allows data samples to come from \emph{continuous, discrete or mixed distribution} (i.e. ties may appear), and is more accurate and generally applicable for \emph{large sample sizes}.

}
\value{
Numeric value corresponding to \eqn{P(V_{m,n}\geq q)}{}, given sample sizes \code{m}, \code{n} and \code{M}. If the value of \code{m}, \code{n} are non-positive, or their least common multiple exceeds the limit 2147483647, then the function returns \code{-1}, the non-permitted value of \code{M} returns \code{-2}, numerically unstable calculation returns \code{-3}.
}

\references{

Maag, U. R., Stephens, M. A. (1968). The \eqn{V_{NM}} Two-Sample Test. The Annals of Mathematical Statistics, \bold{39}(3), 923-935.

Hirakawa, K. (1973). The two-sample Kuiper test. TRU Mathematics, \bold{9}, 99-118.

Nikiforov, A. M. (1994). "Algorithm AS 288: Exact Smirnov Two-Sample Tests for Arbitrary Distributions." Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{43}(1), 265–270.

Viehmann, T. (2021). Numerically more stable computation of the p-values for the two-sample Kolmogorov-Smirnov test. \emph{arXiv preprint} arXiv:2102.08037.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests". \emph{submitted}
}
\examples{
## Computing the unweighted two-sample Kolmogorov-Smirnov test
## Example see in Nikiforov (1994)

m <- 120
n <- 150
q <- 0.183333333
M <- c(80,70,40,80)
Kuiper2sample_Rcpp(m, n, M, q)
}

