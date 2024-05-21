\name{KS2sample_c_Rcpp}
\alias{KS2sample_c_Rcpp}
\title{
R function calling the C++ routines that compute the complementary p-value for a (weighted) two-sample Kolmogorov-Smirnov (KS) test, given an arbitrary positive weight function and arbitrary data samples with possibly repeated observations (i.e. ties)
}
\description{
Function calling directly the C++ routines that compute the exact complementary p-value \eqn{P(D_{m,n} < q)}{} for the (weighed) two-sample one- or two-sided Kolmogorov-Smirnov statistic, at a fixed \eqn{q}{}, \eqn{q\in [0,1]}{}, given the sample sizes \code{m} and \code{n}, the vector of weights \code{w_vec} and the vector \code{M} containing the number of times each distinct observation is repeated in the pooled sample.
}
\usage{KS2sample_c_Rcpp(m, n, kind, M, q, w_vec, tol)}
\arguments{
  \item{m}{
  the sample size of first tested sample.
}
  \item{n}{
  the sample size of second tested sample.
}
  \item{kind}{
an integer value (= 1,2 or 3) which specified the alternative hypothesis. When = 1, the test is two-sided. When = 2 or 3, the test is one-sided. See \sQuote{Details} for the meaning of the possible values. Other value is invalid.
  }
  \item{M}{
an integer-valued vector with \eqn{k}{} cells, where \eqn{k}{} denotes the number of distinct values in the ordered pooled sample of tested pair of samples(i.e. \eqn{a_1<a_2<\ldots<a_k}{}). \code{M[i]} is the number of times that \eqn{a_i}{} is repeated in the pooled sample. A valid \code{M} must have strictly positive integer values and have the sum of all cells equals to \code{m+n}.
  }
  \item{q}{
  numeric value between 0 and 1, at which the p-value \eqn{P(D_{m,n}< q)}{} is computed.
  }
  \item{w_vec}{a vector with \code{m+n-1} cells, giving weights to each observation in the pooled sample. Valid \code{w_vec} must have \code{m+n-1} cells and strictly positive value. See \sQuote{Details} for the meaning of values in each cell.
  }
  \item{tol}{the value of \eqn{\epsilon}{} for computing \eqn{P(D_{m,n} \leq{}q- \epsilon)}{}, which is equivalent to \eqn{P(D_{m,n} < q)}{}. Non-positive input (\code{tol} \eqn{\leq 0}{}) or large input (\code{tol} \eqn{>}{}\code{1e-6}) are replaced by \code{tol=1e-6}. In cases when \code{m} and \code{n} have large least common multiple, a smaller value is highly recommended.
  }
}
\details{
Given a pair of random samples \eqn{\bm{X}_m=(X_{1},..., X_{m})}{} and \eqn{\bm{Y}_n=(Y_{1},..., Y_{n})}{} of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from some unknown cdfs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. The task is to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{},  either against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}, which corresponds to the two-sided test, or against \eqn{H_1: F(x)> G(x)}{}  and \eqn{H_1: F(x)< G(x)}{}   for at least one \eqn{x}{}, which corresponds to the two one-sided tests. The (weighted) two-sample Kolmogorov-Smirnov goodness-of-fit statistics that are used to test these hypotheses are generally defined as:
\deqn{\Delta_{m,n} = \sup |F_{m}(t) - G_n(t)|W(E_{m+n}(t), \textnormal{ to test against the alternative } H_1: F(x)\neq G(x)}{} 
\deqn{\Delta_{m,n}^{+} = \sup [F_{m}(t) - G_n(x)]W(E_{m+n}(t)), \textnormal{ to test against the alternative } H_1: F(x)> G(x)}{} 
\deqn{\Delta_{m,n}^{-} = \sup [G_n(t) - F_{m}(x)]W(E_{m+n}(t)), \textnormal{ to test against the alternative } H_1: F(x)< G(x), }{}


where \eqn{E_{m+n}(t)}{} is the empirical cdf of the pooled sample \eqn{\bm{Z}_{m,n}=(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, \eqn{W( )}{} is a strictly positive weight function defined on \eqn{[0,1]}{}. 
	
	
	\code{w_vec[i]} (0\eqn{<}{}\code{i}\eqn{<}{}\eqn{m+n}{}) is then equal to \eqn{W(Z_i)=W(\frac{i}{m+n})}{}(\eqn{Z_i}{} is the i-th smallest observation in the pooled sample \eqn{\bm{Z}_{m,n}}{}). 
Different value of \code{w_vec} specifies the weighted Kolmogorov-Smirnov test differently. For example, when \code{w_vec=rep(1,m+n-1)}, \code{\link{KS2sample_Rcpp}} calculates the p-value of the unweighted two-sample Kolmogorov-Smirnov test, when \code{w_vec} = \code{((1:(m+n-1))*((m+n-1):1))^(-1/2)}, it calculates the p-value for the weighted two-sample Kolmogorov-Smirnov test with Anderson-Darling weight \eqn{W(t) = 1/[t(1-t)]^{1/2}}{}.
	
  Possible values of \code{kind} are 1,2 and 3, which specify the alternative hypothesis, i.e. specify the test statistic to be either \eqn{\Delta_{m,n}}{},  \eqn{\Delta_{m,n}^{+}}{} or \eqn{\Delta_{m,n}^{-}}{} respectively.
 
The numeric array \code{M} specifies the number of \emph{repeated observations} in the pooled sample. For a particular realization of the pooled sample \eqn{\bm{Z}_{m,n}=(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, let there be \eqn{k}{} distinct values, \eqn{a_1<a_2<...<a_k}{}, in the ordered, pooled sample \eqn{(z_1\leq z_2\leq \ldots \leq z_{m+n})}{}, where \eqn{k\leq m+n}{}, and where \eqn{m_i}{}=\code{M[i]} is the number of times \eqn{a_i}{}, \eqn{i=1,\ldots,k}{} appears in the pooled sample. The calculated complementary p-value is the conditional probability:
\deqn{P(D_{m,n}< q)}{}
where \eqn{D_{m,n}}{} is the two-sample Kolmogorov-Smirnov test statistic defined according to the value of \code{weight} and \code{alternative}, for two samples \eqn{\bm{X}'_m}{} and  \eqn{\bm{Y}'_n}{} of sizes \eqn{m}{} and \eqn{n}{}, \emph{randomly drawn from the pooled sample without replacement}, i.e. \eqn{D_{m,n}}{} is defined on the space \eqn{\Omega}{} (see further details in \code{\link{KS2sample}}), and \eqn{q}{} is the observed value of \eqn{\Delta_{m,n}}{} defined according to the value of \code{weight} and \code{alternative} with respect to the two observed samples \eqn{\bm{X}_m}{} and \eqn{\bm{Y}_n}{}.

 

	\code{\link{KS2sample_c_Rcpp}} implements an exact algorithm, extending the Fortran 77 subroutine due to Nikiforov (1994), an extended functionality by allowing more flexible choice of weight, as well as for \emph{large sample sizes}. This leads to faster computation time, as well as, relatively high accuracy for very large \code{m} and \code{n} (less accurate than \code{\link{KS2sample_Rcpp}}). Compared with other known algorithms, it allows data samples come from \emph{continuous, discrete or mixed distribution}(i.e. ties may appear), and it is more efficient and more generally applicable for \emph{large sample sizes}. This algorithm ensures a total worst-case run-time of order \eqn{O(nm)}{}.
}

\value{
Numeric value corresponding to \eqn{P(D_{m,n}< q)}{}, given sample sizes \code{m}, \code{n}, \code{M} and \code{w_vec}. If the value of \code{m}, \code{n} are non-positive, or if the length of \code{w_vec} is not equal to \code{m+n-1}, then the function returns \code{-1}, the non-permitted value of \code{M} or non-permitted value inside \code{w_vec} returns \code{-2}, numerically unstable calculation returns \code{-3}.
}

\source{
Based on the Fortran subroutine by Nikiforov (1994). See also Dimitrova, Jia, Kaishev (2024).
}
\references{
Paul L. Canner (1975). "A Simulation Study of One- and Two-Sample Kolmogorov-Smirnov Statistics with a Particular Weight Function". Journal of the American Statistical Association, \bold{70}(349), 209-211.

Nikiforov, A. M. (1994). "Algorithm AS 288: Exact Smirnov Two-Sample Tests for Arbitrary Distributions." Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{43}(1), 265–270.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests".  \emph{submitted}
}
\examples{
## Computing the unweighted two-sample Kolmogorov-Smirnov test
## Example see in Nikiforov (1994)

m <- 120
n <- 150
kind <- 1
q <- 0.1
M <- c(80,70,40,80)
w_vec <- rep(1,m+n-1)
tol <- 1e-6
KS2sample_c_Rcpp(m, n, kind, M, q, w_vec, tol)

kind <- 2
KS2sample_c_Rcpp(m, n, kind, M, q, w_vec, tol)

## Computing the weighted two-sample Kolmogorov-Smirnov test
## with Anderson-Darling weight
kind <- 3
w_vec <- ((1:(m+n-1))*((m+n-1):1))^(-1/2)
KS2sample_c_Rcpp(m, n, kind, M, q, w_vec, tol)
}
