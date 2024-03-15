\name{KS2sample}
\alias{KS2sample}
\title{
Computes the p-value for a (weighted) two-sample Kolmogorov-Smirnov test, given an arbitrary positive weight function and arbitrary data samples with possibly repeated observations (i.e. ties)
}
\description{
Computes the p-value \eqn{P(D_{m,n} \ge q)}{}, where \eqn{D_{m,n}}{} is the one- or two-sided two-sample Kolmogorov-Smirnov test statistic with weight function \code{weight}, \eqn{q}{} is the observed value of KS statistic computed based on two data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{} that may come from continuous, discrete or mixed distribution, i.e. they may have repeated observations (ties).
}
\usage{
KS2sample(x, y, alternative = c("two.sided", "less", "greater"),
conservative = F, weight = 0, tol = 1e-08)
}
\arguments{
  \item{x}{
     a numeric vector of data sample values \eqn{\{x_{1}, ..., x_{m}\}}{}.
}
  \item{y}{
     a numeric vector of data sample values \eqn{\{y_{1}, ..., y_{n}\}}{}
}
  \item{alternative}{
     Indicates the alternative hypothesis and must be one of "two.sided" (default), "less", or "greater". You can specify just the initial letter of the value, but the argument name must be given in full.  See \sQuote{Details} for the meaning of the possible values.
}
  \item{conservative}{
     logical indicating whether ties should be considered. See \sQuote{Details} for the meaning.
}
  \item{weight}{
     either a numerical value between 0 and 1 which specifies form weight function from a class of pre-defined functions, or a self-defined strictly positive function containing one variable. By default, no weight function is assumed. See \sQuote{Details} for the meaning of the possible values.
}
  \item{tol}{ 
     the value of \eqn{\epsilon}{} for computing \eqn{P(D_{m,n} \geq{}q+ \epsilon)}{}. Non-positive value and large value (\code{tol>1e-6}) are replaced by the default \code{tol=1e-6}. In cases when \code{m} and \code{n} have large least common multiple, a smaller value is highly recommended.
  }
}
\details{
Given a pair of random samples \eqn{(X_{1},..., X_{m})}{} and \eqn{(Y_{1},..., Y_{n})}{} of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from unknown CDFs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. We want to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{},  either against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}, which corresponds to the two-sided test, or against \eqn{H_1: F(x)> G(x)}{}  and \eqn{H_1: F(x)< G(x)}{}   for at least one \eqn{x}{}, which corresponds to the two one-sided tests. The (weighted) two-sample Kolmogorov-Smirnov goodness-of-fit statistics that are used to test these hypotheses are generally defined as:
\deqn{\Delta_{m,n} = \sup |F_{m}(t) - G_n(t)|W(E_{m+n}(t), \textnormal{ to test against the alternative } H_1: F(x)\neq G(x)}{} 
\deqn{\Delta_{m,n}^{+} = \sup [F_{m}(t) - G_n(x)]W(E_{m+n}(t)), \textnormal{ to test against the alternative } H_1: F(x)> G(x)}{} 
\deqn{\Delta_{m,n}^{-} = \sup [G_n(t) - F_{m}(x)]W(E_{m+n}(t)), \textnormal{ to test against the alternative } H_1: F(x)< G(x), }{}

	where \eqn{E_{m+n}(t)}{} is the empirical cdf of the pooled sample \eqn{(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, \eqn{W( )}{} is a strictly positive weight function defined on \eqn{[0,1]}{}. 
	
	
	Possible values of \code{alternative} are \code{"two.sided"}, \code{"greater"} and \code{"less"} which specify the alternative hypothesis, i.e. specify the test statistics to be either \eqn{\Delta_{m,n}}{}, \eqn{\Delta_{m,n}^{+}}{} or  \eqn{\Delta_{m,n}^{-}}{} respectively. 

	
  When \code{weight} is assigned with a numeric value \eqn{\nu}{} between \code{0} and \code{1}, the test statistic is specified as the weighted two-sample Kolmorogov-Smirnov test with generalized Anderson-Darling weight \eqn{W(t)=1/[t(1-t)]^{\nu}}{} (see Finner and Gontscharuk 2018). Then for example, the two-sided two-sample Kolmogorov-Smirnov statistic has the following form:
\deqn{\Delta_{m,n}=\sup\limits_{t} \frac{|F_m(t)-G_n(t)|}{[E_{m+n}(t)(1-E_{m+n}(t))]^{\nu}}}{}
The latter specification defines a family of weighted Kolmogorov-Smirnov tests, covering the unweighted test (when \code{weight = 0}), and the widely-known weighted Kolmogorov-Smirnov test with Anderson-Darling weight (when \code{weight = 0.5}, see definition of this statistic also in Canner 1975).
If one wants to implement a weighted test with a user-specified weight function, for example, \eqn{W(t)=1/[t(2-t)]^{1/2}}{} suggested by Buning (2001), which ensures higher power when both \code{x} and \code{y} coming from distributions that are left-skewed and heavy-tailed, one can directly assign a univariate function with output value \code{1/sqrt(t*(2-t))} to \code{weight}.  See `Examples' for this demonstration.
	
  Given two data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{}, let there be \eqn{K}{} (\eqn{K\leq }{}\code{m+n}) distinct observation, \eqn{A_1<A_2<\ldots<A_K}{}, in the pooled sample \eqn{\{x_{1},..., x_{m},y_{1},..., y_{n}\}}{} and denote by \eqn{M_i}{} the number of times \eqn{A_i}{} is repeated in the pooled sample, and by \eqn{r_i}{}  and \eqn{l_i}{}  the number of times \eqn{A_i}{} is repeated in the first and second samples respectively, the p-value is defined as the conditional probability:
\deqn{P(D_{m,n}\geq q|M_1=r_1+l_1,\ldots,M_K=r_K+l_k),}{}
where \eqn{D_{m,n}}{} is the two-sample Kolmogorov-Smirnov test statistic defined according to the value of \code{weight} and \code{alternative}, for two samples of sizes \code{m} and \code{n}, \emph{randomly drawn from the pooled sample without replacement}, as before, \code{q} is its observed value of statistic calculated based on the samples \code{x} and \code{y}.

  Note that \eqn{D_{m,n}}{} is defined on the space \eqn{\Omega}{} of all possible \eqn{\frac{(m+n)!}{m!n!}}{}  pairs of samples,  \eqn{(X_1,\ldots,X_m)}{}  and \eqn{(Y_1,\ldots,Y_n)}{}  of sizes \eqn{m}{}  and \eqn{n}{}, with empirical cdfs \eqn{F_m(x)}{}  and \eqn{G_n(x)}{}  respectively, that are randomly drawn without replacement from the pooled sample \eqn{\{x_{1},..., x_{m},y_{1},..., y_{n}\}}{}.
	
	\code{conservative} is a logical variable whether the test should be conducted conservatively. By default, \code{conservative = FALSE}, \code{\link{KS2sample}} returns the p-value that is defined through the conditional probability above. However, when the user has a priori knowledge that both samples are from a continuous distribution even if ties are present, for example, repeated observations are caused by rounding errors, the value \code{conservative = TRUE} should be assigned, since the conditional probability is no longer relevant. In this case, \code{\link{KS2sample}} computes p-values for the Kolmogorov-Smirnov test assuming no tie is present, and returns a p-value which is an upper bound of the true p-value. If the calculated upper bound p-value rejects the null hypothesis, the true p-value should also reject the null hypothesis. 
	
	\code{\link{KS2sample}} implements an exact algorithm which generalizes the method due to Nikiforov (1994), and calculates the exact p-value of the KS test. \code{\link{KS2sample}} ensures a total worst-case run-time of order \eqn{O(nm)}{}. In comparison with other known algorithms, it not only allows the flexible choice of weights which in some cases improve the statistical power (see Dimitrova, Jia, Kaishev 2024), but also is more efficient and generally applicable for \emph{large sample sizes}.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{data.name}{a character string giving names of the data.}
}
\source{
	Based on the Fortran subroutine by Nikiforov (1994). See also Dimitrova, Jia, Kaishev (2024).
}
\references{
Buning H (2001). "Kolmogorov-Smirnov- and Cramer-von Mises Type Two-sample Tests With Various Weight Functions." Communications in Statistics - Simulation and Computation, \bold{30}(4), 847-865.

Finner H, Gontscharuk V (2018). "Two-sample Kolmogorov-Smirnov-type tests revisited: Old and new tests in terms of local levels." The Annals of Statistics, \bold{46}(6A), 3014-3037.

Paul L. Canner (1975). "A Simulation Study of One- and Two-Sample Kolmogorov-Smirnov Statistics with a Particular Weight Function". Journal of the American Statistical Association, \bold{70}(349), 209-211.

Nikiforov, A. M. (1994). "Algorithm AS 288: Exact Smirnov Two-Sample Tests for Arbitrary Distributions." Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{43}(1), 265-270.

   Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests". \emph{submitted}
}
\examples{
##Computes p-value of two-sided unweighted test for continuous data 
data1 <- rexp(750, 1)
data2 <- rexp(800, 1)
KS2sample(data1, data2)
##Computes p-value of one-sided test with Anderson-Darling weight function
KS2sample(data1, data2, alternative = "greater", weight = 0.5)

##Computes p-values of two-sided test with Buning's weight function for discrete data
data3 <- rnbinom(100, size = 3, prob = 0.6)
data4 <- rpois(120, lambda = 2)
f <- function(t) 1 / sqrt( t * (2 - t) ) 
KS2sample(data3, data4, weight = f)
}
