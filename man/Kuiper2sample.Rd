\name{Kuiper2sample}
\alias{Kuiper2sample}
\title{
Computes the p-value for a two-sample Kuiper test, given arbitrary data samples on the real line or on the circle with possibly repeated observations (i.e. ties)
}
\description{
Computes the p-value \eqn{P(V_{m,n} \ge q)}{}, where \eqn{V_{m,n}}{} is the two-sample Kuiper test statistic, \eqn{q}{} is the observed value of Kuiper statistic computed based on two data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{} that may come from continuous, discrete or mixed distribution, i.e. they may have repeated observations (ties).
}
\usage{Kuiper2sample(x, y, conservative = F, tol = 1e-08)}
\arguments{
  \item{x}{
     a numeric vector of data sample values \eqn{\{x_{1}, ..., x_{m}\}}{}
}
  \item{y}{
     a numeric vector of data sample values \eqn{\{y_{1}, ..., y_{n}\}}{}
}
  \item{conservative}{
     logical indicating whether ties should be considered. See \sQuote{Details} for the meaning.
}
  \item{tol}{the value of \eqn{\epsilon}{} for computing \eqn{P(V_{m,n}\geq q+ \epsilon)}{}, non-positive value and large value (\code{tol>1e-6}) is replaced by default \code{tol=1e-6}. In cases when \code{m} and \code{n} have large least common multiple, a smaller value is highly recommended.
  }
}
\details{
Given a pair of random samples, either on the real line or the circle, denoted by \eqn{(X_{1},..., X_{m})}{} and \eqn{(Y_{1},..., Y_{n})}{}, of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from unknown CDFs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. We want to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{}, against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}. The two-sample Kuiper goodness-of-fit statistic that is used to test this hypothesis is defined as:
\deqn{\varsigma_{m,n} = \sup [F_{m}(t) - G_n(t)] - \inf [F_{m}(t) - G_n(t)].}{}
	
	
Given two data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{}, let there be \eqn{K}{} (\eqn{K\leq }{}\code{m+n}) distinct observation, \eqn{A_1<A_2<\ldots<A_K}{}, in the pooled sample \eqn{\{x_{1},..., x_{m},y_{1},..., y_{n}\}}{} and denote by \eqn{M_i}{} the number of times \eqn{A_i}{} is repeated in the pooled sample, and by \eqn{r_i}{}  and \eqn{l_i}{}  the number of times \eqn{A_i}{} is repeated in the first and second samples respectively, the p-value is defined as the conditional probability:
\deqn{P(V_{m,n}\geq q|M_1=r_1+l_1,\ldots,M_K=r_K+l_k)}{}
where \eqn{V_{m,n}}{} is the two-sample Kuiper test statistic defined as \eqn{\varsigma_{m,n}}{}, for two samples of sizes \code{m} and \code{n}, \emph{randomly drawn from the pooled sample without replacement}, as before, \code{q} is the value of two-sample Kuiper statistic computed based on the  data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{}.


Note that \eqn{V_{m,n}}{} is defined on the space \eqn{\Omega}{} of all possible \eqn{\frac{(m+n)!}{m!n!}}{}  pairs of samples,  \eqn{(X_1,\ldots,X_m)}{}  and \eqn{(Y_1,\ldots,Y_n)}{}  of sizes \eqn{m}{}  and \eqn{n}{}, with empirical cdfs \eqn{F_m(x)}{}  and \eqn{G_n(x)}{}  respectively, that are randomly drawn without replacement from the pooled sample \eqn{\{x_{1},..., x_{m},y_{1},..., y_{n}\}}{}.
	
	
\code{conservative} is a logical variable whether the test should be conducted conservatively. By default, \code{conservative = FALSE}, \code{\link{Kuiper2sample}} returns the p-value that is defined through the conditional probability above. However, when the user has a priori knowledge that both samples are from a continuous distribution even if ties are present, for example, repeated observations are caused by rounding errors, the value \code{conservative = TRUE} should be assigned, since the conditional probability is no longer relevant. In this case, \code{\link{Kuiper2sample}} computes p-values for the Kolmogorov-Smirnov test assuming no tie is present, and returns a p-value which is an upper bound of the true p-value. If the calculated upper bound p-value rejects the null hypothesis, the true p-value should also reject the null hypothesis. 
	
	
	\code{\link{Kuiper2sample}} implements an algorithm that is based on extending the algorithm provided by Nikiforov (1994) and generalizing the method due to  Maag and Stephens (1968) and Hirakawa (1973). It is accurate and valid for \emph{arbitrary (possibly large) sample sizes}. This algorithm ensures a total worst-case run-time of order \eqn{O((mn)^{2})}{}. When  \code{m} and \code{n} have large greatest common divisor (an extreme case is \code{m} = \code{n}), it ensures a total worst-case run-time of order \eqn{O((m)^{2}n)}{}. 

	\code{\link{Kuiper2sample}} is accurate and fast compared with the function based on the Monte Carlo simulation. Compared to the implementation using asymptotic method, \code{\link{Kuiper2sample}} allows data samples to come from \emph{continuous, discrete or mixed distribution} (i.e. ties may appear), and is more accurate than asymptotic method when sample sizes are small.

}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{data.name}{a character string giving names of the data.}
}

\references{

Maag, U. R., Stephens, M. A. (1968). The \eqn{V_{NM}} Two-Sample Test. The Annals of Mathematical Statistics, \bold{39}(3), 923-935.

Hirakawa, K. (1973). The two-sample Kuiper test. TRU Mathematics, \bold{9}, 99-118.

Nikiforov, A. M. (1994). "Algorithm AS 288: Exact Smirnov Two-Sample Tests for Arbitrary Distributions." Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{43}(1), 265–270.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests". \emph{submitted}
}

\examples{
##Computes discrete circular data
data1 <- c(rep(pi/2,30),rep(pi,30),rep(3*pi/2,30),rep(2*pi,30))
data2 <- c(rep(pi/2,50),rep(pi,40),rep(3*pi/2,10),rep(2*pi,50))
Kuiper2sample(data1, data2)

##The calculated p-value does not change with choice of original point
data3 <- c(rep(pi/2,30),rep(pi,30),rep(3*pi/2,30),rep(2*pi,30))
data4 <- c(rep(pi/2,50),rep(pi,50),rep(3*pi/2,40),rep(2*pi,10))
Kuiper2sample(data3, data4)
}

