\name{Kuiper2sample}
\alias{Kuiper2sample}
\title{
Computes the p-value for a two-sample Kuiper test, given arbitrary data samples on the real line or on the circle with possibly repeated observations (i.e. ties)
}
\description{
Computes the p-value, \eqn{P(V_{m,n} \geq q)}{}, where \eqn{V_{m,n}}{} is the two-sample Kuiper test statistic, \eqn{q}{} is the observed value of the Kuiper statistic, computed based on two data samples \eqn{\{x_{1},..., x_{m}\}}{} and \eqn{\{y_{1},..., y_{n}\}}{} that may come from continuous, discrete or mixed distribution, i.e. they may have repeated observations (ties).
}
\usage{Kuiper2sample(x, y, conservative = F, tail = T)}
\arguments{
  \item{x}{
     a numeric vector of data sample values \eqn{\{x_{1}, ..., x_{m}\}}{}
}
  \item{y}{
     a numeric vector of data sample values \eqn{\{y_{1}, ..., y_{n}\}}{}
}
  \item{conservative}{
     logical variable indicating whether ties should be considered. See \sQuote{Details} for the meaning.
}  
  \item{tail}{
     logical variable indicating whether a p-value, \eqn{P(V_{m,n} \ge q)}{} or one minus the p-value, \eqn{P(V_{m,n} < q)}{}, should be computed. By default, the p-value \eqn{P(V_{m,n} \ge q)}{} is computed. See \sQuote{Details} for the meaning.
}
}

\details{
Given a pair of random samples, either on the real line or the circle, denoted by \eqn{\bm{X}_m=(X_{1},..., X_{m})}{} and \eqn{\bm{Y}_n=(Y_{1},..., Y_{n})}{}, of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from some unknown cdfs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. The task is to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{}, against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}. The two-sample Kuiper goodness-of-fit statistic that is used to test this hypothesis is defined as:
\deqn{\varsigma_{m,n} = \sup [F_{m}(t) - G_n(t)] - \inf [F_{m}(t) - G_n(t)].}{}
	
	
For a particular realization of the pooled sample \eqn{\bm{Z}_{m,n}=(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, let there be \eqn{k}{} distinct values, \eqn{a_1<a_2<...<a_k}{}, in the ordered, pooled sample \eqn{(z_1\leq z_2\leq \ldots \leq z_{m+n})}{}, where \eqn{k\leq m+n}{}, and where \eqn{m_i}{} is the number of times \eqn{a_i}{}, \eqn{i=1,\ldots,k}{} appears in the pooled sample. The p-value is then defined as the probability
\deqn{p=P\left(V_{m,n}\geq q\right),}{}
where \eqn{V_{m,n}}{} is the two-sample Kuiper test statistic defined as \eqn{\varsigma_{m,n}}{}, for two samples \eqn{\bm{X}'_m}{} and  \eqn{\bm{Y}'_n}{} of sizes \eqn{m}{} and \eqn{n}{}, \emph{randomly drawn from the pooled sample without replacement} and \eqn{q}{} is the observed value of the statistic calculated based on the user provided data samples \code{x} and \code{y}. By default \code{tail = T}, the p-value is returned, otherwise \eqn{1-p} is returned.

Note that, \eqn{V_{m,n}}{} is defined on the space \eqn{\Omega}{} of all possible pairs,  \eqn{C = \frac{(m+n)!}{m!n!}}{}  of edfs \eqn{F_m(x,\omega)}{} and \eqn{G_n(x,\omega)}{}, \eqn{\omega \in \Omega}{}, that correspond to the pairs of samples \eqn{\bm{X}'_m}{} and  \eqn{\bm{Y}'_n}{}, randomly drawn from, \eqn{\bm{Z}_{m+n}}{}, as follows.  First, \eqn{m}{} observations are drawn at random without replacement, forming the first sample \eqn{\bm{X}'_m}{}, with corresponding edf, \eqn{F_m(x,\omega)}{}. The remaining \eqn{n} observations are then assigned to the second sample \eqn{\bm{Y}'_n}{}, with corresponding edf \eqn{G_n(x,\omega)}{}. Observations are then replaced back in \eqn{\bm{Z}_{m+n}}{} and re-sampling is continued until the occurrence of all the \eqn{C}{} possible pairs of edfs \eqn{F_m(x,\omega)}{} and \eqn{G_n(x,\omega)}{},  \eqn{\omega \in \Omega}. The pairs of edf's may be coincident if there are ties in the data and each pair, \eqn{F_m(x,\omega)}{} and \eqn{G_n(x,\omega)}{} occurs with probability \eqn{1/C}.
	
	
\code{conservative} is a logical variable whether the test should be conducted conservatively. By default, \code{conservative = F}, \code{\link{Kuiper2sample}} returns the p-value that is defined through the conditional probability above. However, when the user has a priori knowledge that both samples are from a continuous distribution even if ties are present, for example, repeated observations are caused by rounding errors, the value \code{conservative = T} should be assigned, since the conditional probability is no longer relevant. In this case, \code{\link{Kuiper2sample}} computes p-values for the Kuiper test assuming no ties are present, and returns a p-value which is an upper bound of the true p-value. Note that, if the null hypothesis is rejected using the calculated upper bound for the p-value, it should also be rejected with the true p-value.
	
	
	\code{\link{Kuiper2sample}} calculates the exact p-value of the Kuiper test using an algorithm from Dimitrova, Jia, Kaishev (2024), which is based on extending the algorithm provided by Nikiforov (1994) and generalizing the method due to  Maag and Stephens (1968) and Hirakawa (1973). If \code{tail = F}, \code{\link{Kuiper2sample}} calculates the complementary p-value \eqn{1-p}. For the purpose, an exact algorithm which generalizes the method due to Nikiforov (1994) is implemented. Alternatively, if \code{tail = T}, a version of the Nikiforov's recurrence proposed recently by Viehmann (2021) is further incorporated, which computes directly the p-value, with up to 4 digits extra accuracy, but at up to 3 times higher computational cost. It is accurate and valid for \emph{arbitrary (possibly large) sample sizes}. This algorithm ensures a total worst-case run-time of order \eqn{O((mn)^{2})}{}. When  \code{m} and \code{n} have large greatest common divisor (an extreme case is \code{m} = \code{n}), it ensures a total worst-case run-time of order \eqn{O((m)^{2}n)}{}. 

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

Viehmann, T. (2021). Numerically more stable computation of the p-values for the two-sample Kolmogorov-Smirnov test. \emph{arXiv preprint} arXiv:2102.08037.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests". \emph{submitted}
}

\examples{
##Computes discrete circular data
data1 <- c(rep(pi/2,30),rep(pi,30),rep(3*pi/2,30),rep(2*pi,30))
data2 <- c(rep(pi/2,50),rep(pi,40),rep(3*pi/2,10),rep(2*pi,50))
Kuiper2sample(data1, data2)

##The calculated p-value does not change with the choice of the original point
data3 <- c(rep(pi/2,30),rep(pi,30),rep(3*pi/2,30),rep(2*pi,30))
data4 <- c(rep(pi/2,50),rep(pi,50),rep(3*pi/2,40),rep(2*pi,10))
Kuiper2sample(data3, data4)
}

