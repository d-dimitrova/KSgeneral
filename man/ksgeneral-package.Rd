\name{KSgeneral-package}
\alias{KSgeneral-package}
\docType{package}
\title{
\packageTitle{KSgeneral}
}
\description{
This package computes p-values of the one-sample and two-sample Kolmogorov-Smirnov (KS) tests and the two-sample Kuiper test.


The one-sample two-sided Kolmogorov-Smirnov (KS) statistic is one of the most popular goodness-of-fit test statistics that is used to measure how well the distribution of a random sample agrees with a prespecified theoretical distribution.
Given a random sample \eqn{\{X_{1},..., X_{n}\}}{} of size \eqn{n}{} with an empirical cdf \eqn{F_{n}(x)}{}, the two-sided KS statistic is defined as
\eqn{D_{n} = \sup | F_{n}(x) - F(x) | }{}, where \eqn{F(x)}{}  is the cdf of the prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}{}, that \eqn{ \{ X_{1},..., X_{n} \} }{} comes from \eqn{F(x)}{}.
The package \pkg{KSgeneral} implements a novel, accurate and efficient Fast Fourier Transform (FFT)-based method, referred as Exact-KS-FFT method to compute the complementary cdf,
\eqn{P(D_{n} \ge q)}{}, at a fixed \eqn{q\in [0, 1]}{} for a given (hypothezied) purely discrete, mixed or continuous underlying cdf \eqn{F(x)}{}, and arbitrary, possibly very large sample size \eqn{n}{}.
A plot of the complementary cdf \eqn{P(D_{n} \ge q)}{}, \eqn{0 \le q \le 1}{}, can also be produced.

In other words, the package computes the p-value, \eqn{P(D_{n} \ge q)}{} for any fixed critical level \eqn{q\in [0, 1]}{}.
If an observed (data) sample, \eqn{\{x_{1},..., x_{n}\}}{}  is supplied, \pkg{KSgeneral} computes the p-value \eqn{P(D_{n} \ge d_{n})}{}, where \eqn{d_{n}}{}  is the value of the KS test statistic computed based on \eqn{\{x_{1},..., x_{n}\}}{}.  One can also compute the (complementary) cdf for the one-sided KS statistics \eqn{D_{n}^{-}}{}  or \eqn{D_{n}^{+}}{}  (cf., Dimitrova, Kaishev, Tan (2020)) by appropriately specifying correspondingly \eqn{A_{i} = 0}{}  for all \eqn{i}{}  or \eqn{B_{i} = 1}{}  for all \eqn{i}{}, in the function \code{\link{ks_c_cdf_Rcpp}}.

The two-sample Kolmogorov-Smirnov (KS) and the Kuiper statistics are widely used to test the null hypothesis (\eqn{H_0}{}) that two data samples come from the same underlying distribution. Given a pair of random samples \eqn{\bm{X}_m=(X_{1},..., X_{m})}{} and \eqn{\bm{Y}_n=(Y_{1},..., Y_{n})}{} of sizes \code{m} and \code{n} with empirical cdfs \eqn{F_{m}(t)}{} and \eqn{G_{n}(t)}{} respectively, coming from unknown CDFs \eqn{F(x)}{} and \eqn{G(x)}{}. It is assumed that \eqn{F(x)}{} and \eqn{G(x)}{} could be either \emph{continuous}, \emph{discrete} or \emph{mixed}, which means that repeated observations are allowed in the corresponding observed samples. We want to test the null hypothesis \eqn{H_0: F(x) = G(x)}{}  for all \eqn{x}{},  either against the alternative hypothesis \eqn{H_1: F(x)\neq G(x)}{}  for at least one \eqn{x}{}, which corresponds to the two-sided test, or against \eqn{H_1: F(x)> G(x)}{}  and \eqn{H_1: F(x)< G(x)}{}   for at least one \eqn{x}{}, which corresponds to the two one-sided tests. The (weighted) two-sample Kolmogorov-Smirnov goodness-of-fit statistics that are used to test these hypotheses are generally defined as:
\deqn{\Delta_{m,n} = \sup |F_{m}(t) - G_n(t)|W(E_{m+n}(t), \textnormal{ to test against the alternative } H_1: F(x)\neq G(x)}{} 
\deqn{\Delta_{m,n}^{+} = \sup [F_{m}(t) - G_n(x)]W(E_{m+n}(t)), \textnormal{ to test against the alternative } H_1: F(x)> G(x)}{} 
\deqn{\Delta_{m,n}^{-} = \sup [G_n(t) - F_{m}(x)]W(E_{m+n}(t)), \textnormal{ to test against the alternative } H_1: F(x)< G(x)}{} 
where \eqn{E_{m+n}(t)}{}  is the empirical cdf of the pooled sample \eqn{\bm{Z}_{m,n}=(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, \eqn{W( )}{}  is a strictly positive weight function defined on \eqn{(0,1)}{}. \pkg{KSgeneral} implements an exact algorithm which is an extension of the Fortran 77 subroutine due to Nikiforov (1994), to calculate the exact p-value \eqn{P(D_{m,n} \ge q)}{}, where \eqn{q}{} is the observed value of \eqn{\Delta_{m,n}}{} with respect to the two observed samples \eqn{\bm{X}_m=\{x_1,\ldots,x_m\}}{} and \eqn{\bm{Y}_n=\{y_1,\ldots,y_n\}}{} and \eqn{D_{m,n}}{}  is the two-sample Kolmogorov-Smirnov goodness-of-fit test defined on the space \eqn{\Omega}{} of all possible \eqn{\frac{(m+n)!}{m!n!}}{}  pairs of samples, \eqn{\bm{X}'_m}{} and  \eqn{\bm{Y}'_n}{} of sizes \eqn{m}{} and \eqn{n}{}, that are \emph{randomly drawn from the pooled sample} \eqn{\bm{Z}_{m+n}}{} \emph{without replacement}. Samples may come from any continuous, discrete or mixed distribution, i.e. the test allows repeated observations to appear in the user provided data samples \eqn{\{x_1,\ldots,x_m\}}{}, \eqn{\{y_1,\ldots,y_n\}}{}  and their pooled sample \eqn{\bm{Z}_{m+n}=\{x_1,\ldots,x_m,y_1,\ldots,y_n\}}{}.

The two-sample (unweighted) Kuiper goodness-of-fit statistic is defined as:
\deqn{\varsigma_{m,n} = \sup [F_{m}(t) - G_n(t)] - \inf [F_{m}(t) - G_n(t)].}{} 
It is widely used when the data samples are periodic or circular (data that are measured in radians). \pkg{KSgeneral} calculates the exact p-value \eqn{P(V_{m,n} \ge q)}{}, where \eqn{V_{m,n}}{}  is the two-sample Kuiper goodness-of-fit test defined on the on the space, \eqn{\Omega}{}, as described above, and \eqn{q}{}  is the observed value of \eqn{\varsigma_{m,n}}{} with respect to the two observed samples \eqn{\bm{X}_m=\{x_1,\ldots,x_m\}}{} and \eqn{\bm{Y}_n=\{y_1,\ldots,y_n\}}{}. Similarly, as for the KS test, the two-sample Kuiper test also allows repeated observations in the user provided data samples \eqn{\{x_1,\ldots,x_m\}}{}, \eqn{\{y_1,\ldots,y_n\}}{}  and their pooled sample \eqn{\bm{Z}_{m+n}=\{x_1,\ldots,x_m,y_1,\ldots,y_n\}}{}.
}
\details{

\bold{One-sample KS test}:


The Exact-KS-FFT method to compute p-values of the one-sample KS test in \pkg{KSgeneral} is based on expressing the p-value \eqn{P(D_{n} \ge q)} in terms of an appropriate rectangle probability with respect to the uniform order statistics, as noted by Gleser (1985) for \eqn{P(D_{n} > q)}.
The latter representation is used to express \eqn{P(D_{n} \ge q)} via a double-boundary non-crossing probability for a homogeneous Poisson process, with intensity \eqn{n}, which is then efficiently computed using FFT, ensuring total run-time of order \eqn{O(n^{2}log(n))} (see Dimitrova, Kaishev, Tan (2020) and also Moscovich and Nadler (2017) for the special case when \eqn{F(x)} is continuous).

The code for the one-sample KS test in \pkg{KSgeneral} represents an R wrapper of the original C++ code due to Dimitrova, Kaishev, Tan (2020) and based on the C++ code developed by Moscovich and Nadler (2017).
The package includes the functions \code{\link{disc_ks_c_cdf}}, \code{\link{mixed_ks_c_cdf}} and \code{\link{cont_ks_c_cdf}} that compute the complementary cdf \eqn{P(D_n \ge q)}, for a fixed \eqn{q}, \eqn{0 \le q \le 1}, when \eqn{F(x)} is purely discrete, mixed or continuous, respectively.
\pkg{KSgeneral} includes also the functions \code{\link{disc_ks_test}}, \code{\link{mixed_ks_test}} and \code{\link{cont_ks_test}} that compute the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a user provided data sample \eqn{\{x_{1}, ..., x_{n}\}}, when \eqn{F(x)} is purely discrete, mixed or continuous, respectively.


The functions \code{\link{disc_ks_test}} and \code{\link{cont_ks_test}} represent accurate and fast (run time \eqn{O(n^{2}log(n))}) alternatives to the functions \code{\link{ks.test}} from the package \pkg{dgof} and the function \code{\link{ks.test}} from the package \pkg{stat}, which compute p-values of \eqn{P(D_{n} \ge d_{n})}, assuming \eqn{F(x)} is purely discrete or continuous, respectively.

The package also includes the function \code{\link{ks_c_cdf_Rcpp}} which gives the flexibility to compute the complementary cdf (p-value) for the one-sided KS test statistics \eqn{D_{n}^{-}} or \eqn{D_{n}^{+}}.
It also allows for faster computation time and possibly higher accuracy in computing \eqn{P(D_{n} \ge q)}.


\bold{Two-sample KS test and Kuiper test}:


The method underlying for computing p-values of the two-sample KS and Kuiper tests in \pkg{KSgeneral} is the extension of the algorithm due to Nikiforov (1994) and is based on expressing the p-value as the probability that a point sequence stays within a certain region in the two-dimensional integer-valued lattice. The algorithm for both tests uses a recursive formula to calculate the total number of point sequences within the region which is divided by the total number of elements in \eqn{\Omega}{}, i.e. \eqn{\frac{(m+n)!}{m!n!}}{}  to obtain the probability.

For a particular realization of the pooled sample \eqn{\bm{Z}_{m,n}=(X_{1},..., X_{m},Y_{1},..., Y_{n})}{}, the p-values calculated by the functions \code{KS2sample} and \code{Kuiper2sample} are the probabilities:
\deqn{P(D_{m,n}\geq q), P(V_{m,n}\geq q),}{}
where \eqn{D_{m,n}}{} and \eqn{V_{m,n}}{}  are the two-sample Kolmogorov-Smirnov and Kuiper test statistics respectively, for two samples \eqn{\bm{X}'_m}{} and  \eqn{\bm{Y}'_n}{} of sizes \eqn{m}{} and \eqn{n}{}, \emph{randomly drawn from the pooled sample without replacement}, i.e. they are defined on the space \eqn{\Omega}{}. As before, \eqn{q}{} is the observed value of the statistic with respect to the two observed samples \eqn{\bm{X}_m}{} and \eqn{\bm{Y}_n}{}.


Both \code{\link{KS2sample}} and \code{\link{Kuiper2sample}} implement algorithms which generalize the method due to Nikiforov (1994), and calculate the exact p-values of the KS test and the Kuiper test respectively. Both of them allow tested data samples to come from continuous, discrete or mixed distributions (ties are also allowed).

\code{\link{KS2sample}} ensures a total worst-case run-time of order \eqn{O(nm)}. Compared with other known algorithms, it not only allows more flexible choices on weights leading to better power (see Dimitrova, Jia, Kaishev 2024), but also is more efficient and more generally applicable for \emph{large sample sizes}. \code{\link{Kuiper2sample}} is accurate and valid for large sample sizes. It ensures a total worst-case run-time of order \eqn{O((mn)^{2})}. When \code{m} and \code{n} have large greatest common divisor (an extreme case is \code{m} = \code{n}), it ensures a total worst-case run-time of order \eqn{O((m)^{2}n)}. 

}
\author{
\packageAuthor{KSgeneral}

Maintainer: \packageMaintainer{KSgeneral}
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

Gleser L.J. (1985). "Exact Power of Goodness-of-Fit Tests of Kolmogorov Type for Discontinuous Distributions". Journal of the American Statistical Association, \bold{80}(392), 954-958.

Moscovich A., Nadler B. (2017). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes". Statistics and Probability Letters, \bold{123}, 177-182.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2024). "The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests".  \emph{submitted}

}

