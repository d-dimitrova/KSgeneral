#ifndef KSGENERAL_H
#define KSGENERAL_H

#include <vector>

/*
 * Public C++ interface to the numerical routines in KSgeneral.
 *
 * This header is intended for users who compile the KSgeneral source directly.
 * It deliberately contains no R or Rcpp types.  The Rcpp-facing wrappers are
 * centralized in KSgeneral.cpp.
 *
 * The implementation still depends on the same numerical sources/libraries as
 * the R package (in particular FFTW for the one-sample Exact-KS-FFT routine).
 */

namespace KSgeneral {

/*
 * Exact-KS-FFT complementary probability from pre-computed boundaries.
 *
 * n : sample size
 * A : lower boundary values A_i
 * B : upper boundary values B_i
 *
 * Returns 1 - P(A_i <= U_(i) <= B_i, i = 1,...,n).
 */
double ks_c_cdf(long n,
                const std::vector<double>& A,
                const std::vector<double>& B);

/*
 * Exact two-sample Kolmogorov-Smirnov p-value P(D_{m,n} >= q).
 *
 * kind = 1 : two-sided
 * kind = 2 or 3 : one-sided alternatives used by KSgeneral
 *
 * M contains multiplicities of the distinct observations in the pooled sample.
 * w_vec contains the m+n-1 positive weights.
 */
double KS2sample(int m,
                 int n,
                 int kind,
                 const std::vector<int>& M,
                 double q,
                 const std::vector<double>& w_vec,
                 double tol);

/*
 * Complementary two-sample Kolmogorov-Smirnov probability P(D_{m,n} < q).
 * Arguments have the same meaning as for KS2sample().
 */
double KS2sample_c(int m,
                   int n,
                   int kind,
                   const std::vector<int>& M,
                   double q,
                   const std::vector<double>& w_vec,
                   double tol);

/*
 * Exact two-sample Kuiper p-value P(V_{m,n} >= q).
 */
double Kuiper2sample(int m,
                     int n,
                     const std::vector<int>& M,
                     double q);

/*
 * Complementary two-sample Kuiper probability P(V_{m,n} < q).
 */
double Kuiper2sample_c(int m,
                       int n,
                       const std::vector<int>& M,
                       double q);

} // namespace KSgeneral

#endif // KSGENERAL_H
