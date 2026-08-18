#include "KSgeneral.h"
#include "k1sample.h"

#include <vector>

/*
 * R interface
 * -----------
 * All Rcpp export attributes are intentionally centralized in this file.
 *
 * The numerical implementations live in k1sample.cpp and k2sample.cpp and
 * expose an ordinary C++ interface through KSgeneral.h.
 */

// [[Rcpp::export(name = ".ks_c_cdf_direct")]]
double ks_c_cdf_direct(double n,
                       std::vector<double> B_steps,
                       std::vector<double> A_steps)
{
    return KSgeneral::ks_c_cdf(
        static_cast<long>(n), A_steps, B_steps);
}


// [[Rcpp::export(name = ".ks_c_cdf_Rcpp_legacy")]]
double ks_c_cdf_Rcpp(double n)
{
    return cont_ks_distribution(static_cast<long>(n));
}


// [[Rcpp::export]]
double KS2sample_c_Rcpp(int m,
                        int n,
                        int kind,
                        std::vector<int> M,
                        double q,
                        std::vector<double> w_vec,
                        double tol)
{
    return KSgeneral::KS2sample_c(m, n, kind, M, q, w_vec, tol);
}


// [[Rcpp::export]]
double Kuiper2sample_Rcpp(int m,
                          int n,
                          std::vector<int> M,
                          double q)
{
    return KSgeneral::Kuiper2sample(m, n, M, q);
}


// [[Rcpp::export]]
double Kuiper2sample_c_Rcpp(int m,
                            int n,
                            std::vector<int> M,
                            double q)
{
    return KSgeneral::Kuiper2sample_c(m, n, M, q);
}


// [[Rcpp::export]]
double KS2sample_Rcpp(int m,
                      int n,
                      int kind,
                      std::vector<int> M,
                      double q,
                      std::vector<double> w_vec,
                      double tol)
{
    return KSgeneral::KS2sample(m, n, kind, M, q, w_vec, tol);
}
