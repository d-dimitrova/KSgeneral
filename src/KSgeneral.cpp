#include <Rcpp.h>

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
Rcpp::List ks_c_cdf_direct(double n,
                           std::vector<double> B_steps,
                           std::vector<double> A_steps)
{
    /*
     * R/ks_c_cdf.R performs only the inexpensive scalar/type/length checks.
     * ks_cdf_impl() performs elementwise validation and crossing detection in
     * one C++ pass before entering the numerical calculation.
     */
    const KSCdfResult result = ks_cdf_impl(
        static_cast<long>(n), B_steps, A_steps);

    return Rcpp::List::create(
        Rcpp::_["value"] = result.probability,
        Rcpp::_["crossing"] = result.crossing_index);
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
