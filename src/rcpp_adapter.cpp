#include <Rcpp.h>

#include "ksg_core.h"

// [[Rcpp::export]]
double ks_c_cdf_Rcpp(double n) {
    return cont_ks_distribution(static_cast<long>(n));
}

// [[Rcpp::export]]
double KS2sample_c_Rcpp(int m, int n, int kind, Rcpp::IntegerVector M,
                        double q, Rcpp::NumericVector w_vec, double tol) {
    return ks2sample_c_cpp(m, n, kind, M.begin(), M.size(), q,
                           w_vec.begin(), w_vec.size(), tol);
}

// [[Rcpp::export]]
double Kuiper2sample_Rcpp(int m, int n, Rcpp::IntegerVector M, double q) {
    return kuiper2sample_cpp(m, n, M.begin(), M.size(), q);
}

// [[Rcpp::export]]
double Kuiper2sample_c_Rcpp(int m, int n, Rcpp::IntegerVector M, double q) {
    return kuiper2sample_c_cpp(m, n, M.begin(), M.size(), q);
}

// [[Rcpp::export]]
double KS2sample_Rcpp(int m, int n, int kind, Rcpp::IntegerVector M,
                      double q, Rcpp::NumericVector w_vec, double tol) {
    return ks2sample_cpp(m, n, kind, M.begin(), M.size(), q,
                         w_vec.begin(), w_vec.size(), tol);
}
