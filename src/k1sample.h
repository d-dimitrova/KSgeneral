#ifndef KSGENERAL_K1SAMPLE_H
#define KSGENERAL_K1SAMPLE_H

#include <vector>

/*
 * Internal one-sample Exact-KS-FFT entry points.
 *
 * ks_cdf_impl() accepts the two boundary vectors in the historical internal
 * order (B first, A second). It validates all boundary elements and records
 * the first crossing in a single pass.
 *
 * cont_ks_distribution() retains the historical file-based implementation
 * solely for the deprecated R function ks_c_cdf_Rcpp().
 */

struct KSCdfResult {
    double probability;
    long crossing_index;  // -1 when the boundaries do not cross
};

KSCdfResult ks_cdf_impl(
    long n,
    const std::vector<double>& B_steps,
    const std::vector<double>& A_steps);

double cont_ks_distribution(long n);

#endif // KSGENERAL_K1SAMPLE_H
