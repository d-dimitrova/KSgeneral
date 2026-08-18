#include "k1sample.h"
#include "KSgeneral.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "read_boundaries_file.h"
#include "two_sided_noncrossing_probability.h"

namespace {

void validate_one_sample_boundaries(long n,
                                    const std::vector<double>& A,
                                    const std::vector<double>& B)
{
    if (n <= 0) {
        throw std::invalid_argument("n must be a positive integer");
    }

    if (A.size() != static_cast<std::size_t>(n) ||
        B.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument("A and B must both have length n");
    }

    const auto valid_boundary = [](const std::vector<double>& x) {
        return std::all_of(x.begin(), x.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0 && value <= 1.0;
        }) && std::is_sorted(x.begin(), x.end());
    };

    if (!valid_boundary(A) || !valid_boundary(B)) {
        throw std::invalid_argument(
            "A and B must be nondecreasing finite vectors with values in [0, 1]");
    }
}

} // namespace


/*
 * Legacy file-based one-sample entry point.
 *
 * This is retained only to back the deprecated R function ks_c_cdf_Rcpp().
 * New C++ code should use KSgeneral::ks_c_cdf().
 */
double cont_ks_distribution(long n)
{
    std::pair<std::vector<double>, std::vector<double> > bounds =
        read_boundaries_file("Boundary_Crossing_Time.txt");

    const std::vector<double>& first_boundary_steps = bounds.first;
    const std::vector<double>& second_boundary_steps = bounds.second;

    const bool use_fft = true;
    return 1.0 - ecdf_noncrossing_probability(
        n, first_boundary_steps, second_boundary_steps, use_fft);
}


/*
 * Internal direct-boundary Exact-KS-FFT calculation.
 *
 * B_steps and A_steps preserve the positional order used by the historical
 * file interface. The public API below accepts the natural (A, B) order.
 */
double ks_cdf_impl(long n,
                   const std::vector<double>& B_steps,
                   const std::vector<double>& A_steps)
{
    const bool use_fft = true;
    return 1.0 - ecdf_noncrossing_probability(n, B_steps, A_steps, use_fft);
}


namespace KSgeneral {

double ks_c_cdf(long n,
                const std::vector<double>& A,
                const std::vector<double>& B)
{
    validate_one_sample_boundaries(n, A, B);

    /* Preserve the historical internal boundary ordering in one place. */
    return ks_cdf_impl(n, B, A);
}

} // namespace KSgeneral
