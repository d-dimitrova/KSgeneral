#include "k1sample.h"
#include "KSgeneral.h"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "read_boundaries_file.h"
#include "two_sided_noncrossing_probability.h"


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
 * file interface.
 *
 * Boundary validation and crossing detection are deliberately combined in a
 * single O(n) pass. This avoids the multiple full-vector scans and temporary
 * logical vectors that equivalent R-level checks would require.
 */
KSCdfResult ks_cdf_impl(
    long n,
    const std::vector<double>& B_steps,
    const std::vector<double>& A_steps)
{
    if (n <= 0) {
        throw std::invalid_argument("'n' must be a positive integer");
    }

    if (A_steps.size() != static_cast<std::size_t>(n) ||
        B_steps.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument(
            "'A' and 'B' must both have length 'n'");
    }

    bool has_nonfinite = false;
    bool out_of_range = false;
    bool not_nondecreasing = false;
    long crossing_index = -1;

    for (long i = 0; i < n; ++i) {
        const double a = A_steps[static_cast<std::size_t>(i)];
        const double b = B_steps[static_cast<std::size_t>(i)];

        if (!std::isfinite(a) || !std::isfinite(b)) {
            has_nonfinite = true;
        }

        if (a < 0.0 || a > 1.0 || b < 0.0 || b > 1.0) {
            out_of_range = true;
        }

        if (i > 0) {
            const std::size_t prev = static_cast<std::size_t>(i - 1);
            if (a < A_steps[prev] || b < B_steps[prev]) {
                not_nondecreasing = true;
            }
        }

        if (crossing_index < 0 && b < a) {
            crossing_index = i;
        }
    }

    /*
     * Keep the same validation priority as the previous R wrapper:
     * finite/non-missing -> range -> monotonicity -> crossing warning.
     */
    if (has_nonfinite) {
        throw std::invalid_argument(
            "'A' and 'B' must contain only finite, non-missing values");
    }

    if (out_of_range) {
        throw std::invalid_argument(
            "'A' and 'B' must contain values in [0, 1]");
    }

    if (not_nondecreasing) {
        throw std::invalid_argument(
            "'A' and 'B' must be nondecreasing");
    }

    KSCdfResult result;
    result.crossing_index = crossing_index;

    if (crossing_index >= 0) {
        /*
         * The empirical process cannot remain between crossed boundaries.
         * Thus the non-crossing probability is zero and the complementary
         * probability returned by ks_c_cdf() is one.
         */
        result.probability = 1.0;
        return result;
    }

    const bool use_fft = true;
    const double noncrossing_probability =
        ecdf_noncrossing_probability_prechecked(
            n, B_steps, A_steps, use_fft);

    result.probability = 1.0 - noncrossing_probability;
    return result;
}


namespace KSgeneral {

double ks_c_cdf(long n,
                const std::vector<double>& A,
                const std::vector<double>& B)
{
    /* Preserve the historical internal boundary ordering in one place. */
    return ks_cdf_impl(n, B, A).probability;
}

} // namespace KSgeneral
