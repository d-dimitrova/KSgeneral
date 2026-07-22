#include "ksgeneral/ksgeneral.h"

#include "ksg_core.h"
#include "two_sided_noncrossing_probability.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <new>
#include <numeric>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;

ksg_status_t set_error(ksg_status_t status, const char* message) noexcept {
    last_error = message ? message : "native calculation failed";
    return status;
}

void clear_error() noexcept { last_error.clear(); }

ksg_status_t convert_multiplicities(const int32_t* input, size_t count,
                                    std::vector<int>& output) {
    if (count == 0) return set_error(KSG_ERANGE, "multiplicity_count must be positive");
    if (input == nullptr) return set_error(KSG_EINVAL, "multiplicities is null");
    if (count > static_cast<size_t>(INT_MAX)) return set_error(KSG_ERANGE, "multiplicity_count exceeds INT_MAX");
    output.clear();
    output.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        if (input[i] <= 0) return set_error(KSG_ERANGE, "multiplicities must be positive");
        output.push_back(static_cast<int>(input[i]));
    }
    return static_cast<ksg_status_t>(KSG_OK);
}

ksg_status_t validate_result(double value) noexcept {
    if (!std::isfinite(value)) return set_error(KSG_ENUMERIC, "calculation produced a non-finite value");
    if (value < 0.0) return set_error(KSG_ENUMERIC, "native calculation reported a numerical or input error");
    return static_cast<ksg_status_t>(KSG_OK);
}

template <typename F>
ksg_status_t guard(double* out, F&& f) noexcept {
    if (out == nullptr) return set_error(KSG_EINVAL, "out_probability is null");
    *out = std::numeric_limits<double>::quiet_NaN();
    try { return f(); }
    catch (const std::bad_alloc&) { return set_error(KSG_ENOMEM, "memory allocation failed"); }
    catch (const std::exception& e) { return set_error(KSG_EINTERNAL, e.what()); }
    catch (...) { return set_error(KSG_EINTERNAL, "unknown native exception"); }
}
}

extern "C" {

uint32_t KSG_CALL ksg_abi_version(void) { return 1u; }

const char* KSG_CALL ksg_version_string(void) { return "KSgeneral C ABI 1"; }

size_t KSG_CALL ksg_last_error(char* buffer, size_t capacity) {
    const size_t required = last_error.size();
    if (buffer != nullptr && capacity > 0) {
        const size_t n = std::min(required, capacity - 1);
        std::memcpy(buffer, last_error.data(), n);
        buffer[n] = '\0';
    }
    return required;
}

ksg_status_t KSG_CALL ksg_ks2_probability_summary(
    int32_t m, int32_t n, ksg_alternative_t alternative,
    const int32_t* multiplicities, size_t multiplicity_count,
    double statistic, const double* weights, size_t weight_count,
    double tolerance, ksg_probability_t probability_type,
    double* out_probability) noexcept {
    return guard(out_probability, [&]() {
        clear_error();
        if (m <= 0 || n <= 0) return set_error(KSG_ERANGE, "m and n must be positive");
        if (alternative < KSG_ALT_TWO_SIDED || alternative > KSG_ALT_LESS) return set_error(KSG_EINVAL, "invalid alternative");
        if (probability_type != KSG_PROB_LT && probability_type != KSG_PROB_GE) return set_error(KSG_EINVAL, "invalid probability_type");
        if (!std::isfinite(statistic)) return set_error(KSG_ERANGE, "statistic must be finite");
        if (weights == nullptr) return set_error(KSG_EINVAL, "weights is null");
        if (weight_count != static_cast<size_t>(m + n - 1)) return set_error(KSG_ERANGE, "weight_count must equal m + n - 1");
        if (weight_count > static_cast<size_t>(INT_MAX)) return set_error(KSG_ERANGE, "weight_count exceeds INT_MAX");
        for (size_t i = 0; i < weight_count; ++i) if (!std::isfinite(weights[i]) || weights[i] <= 0.0) return set_error(KSG_ERANGE, "weights must be positive finite values");
        std::vector<int> native_m;
        ksg_status_t status = convert_multiplicities(multiplicities, multiplicity_count, native_m);
        if (status != KSG_OK) return status;
        if (std::accumulate(native_m.begin(), native_m.end(), 0) != m + n) return set_error(KSG_ERANGE, "sum(multiplicities) must equal m + n");
        const double lt = ks2sample_c_cpp(m, n, alternative, native_m.data(), static_cast<int>(native_m.size()), statistic, weights, static_cast<int>(weight_count), tolerance);
        status = validate_result(lt);
        if (status != KSG_OK) return status;
        *out_probability = (probability_type == KSG_PROB_LT) ? lt : 1.0 - lt;
        return static_cast<ksg_status_t>(KSG_OK);
    });
}

ksg_status_t KSG_CALL ksg_kuiper2_probability_summary(
    int32_t m, int32_t n, const int32_t* multiplicities,
    size_t multiplicity_count, double statistic,
    ksg_probability_t probability_type, double* out_probability) noexcept {
    return guard(out_probability, [&]() {
        clear_error();
        if (m <= 0 || n <= 0) return set_error(KSG_ERANGE, "m and n must be positive");
        if (probability_type != KSG_PROB_LT && probability_type != KSG_PROB_GE) return set_error(KSG_EINVAL, "invalid probability_type");
        std::vector<int> native_m;
        ksg_status_t status = convert_multiplicities(multiplicities, multiplicity_count, native_m);
        if (status != KSG_OK) return status;
        if (std::accumulate(native_m.begin(), native_m.end(), 0) != m + n) return set_error(KSG_ERANGE, "sum(multiplicities) must equal m + n");
        const double lt = kuiper2sample_c_cpp(m, n, native_m.data(), static_cast<int>(native_m.size()), statistic);
        status = validate_result(lt);
        if (status != KSG_OK) return status;
        *out_probability = (probability_type == KSG_PROB_LT) ? lt : 1.0 - lt;
        return static_cast<ksg_status_t>(KSG_OK);
    });
}

ksg_status_t KSG_CALL ksg_one_sample_boundary_probability(
    int32_t n, const double* g_steps, size_t g_count,
    const double* h_steps, size_t h_count, int32_t use_fft,
    double* out_probability) noexcept {
    return guard(out_probability, [&]() {
        clear_error();
        if (n <= 0) return set_error(KSG_ERANGE, "n must be positive");
        if (!g_steps || !h_steps) return set_error(KSG_EINVAL, "boundary arrays must not be null");
        std::vector<double> g(g_steps, g_steps + g_count), h(h_steps, h_steps + h_count);
        const double p = 1.0 - ecdf_noncrossing_probability(n, g, h, use_fft != 0);
        ksg_status_t status = validate_result(p);
        if (status != KSG_OK) return status;
        *out_probability = p;
        return static_cast<ksg_status_t>(KSG_OK);
    });
}

ksg_status_t KSG_CALL ksg_continuous_ks_ccdf(int32_t n, double q, double* out_probability) noexcept {
    return guard(out_probability, [&]() {
        clear_error();
        if (n <= 0) return set_error(KSG_ERANGE, "n must be positive");
        if (!std::isfinite(q) || q < 0.0 || q > 1.0) return set_error(KSG_ERANGE, "q must be in [0, 1]");
        std::vector<double> row0(static_cast<size_t>(n));
        std::vector<double> row1(static_cast<size_t>(n));
        for (int32_t i = 1; i <= n; ++i) {
            row0[static_cast<size_t>(i - 1)] = std::min(1.0, (i - 1.0) / n + q);
            row1[static_cast<size_t>(i - 1)] = std::max(0.0, i / static_cast<double>(n) - q);
        }
        const double p = 1.0 - ecdf_noncrossing_probability(n, row0, row1, true);
        ksg_status_t status = validate_result(p);
        if (status != KSG_OK) return status;
        *out_probability = p;
        return static_cast<ksg_status_t>(KSG_OK);
    });
}

}
