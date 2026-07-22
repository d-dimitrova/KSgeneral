#ifndef KSGENERAL_H
#define KSGENERAL_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#  if defined(KSG_BUILD_DLL)
#    define KSG_API __declspec(dllexport)
#  elif defined(KSG_USE_DLL)
#    define KSG_API __declspec(dllimport)
#  else
#    define KSG_API
#  endif
#  define KSG_CALL __cdecl
#else
#  define KSG_API __attribute__((visibility("default")))
#  define KSG_CALL
#endif

#ifdef __cplusplus
#  define KSG_NOEXCEPT noexcept
extern "C" {
#else
#  define KSG_NOEXCEPT
#endif

typedef int32_t ksg_status_t;
typedef int32_t ksg_alternative_t;
typedef int32_t ksg_probability_t;

enum { KSG_OK = 0, KSG_EINVAL = 1, KSG_ERANGE = 2, KSG_ENUMERIC = 3, KSG_ENOMEM = 4, KSG_EINTERNAL = 5 };
enum { KSG_ALT_TWO_SIDED = 1, KSG_ALT_GREATER = 2, KSG_ALT_LESS = 3 };
enum { KSG_PROB_LT = 1, KSG_PROB_GE = 2 };

KSG_API uint32_t KSG_CALL ksg_abi_version(void);
KSG_API const char* KSG_CALL ksg_version_string(void);
KSG_API size_t KSG_CALL ksg_last_error(char* buffer, size_t capacity);

KSG_API ksg_status_t KSG_CALL ksg_ks2_probability_summary(
    int32_t m, int32_t n, ksg_alternative_t alternative,
    const int32_t* multiplicities, size_t multiplicity_count,
    double statistic, const double* weights, size_t weight_count,
    double tolerance, ksg_probability_t probability_type,
    double* out_probability) KSG_NOEXCEPT;

KSG_API ksg_status_t KSG_CALL ksg_kuiper2_probability_summary(
    int32_t m, int32_t n, const int32_t* multiplicities,
    size_t multiplicity_count, double statistic,
    ksg_probability_t probability_type, double* out_probability) KSG_NOEXCEPT;

KSG_API ksg_status_t KSG_CALL ksg_one_sample_boundary_probability(
    int32_t n, const double* g_steps, size_t g_count,
    const double* h_steps, size_t h_count, int32_t use_fft,
    double* out_probability) KSG_NOEXCEPT;

KSG_API ksg_status_t KSG_CALL ksg_continuous_ks_ccdf(
    int32_t n, double q, double* out_probability) KSG_NOEXCEPT;

#ifdef __cplusplus
}
#endif

#endif
