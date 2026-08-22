#ifndef __aligned_mem__
#define __aligned_mem__

#include <complex>
#include <fftw3.h>

inline double* allocate_aligned_doubles(int n)
{
    return static_cast<double*>(fftw_malloc(static_cast<size_t>(n) * sizeof(double)));
}

inline std::complex<double>* allocate_aligned_complexes(int n)
{
    return static_cast<std::complex<double>*>(
        fftw_malloc(static_cast<size_t>(n) * sizeof(std::complex<double>))
    );
}

inline void free_aligned_mem(void* p)
{
    fftw_free(p);
}

#endif
