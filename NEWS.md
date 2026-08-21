# KSgeneral 2.1.0

## New interface

- Added `ks_c_cdf(n, A, B)`, a documented low-level interface to the
  Exact-KS-FFT boundary-crossing calculation. It accepts the lower (`A`) and
  upper (`B`) boundary vectors directly and does not use an intermediate file.

## C++ source interface

- Added `src/KSgeneral.h` and `src/KSgeneral.cpp` as the source-level C++
  interface. The header uses only standard C++ types and can be included by
  users compiling the KSgeneral numerical sources directly.
- Centralized all `Rcpp::export` attributes in `src/KSgeneral.cpp`; the
  numerical implementation files no longer contain Rcpp wrappers.
- Organized the numerical entry points as `src/k1sample.cpp`/`src/k1sample.h`
  for the one-sample implementation and `src/k2sample.cpp` for the two-sample
  KS and Kuiper implementations.

## Internal changes

- Refactored the one-sample Exact-KS-FFT calculation so package functions pass
  boundary vectors directly to C++ instead of writing and reading
  `Boundary_Crossing_Time.txt`.
- Removed the package-internal working-directory side effect and the shared
  temporary-file dependency from continuous, discrete, and mixed one-sample
  calculations.
- Deprecated the historical file-based `ks_c_cdf_Rcpp(n)` interface. It is
  retained temporarily for backward compatibility and now issues a deprecation
  warning directing users to `ks_c_cdf(n, A, B)`.
