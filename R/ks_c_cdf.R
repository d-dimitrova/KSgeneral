############################################################
## Direct boundary interface for the one-sample Exact-KS-FFT method

ks_c_cdf <- function(n, A, B)
{
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || !is.finite(n) ||
      n <= 0 || n != floor(n)) {
    stop("'n' must be a positive integer")
  }

  if (!is.numeric(A) || !is.numeric(B)) {
    stop("'A' and 'B' must be numeric vectors")
  }

  if (length(A) != n || length(B) != n) {
    stop("'A' and 'B' must both have length 'n'")
  }

  if (anyNA(A) || anyNA(B) || any(!is.finite(A)) || any(!is.finite(B))) {
    stop("'A' and 'B' must contain only finite, non-missing values")
  }

  if (any(A < 0 | A > 1) || any(B < 0 | B > 1)) {
    stop("'A' and 'B' must contain values in [0, 1]")
  }

  if (is.unsorted(A) || is.unsorted(B)) {
    stop("'A' and 'B' must be nondecreasing")
  }

  # The internal C++ routine preserves the historical boundary order
  # Preserve the historical boundary ordering: upper B first, lower A second.
  .ks_c_cdf_direct(n, as.numeric(B), as.numeric(A))
}
