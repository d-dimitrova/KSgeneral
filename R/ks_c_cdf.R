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

  ## Elementwise validation (finite values, range and monotonicity) and
  ## first-crossing detection are combined in one C++ pass.
  ##
  ## Preserve the historical internal boundary ordering:
  ## upper boundary B first, lower boundary A second.
  ans <- .ks_c_cdf_direct(n, B, A)

  if (ans$crossing >= 0L) {
    warning(
      sprintf(
        "The lower and upper boundaries cross! i=%d.",
        ans$crossing
      ),
      call. = FALSE
    )
  }

  ans$value
}
