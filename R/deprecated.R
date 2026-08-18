############################################################
## Deprecated interfaces retained for backward compatibility

ks_c_cdf_Rcpp <- function(n)
{
  .Deprecated("ks_c_cdf", package = "KSgeneral")
  .ks_c_cdf_Rcpp_legacy(n)
}
