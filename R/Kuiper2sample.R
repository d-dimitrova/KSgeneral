Kuiper2sample <- function(x, y, conservative = F, tol = 1e-08) {
  Varname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  Method2 <- NULL
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  Nx <- length(x)
  Ny <- length(y)
  if (Nx < 1 || Ny < 1) {
    stop("not enough 'x' or 'y' data")
  }
  
  EDF1 <- ecdf(x)
  EDF2 <- ecdf(y)
  N <- Nx + Ny
  Joint <- sort(c(x, y))
  Crit <- unique(Joint)
  M <- as.numeric(table(Joint))
  
  z <- EDF1(Crit) - EDF2(Crit)
  DSTAT <- max(z) + max(-z)
  
  
  if (length(M) == N) {
    Method1 <- "Without Tie"
  } else {
    Method1 <- "With Ties"
    if (conservative) {
      Method1 <- paste(Method1, "(Upper Bound P-value)")
      M <- rep(1, N)
    }
  }
  
  names(DSTAT) <- "D"
  Method <- paste("Two-sample Kuiper Test", Method1, Method2)
  result <- Kuiper2sample_Rcpp(Nx, Ny, M, DSTAT, tol)
  
  if(result < -2.5){
    stop("Calculation unstable")
  }
  
  Kuiper2samp <- list(p.value = result, method = Method, statistic = DSTAT, alternative = "two-sided",
                     data.name = Varname)
  class(Kuiper2samp) = "htest"
  return(Kuiper2samp)
}
