KS2sample <- function(x, y, alternative = c("two.sided", "less", "greater"), 
                     conservative = F, weight = 0, tol = 1e-08) {
  Varname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  Method2 <- NULL
  if (is.numeric(weight)) {
    if (weight > 0 && weight <= 1) {
      WEIGHTEDFUN <- function(t) {
        return((t * (1 - t))^(-weight))
      }
      Method <- "Weighted"
      Method2 <- paste("(v=", as.character(weight), ")", sep = "")
    } else if (weight == 0) {
      WEIGHTEDFUN = function(t) {
        return(1)
      }
      Method <- "Unweighted"
    } else {
      stop("Please enter legal index of weight function")
    }
  } else if (is.function(weight)) {
    if (length(formals(weight)) != 1) {
      stop("The weight function must be a unary function")
    } else {
      WEIGHTEDFUN <- weight
      Method <- "User Defined Weighted"
    }
  } else {
    stop("Please enter legal index of weight function 
  or enter strictly positive binary function as weight")
  }
  
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
  
  W_vec <- unlist(lapply(1:(N - 1)/N, WEIGHTEDFUN))
  if (length(W_vec) < length(W_vec[W_vec > 0])) {
    stop("Weight function is not strictly positive")
  }
  if (length(W_vec) < length(W_vec[!is.na(W_vec)])) {
    stop("Weight function is not well defined")
  }
  
  M <- NULL
  Joint <- sort(c(x, y))
  Crit <- unique(Joint)
  M <- as.numeric(table(Joint))
  
  z <- EDF1(Crit) - EDF2(Crit)
  z <- z[-length(z)]
  z <- z * W_vec[cumsum(M)[-length(M)]]
  
  alternative <- match.arg(alternative)
  
  if (alternative == "two.sided") {
    KIND <- 1
    DSTAT <- max(abs(z))
    Alter <- "two-sided"
  } else if (alternative == "greater") {
    DSTAT <- max(z)
    KIND <- 2
    Alter <- "greater"
    if (Nx != min(Nx, Ny)) {
      KIND <- 3
      p <- Nx
      Nx <- Ny
      Ny <- p
    }
  } else if (alternative == "less") {
    DSTAT <- max(-z)
    KIND <- 3
    Alter <- "less"
    if (Nx != min(Nx, Ny)) {
      KIND <- 2
      p <- Nx
      Nx <- Ny
      Ny <- p
    }
  }
  
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
  Method <- paste(Method, "Two-sample Kolmogorov-Smirnov Test", Method2, Method1)  
  result <- KS2sample_Rcpp(Nx, Ny, KIND, M, DSTAT, W_vec, tol)
  if(result < -2.5){
    stop("Calculation unstable")
  }
  KS2samp <- list(p.value = result, method = Method, statistic = DSTAT, alternative = Alter,
                 data.name = Varname)
  class(KS2samp) <- "htest"
  return(KS2samp)
}





