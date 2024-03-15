![downloads](https://cranlogs.r-pkg.org/badges/grand-total/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/KSgeneral)
![downloads](https://cranlogs.r-pkg.org/badges/last-week/KSgeneral)
[![Rdoc](http://www.rdocumentation.org/badges/version/KSgeneral)](http://www.rdocumentation.org/packages/KSgeneral)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/KSgeneral)](https://cran.r-project.org/package=KSgeneral)

# KSgeneral
Computes p-values for the one-sample and two-sample Kolmogorov-Smirnov (KS) tests and the two-sample Kuiper test for any fixed critical level and arbitrary (possibly very large) sample sizes. For the one-sample KS test, it allows the pre-specified cumulative distribution function under the null hypothesis to be continuous, purely discrete or mixed. For the two-sample test, it is assumed that both samples come from an unspecified (unknown) continuous, purely discrete or mixed distribution, i.e. ties (repeated observations) are allowed.

If a data sample is supplied, 'KSgeneral' (also available from  https://CRAN.R-project.org/package=KSgeneral) computes the p-value corresponding to the value of the KS test statistic computed based on the user provided data sample.

The functions in 'KSgeneral' for the one-sample KS test implement a novel, accurate and efficient method named Exact-KS-FFT, developed by Dimitrova, Kaishev, Tan (2017), available together with the underlying C++ code from http://openaccess.city.ac.uk/18541.

The functions in 'KSgeneral' for the two-sample test implement algorithms which generalize the method due to Nikiforov (1994), and calculate the exact p-values of the KS test and the Kuiper test respectively. Both of them allow tested data samples to come from continuous, discrete or mixed distributions (ties are also allowed).

**To cite this package in publication: (for the use of the one-sample KS test) Dimitrina S. Dimitrova, Vladimir K. Kaishev, and Senren Tan. Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed, or Continuous. *Journal of Statistical Software*. 2020, 95(10): 1–42. <doi:10.18637/jss.v095.i10>, 
(for the use of the two-sample KS and Kuiper tests) Dimitrina S. Dimitrova, Vladimir K. Kaishev, and Yun Jia (2024). The R functions KS2sample and Kuiper2sample: Efficient Exact Calculation of P-values of the Two-sample Kolmogorov-Smirnov and Kuiper Tests. *submitted*.** 

The p-value for the one-sample KS test is expressed as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using Fast Fourier Transform (FFT). The p-values for the two-sample KS and Kuiper tests are expressed as the ratio of the total numbers for point sequences defined on an integer-valued grid stay wholly in a subset to a combinatorial number.

The package can also be used to compute and plot the complementary cdf of the one-sample KS statistic which is known to depend on the hypothesized distribution when the latter is discontinuous (i.e. purely discrete or mixed).


# Installation
In order to build the KSgeneral package from source, a C++ compiler is required. 

The latter is contained in the Windows Rtools, available from https://cran.r-project.org/bin/windows/Rtools/, or under MacOS in Xcode, downloadable from the App Store.

The package KSgeneral uses Rcpp in R, and utilizes the C++ code that efficiently computes the complementary cdf using the Exact-KS-FFT method developed by Dimitrova, Kaishev, Tan (2017), available together with the underlying C++ code from http://openaccess.city.ac.uk/18541 and the C++ code which generalize the Fortran subroutine due to Nikiforov (1994).

Since the Exact-KS-FFT method requires computation of Fast Fourier Transform (FFT), the FFTW3 library developed by Matteo Frigo and Steven G.Johnson needs to be installed from http://www.fftw.org/index.html. 

It should be noted that the Rtools and FFTW3 should be installed in the system PATH.

For Windows users, The FFTW3 library (static library, with a ".a" extension) for Windows (32-bit or 64-bit) can be found in the local323.zip file, available from http://www.stats.ox.ac.uk/pub/Rtools/libs.html.

For Mac or Unix users, it is straightforward to install the FFTW3 library from the command line, following the instructions from http://www.fftw.org/index.html.
