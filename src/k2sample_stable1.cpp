#include<Rcpp.h>
#include<algorithm>
#include<numeric>
#include<cmath>
 


const double one = 1.0, zero = 0.0, small = 1e-35, smalln = - 80.5904782547916, chknum = 1e32, aln2 = 0.69314718056, smallrev = 1e35;
const int iterup = 116, int_limit = 2147483647;

using namespace std;
using namespace Rcpp;
double ks2sample_c_cpp(int nx, int ny, int kind, int M[], int lengthM, double q, double w_vec[], int lengthw, double tol){
// Based on ALGORITHM AS 288 APPL.STATIST. (1994), VOL.43, NO.1
// P-value calculation for the weighted generalized two-sample Kolmogorov-Smirnov tests.
// The tests are conditional on ties in the pooled sample.
   int i(0), ic(0), icat(1), ile(0), ile2(0), ileft(0), iles(0), iri(0), iri2(0), iright(0), iris(0), jlow(0), jupp(0), l(0), l2(0), n(0), nofdiv(0), nties(0);
   double slope(0), delta(0), deviat(0), dl(0), scl(0), pval(0), eps(1e-6);
   bool newrct = true;

   double* P = new double[nx+2]();
   eps = tol;
   if((tol<=0) || (tol>1e-6)) eps = 1e-6;
   scl = one;
   delta = q - eps;
// Check illegal output
   if((nx<1) || (ny<1) || (kind<1) || (kind>3) || (lengthw!=nx+ny-1)) return -1;
   if((accumulate(M,M+lengthM,0) != nx+ny) || (*min_element(M,M+lengthM)<1) || (*min_element(w_vec,w_vec+lengthw)<eps)) return -2.0;
   if(delta<zero) return 1.0;
   P[0] = one;
// Parameters to define a set for trajectories to lie within it
   n = nx + ny;
   slope = nx / (double)(n);
   delta = slope * delta * ny;
   nties = M[0];
// Variables to prevent from overflows in P
   ic = iterup;
   scl = one;
   for(l=1; l<n; l++)
   {
   		if(nties == 1)
		{
//Calculate boundaries for current L (if Lth value in the pooled sample is unique or `last' in a series of tied values)
   			dl = l * slope;
   			deviat = delta / w_vec[l-1];
   			iri = min(nx, min( int(dl + deviat), l));
   			ile = max(0, max( int(dl - deviat + one), l - ny));
   			nties = M[icat];
   			icat = icat + 1;
   			newrct = true;
	   }
	   else
	   {
//Calculations for tied observations
	   	    nties = nties-1;
//If we have the first observation with the new value (that is not unique), then determine the ICATth `rectangle'
	   	    if(newrct)
	   	    {
	   	    	newrct =  false;
	   	    	l2 = l + nties;
				dl = l2 * slope;
				if(l2==n)
				{
					iri2 = nx;
					ile2 = nx;
				}
				else
				{
					deviat = delta / w_vec[l2-1];
// X axis boundaries of the subset of `line' L(l+1) within ICATth rectangle
					iri2 = min( min( int(dl + deviat), l2), nx);
					ile2 = max( max( int(dl - deviat + one) , l2 - ny), 0);
				}
// Four sides of the rectangle on the X and Y axes
				ileft = ile;
				iright = iri2;
				jupp = l2 - ile2;
				jlow = l - iri - 1;
			}
// Calculate boundaries for current L (Lth value is tied)
			ile = max(ileft, l - jupp);
			iri = min(iright, l - jlow);
	   }
// Set the left (right) boundary for the one-sided test
		if(kind == 3) {
			iri = min(nx , l);
		} else if (kind == 2) {
			ile = max(0 , l - ny);
		}
// Calculate the number of trajectories p(i,j) for current L
		iles = max(1, ile);
		iris = min(l-1, iri); 
		for(i = iris; i >= iles; i--) P[i] = P[i] + P[i-1];
// Check whether elements of P are large enough to multiply them by SMALL 
		ic = ic -1;
		if(ic <= 0)
		{
			dl = zero;
			for(i = iles; i <= iris; i++){
				if(P[i] < 0) return -3.0;
				dl = max(P[i],dl);
			} 
			if(dl == zero) return 1.0;
			if(dl > chknum)
			{
				for(i = iles; i<= iris; i++) P[i] = P[i] * small;
				ic = iterup;
				nofdiv = nofdiv + 1;
				scl = scl * small;
			}
			else
			{ 
// Estimate the number of iterations for DL=Pmax to became of order 1/SMALL
				ic = (-smalln - log(dl)) / aln2;
			}
		}
// Define, whether boundaries lie on the left and lower sides of the rectangle R and define boundary values for the next iterations
		P[iles - 1] = (ile == 0) ? scl : zero;
		P[iris + 1] = (iri == l) ? scl : zero;
   }
   dl = P[nx] + P[nx - 1];
   delete [] P;
   if (dl == zero) return 1.0;
   if(nofdiv == 0)
   {
   		for(l = min(nx,ny); l > 0 ; l--) dl *= double(l)/(l + max(nx,ny) );
   } else {
   		for(l = min(nx,ny); l > 0 ; l--) 
		{
			dl *= double(l)/(l + max(nx,ny) );
			if( dl < 1 && nofdiv > 0)
			{
				dl *= smallrev;
				nofdiv--;
			} 
		}
		for(;nofdiv > 0; nofdiv--) dl *= smallrev;
   }
   pval = dl;
   
   return pval; 
   
}

double ks2sample_cpp(int nx, int ny, int kind, int M[], int lengthM, double q, double w_vec[], int lengthw, double tol){
  // Based on ALGORITHM AS 288 APPL.STATIST. (1994), VOL.43, NO.1
  // P-value calculation for the weighted generalized two-sample Kolmogorov-Smirnov tests.
  // The tests are conditional on ties in the pooled sample.
  int i(0), icat(1), ile(0), ile2(0), ileft(0), iles(0), iri(0), iri2(0), iright(0), iris(0), jlow(0), jupp(0), l(0), l2(0), n(0), nties(0);
  double slope(0), delta(0), deviat(0), dl(0), pval(0), eps(1e-6);
  bool newrct = true;
  
  double* P = new double[nx+2]();
  double* Q = new double[nx+2]();
  eps = tol;
  if((tol<=0) || (tol>1e-6)) eps = 1e-6;
  delta = q - eps;
  // Check illegal output
  if((nx<1) || (ny<1) || (kind<1) || (kind>3) || (lengthw!=nx+ny-1)) return -1;
  if((accumulate(M,M+lengthM,0) != nx+ny) || (*min_element(M,M+lengthM)<1) || (*min_element(w_vec,w_vec+lengthw)<eps)) return -2.0;
  if(delta<zero) return 1.0;
  P[0] = one;
  // Parameters to define a set for trajectories to lie within it
  n = nx + ny;
  slope = nx / (double)(n);
  delta = slope * delta * ny;
  nties = M[0];
  // Variables to prevent from overflows in P


  for(l=1; l<n; l++)
  {
    if(nties == 1)
    {
      //Calculate boundaries for current L (if Lth value in the pooled sample is unique or `last' in a series of tied values)
      dl = l * slope;
      deviat = delta / w_vec[l-1];
      iri = min(nx, min( int(dl + deviat), l));
      ile = max(0, max( int(dl - deviat + one), l - ny));
      nties = M[icat];
      icat = icat + 1;
      newrct = true;
    }
    else
    {
      //Calculations for tied observations
      nties = nties-1;
      //If we have the first observation with the new value (that is not unique), then determine the ICATth `rectangle'
      if(newrct)
      {
        newrct =  false;
        l2 = l + nties;
        dl = l2 * slope;
        if(l2==n)
        {
          iri2 = nx;
          ile2 = nx;
        }
        else
        {
          deviat = delta / w_vec[l2-1];
          // X axis boundaries of the subset of `line' L(l+1) within ICATth rectangle
          iri2 = min( min( int(dl + deviat), l2), nx);
          ile2 = max( max( int(dl - deviat + one) , l2 - ny), 0);
        }
        // Four sides of the rectangle on the X and Y axes
        ileft = ile;
        iright = iri2;
        jupp = l2 - ile2;
        jlow = l - iri - 1;
      }
      // Calculate boundaries for current L (Lth value is tied)
      ile = max(ileft, l - jupp);
      iri = min(iright, l - jlow);
    }
    // Set the left (right) boundary for the one-sided test
    if(kind == 3) {
      iri = min(nx , l);
    } else if (kind == 2) {
      ile = max(0 , l - ny);
    }
    // Calculate the proportion of trajectories J(i,j) for current L
    iles = max(1, ile);
    iris = min(l-1, iri); 
    for(i = iris; i >= iles; i--) Q[i - iles] = ((l - i)*P[i] +  i * P[i - 1])/double(l);
    for(i = iris; i >= iles; i--) P[i] = Q[i - iles];
    // Define, whether boundaries lie on the left and lower sides of the rectangle R and define boundary values for the next iterations
    P[iles - 1] = (ile == 0) ? zero : one;
    P[iris + 1] = (iri == l) ? zero : one;
  }
  delete [] Q;
  pval = (ny * P[nx] + nx * P[nx - 1])/ (double)(n);
  delete [] P;
  return pval; 
}


double kuiperks_p(int nx, int ny, int M[], int lengthM, double dstatup, double dstatdown , double tol){
// Based on ALGORITHM AS 288 APPL.STATIST. (1994), VOL.43, NO.1
// P-value calculation for the weighted generalized two-sample Kolmogorov-Smirnov tests.
// The tests are conditional on ties in the pooled sample.
   int i(0), icat(1), ile(0), ile2(0), ileft(0), iles(0), iri(0), iri2(0), iright(0), iris(0), jlow(0), jupp(0), l(0), l2(0), n(0), nties(0);
   double slope(0), deltalow(0), deltaup(0), deviatlow(0), deviatup(0), pval(0), dl(0);
   bool newrct = true;
   double* P = new double[nx+2]();
   double* Q = new double[nx+2]();
  
   P[0] = one;
   n = nx + ny;
   deltaup = dstatup + tol;
   deltalow = dstatdown + tol;
// Parameters to define a set for trajectories to lie within it
   slope = nx / (double)(n);
   deviatup = slope * deltaup * ny;
   deviatlow = slope * deltalow * ny;
   nties = M[0];
// Variables to prevent from overflows in P
   for(l=1; l<n; l++)
   {
   		if(nties == 1)
		{
//Calculate boundaries for current L (if Lth value in the pooled sample is unique or `last' in a series of tied values)
   			dl = l * slope;
   			iri = min(nx, min( int(dl + deviatlow), l));
   			ile = max(0, max( int(dl - deviatup + one), l - ny));
   			nties = M[icat];
   			icat = icat + 1;
   			newrct = true;
	   }
	   else
	   {
//Calculations for tied observations
	   	    nties = nties-1;
//If we have the first observation with the new value (that is not unique), then determine the ICATth `rectangle'
	   	    if(newrct)
	   	    {
	   	    	newrct =  false;
	   	    	l2 = l + nties;
				dl = l2 * slope;
// X axis boundaries of the subset of `line' L(l+1) within ICATth rectangle
				iri2 = min( min( int(dl + deviatlow), l2), nx);
				ile2 = max( max( int(dl - deviatup + one) , l2 - ny), 0);
// Four sides of the rectangle on the X and Y axes
				ileft = ile;
				iright = iri2;
				jupp = l2 - ile2;
				jlow = l - iri - 1;
			}
// Calculate boundaries for current L (Lth value is tied)
			ile = max(ileft, l - jupp);
			iri = min(iright, l - jlow);
	   }
// Calculate the proportion of trajectories J(i,j) for current L
		iles = max(1, ile);
		iris = min(l-1, iri); 
		for(i = iris; i >= iles; i--) Q[i - iles] = ((l - i)*P[i] +  i * P[i - 1])/double(l);
    	for(i = iris; i >= iles; i--) P[i] = Q[i - iles];
// Define, whether boundaries lie on the left and lower sides of the rectangle R and define boundary values for the next iterations
    	P[iles - 1] = (ile == 0) ? zero : one;
    	P[iris + 1] = (iri == l) ? zero : one;
   }
  delete [] Q;
  pval = (ny * P[nx] + nx * P[nx - 1])/ (double)(n);
  delete [] P;
  return pval; 
}


double kuiperks_n(int nx, int ny, int M[], int lengthM, double dstatup, double dstatdown , double tol){
// Based on ALGORITHM AS 288 APPL.STATIST. (1994), VOL.43, NO.1
// counts the number of trajectories with their one-sided KS statistics (supremum and infimum) equal or smaller than dstatup and dstatdown respectively;
   int i(0), ic(0), icat(1), ile(0), ile2(0), ileft(0), iles(0), iri(0), iri2(0), iright(0), iris(0), jlow(0), jupp(0), l(0), l2(0), n(0), nofdiv(0), nties(0);
   double slope(0), deltalow(0), deltaup(0), deviatlow(0), deviatup(0), scl(0), nval(0), dl(0);
   bool newrct = true;
   double* P = new double[nx+2]();
   
   P[0] = one;
   n = nx + ny;
   deltaup = dstatup + tol;
   deltalow = dstatdown + tol;
// Parameters to define a set for trajectories to lie within it
   slope = nx / (double)(n);
   deviatup = slope * deltaup * ny;
   deviatlow = slope * deltalow * ny;
   nties = M[0];
// Variables to prevent from overflows in P
   ic = iterup;
   scl = one;
   for(l=1; l<n; l++)
   {
   		if(nties == 1)
		{
//Calculate boundaries for current L (if Lth value in the pooled sample is unique or `last' in a series of tied values)
   			dl = l * slope;
   			iri = min(nx, min( int(dl + deviatlow), l));
   			ile = max(0, max( int(dl - deviatup + one), l - ny));
   			nties = M[icat];
   			icat = icat + 1;
   			newrct = true;
	   }
	   else
	   {
//Calculations for tied observations
	   	    nties = nties-1;
//If we have the first observation with the new value (that is not unique), then determine the ICATth `rectangle'
	   	    if(newrct)
	   	    {
	   	    	newrct =  false;
	   	    	l2 = l + nties;
				dl = l2 * slope;
// X axis boundaries of the subset of `line' L(l+1) within ICATth rectangle
				iri2 = min( min( int(dl + deviatlow), l2), nx);
				ile2 = max( max( int(dl - deviatup + one) , l2 - ny), 0);
// Four sides of the rectangle on the X and Y axes
				ileft = ile;
				iright = iri2;
				jupp = l2 - ile2;
				jlow = l - iri - 1;
			}
// Calculate boundaries for current L (Lth value is tied)
			ile = max(ileft, l - jupp);
			iri = min(iright, l - jlow);
	   }
// Set the left (right) boundary for the one-sided test
// Calculate the number of trajectories p(i,j) for current L
		iles = max(1, ile);
		iris = min(l-1, iri); 
		for(i = iris; i >= iles; i--) P[i] = P[i] + P[i-1];
// Check whether elements of P are large enough to multiply them by SMALL 
		ic = ic -1;
		if(ic <= 0)
		{
			dl = zero;
			for(i = iles; i <= iris; i++){
				if(P[i] < 0) return -3.0;
				dl = max(P[i],dl);
			} 
			if(dl == zero) return -1.0;
			if(dl > chknum)
			{
				for(i = iles; i<= iris; i++) P[i] = P[i] * small;
				ic = iterup;
				nofdiv = nofdiv + 1;
				scl = scl * small;
			}
			else
			{ 
// Estimate the number of iterations for DL=Pmax to became of order 1/SMALL
				ic = (-smalln - log(dl)) / aln2;
			}
		}
// Define, whether boundaries lie on the left and lower sides of the rectangle R and define boundary values for the next iterations
		P[iles - 1] = (ile == 0) ? scl : zero;
		P[iris + 1] = (iri == l) ? scl : zero;
   }
   dl = P[nx] + P[nx - 1];
   if (dl == zero) return -2.0;
   delete [] P;
   nval = log(dl) - nofdiv * smalln;
   return nval;
}


int lcm(int m, int n) 
{
    return m % n == 0 ? n : gcd(n, m % n);
}

double kuiper2sample_c_cpp(int nx, int ny, int M[], int lengthM, double q)
{
	double nxydouble(0), eps(1e-6), pval(0);
	long double lfacto(0), neg(0), pos(0);
 	long double minmum(0);
	int p(0), C(0), i(0), nxy(0);
// Check illegal output
	if((nx<1) || (ny<1)) return -1.0;
    if((accumulate(M,M+lengthM,0)!=nx+ny) || (*min_element(M,M+lengthM)<1)) return -2.0;
	if(q >= 2.0) return 0.0;
// Calculate the least common multiple of nx and ny
	nxy = ny / gcd(nx, ny);
	if(int_limit / nxy <= nx) return -1.0;
	nxy *= nx;
	nxydouble = (double)nxy;
	eps = 0.5 / nxydouble;
	p = nx;
	nx = (nx <= ny) ? nx : ny;
	ny = (p <= ny) ? ny : p;
//	lfacto = lgamma(nx+1) + lgamma(ny+1)- lgamma(nx + ny +1);
	for(i = 1; i <= nx; i++) lfacto += log((long double)i) - log((long double)(ny + i));
	
	C = ceil(q * nxydouble - 1);
	if(C <= 0)
	{
		pval = kuiperks_n(nx, ny, M, lengthM, (double)(0.0), (double)(0.0), eps);
		if(pval < -1.5) return 1.0;
		return exp(pval + lfacto);
	}
	else if((C > 0) && (C <= nxy - 1))
	{

		long double* pval_vec_pos = new long double[C+1]();
		long double* pval_vec_neg = new long double[C]();
		
		pval_vec_pos[C] = kuiperks_n(nx, ny, M, lengthM, C / nxydouble, zero, eps);
		if(pval_vec_pos[C] < -2.5) return -3.0;
        minmum = pval_vec_pos[C];
		for(i = 0; i <= C-1; i++)
		{
// Find  the log number of trajectories with specific bound of Kolmogorov-Smirnov statistic
			pval_vec_pos[i] = kuiperks_n(nx, ny, M, lengthM, i / nxydouble, (C - i) / nxydouble , eps);
			pval_vec_neg[i] = kuiperks_n(nx, ny, M, lengthM, i / nxydouble, (C - i - 1) / nxydouble , eps);
			if((pval_vec_pos[i] < -2.5) || (pval_vec_neg[i] < -2.5)) return -3.0;
// Check whether there are overflows in calculation
            minmum = min(minmum, min( ((pval_vec_pos[i]>0) ? pval_vec_pos[i] : minmum) , ((pval_vec_neg[i]>0) ? pval_vec_neg[i] : minmum)) ) ;
		}
			
// Sum the number of trajectories
		
		lfacto = lfacto + minmum;
		
		pos = exp(pval_vec_pos[C] - minmum);

		for(i = 0; i < (C+1)/2.0 - 1; i++)
		{

		    pos += ((pval_vec_pos[i] >= 0) ? exp(pval_vec_pos[i] - minmum) : zero) + ((pval_vec_pos[C - i - 1] >= 0) ? exp(pval_vec_pos[C - i - 1] - minmum) : zero);
			neg += ((pval_vec_neg[i] >= 0) ? exp(pval_vec_neg[i] - minmum) : zero) + ((pval_vec_neg[C - i - 1] >= 0) ? exp(pval_vec_neg[C - i - 1] - minmum) : zero);
		}
 
		if(C%2!=0){
			pos += (pval_vec_pos[C/2] >= 0) ? exp(pval_vec_pos[C/2] - minmum) : zero;
			neg += (pval_vec_neg[C/2] >= 0) ? exp(pval_vec_neg[C/2] - minmum) : zero;
		}
		delete [] pval_vec_pos;
		delete [] pval_vec_neg;

	}
	else
	{
		p = 2 * nxy - C;
		long double* pval_vec_pos = new long double[p + 1]();
		long double* pval_vec_neg = new long double[p]();
// Find  the log number of trajectories with specific bound of Kolmogorov-Smirnov statistic		
		pval_vec_pos[p] = kuiperks_n(nx, ny, M, lengthM, double(1.0), (C - nxy) / nxydouble , eps);
		if(pval_vec_pos[p] < -2.5) return -3.0;
        minmum = pval_vec_pos[p];
		for(i = 0; i <= p-1; i++)
		{
			pval_vec_pos[i] = kuiperks_n(nx, ny, M, lengthM, (C - nxy + i) / nxydouble, (nxy - i) / nxydouble , eps);
			pval_vec_neg[i] = kuiperks_n(nx, ny, M, lengthM, (C - nxy + i) / nxydouble, (nxy - i - 1) / nxydouble , eps);
// Check whether there are overflows in calculation
			if((pval_vec_pos[i] < -2.5) || (pval_vec_neg[i] < -2.5)) return -3.0;
			minmum = min(minmum, min(pval_vec_pos[i],pval_vec_neg[i]));
		}
		C = p;
		
		lfacto = lfacto + minmum;
// Sum the number of trajectories		
		pos = exp(pval_vec_pos[C] - minmum);

		for(i = 0; i < (C+1)/2.0 - 1; i++)
		{
		    pos += ((pval_vec_pos[i] >= 0) ? exp(pval_vec_pos[i] - minmum) : zero) + ((pval_vec_pos[C - i - 1] >= 0) ? exp(pval_vec_pos[C - i - 1] - minmum) : zero);
			neg += ((pval_vec_neg[i] >= 0) ? exp(pval_vec_neg[i] - minmum) : zero) + ((pval_vec_neg[C - i - 1] >= 0) ? exp(pval_vec_neg[C - i - 1] - minmum) : zero);
		}
 
		if(C%2!=0){
			pos += (pval_vec_pos[C/2] >= 0) ? exp(pval_vec_pos[C/2] - minmum) : zero;
			neg += (pval_vec_neg[C/2] >= 0) ? exp(pval_vec_neg[C/2] - minmum) : zero;
		}
		delete [] pval_vec_pos;
		delete [] pval_vec_neg;
		
	}
	
	if( pval < 0){
		return -3;
	}
	
		
	pval = (pos-neg)*exp(lfacto);
	return pval;
}


double kuiper2sample_cpp(int nx, int ny, int M[], int lengthM, double q)
{
	double nxydouble(0), eps(1e-6), pval(0);
	int p(0), C(0), i(0), nxy(0);
// Check illegal output
	if((nx<1) || (ny<1)) return -1.0;
    if((accumulate(M,M+lengthM,0)!=nx+ny) || (*min_element(M,M+lengthM)<1)) return -2.0;

	if(q >= 2.0) return 0.0;
// Calculate the least common multiple of nx and ny
	nxy = ny / gcd(nx, ny);
	if(int_limit / nxy <= nx) return -1.0;
	nxy *= nx;
	nxydouble = (double)nxy;
	eps = 0.5 / nxydouble;
	p = nx;
	nx = (nx <= ny) ? nx : ny;
	ny = (p <= ny) ? ny : p;
	
	C = ceil(q * nxydouble - 1);
	if(C <= 0)
	{
		pval = kuiperks_p(nx, ny, M, lengthM, (double)(0.0), (double)(0.0), eps);
	}
	else if((C > 0) && (C <= nxy - 1))
	{
// Find and sum the probability with specific bound of Kolmogorov-Smirnov statistic
		for(i = 0; i <= C-1; i++)
		{
			pval += kuiperks_p(nx, ny, M, lengthM, i / nxydouble, (C - i) / nxydouble , eps) - kuiperks_p(nx, ny, M, lengthM, i / nxydouble, (C - i - 1) / nxydouble , eps);
		}
		pval += kuiperks_p(nx, ny, M, lengthM, C / nxydouble, zero, eps);
		return pval;
	}
	else
	{
// Find and sum the probability with specific bound of Kolmogorov-Smirnov statistic
		p = 2 * nxy - C;
		for(i = 0; i <= p-1; i++)
		{
			pval += kuiperks_p(nx, ny, M, lengthM, (C - nxy + i) / nxydouble, (nxy - i) / nxydouble , eps) - kuiperks_p(nx, ny, M, lengthM, (C - nxy + i) / nxydouble, (nxy - i - 1) / nxydouble , eps);

		}
		pval += kuiperks_p(nx, ny, M, lengthM, double(1.0), (C - nxy) / nxydouble , eps);
	}
	

	if( pval < 0){
		return -3;
	}
	return pval;
		

}
 

// [[Rcpp::export]]
double KS2sample_c_Rcpp(int m, int n, int kind, Rcpp::IntegerVector M, double q, Rcpp::NumericVector w_vec, double tol){
	  const int lengthM=M.size();
    const int lengthw=w_vec.size();
    double pval=0;
    pval = ks2sample_c_cpp(m,n,kind,M.begin(),lengthM,q,w_vec.begin(),lengthw,tol);
    return(pval);
}

// [[Rcpp::export]]
double Kuiper2sample_Rcpp(int m, int n, Rcpp::IntegerVector M, double q){
    const int lengthM=M.size();
    double pval = 0;
    pval = kuiper2sample_cpp(m,n,M.begin(),lengthM,q);
    return(pval);
}

// [[Rcpp::export]]
double Kuiper2sample_c_Rcpp(int m, int n, Rcpp::IntegerVector M, double q){
    const int lengthM=M.size();
    double pval = 0;
    pval = kuiper2sample_c_cpp(m,n,M.begin(),lengthM,q);
    return(pval);
}
 
// [[Rcpp::export]]
double KS2sample_Rcpp(int m, int n, int kind, Rcpp::IntegerVector M, double q, Rcpp::NumericVector w_vec, double tol){
    const int lengthM=M.size();
    const int lengthw=w_vec.size();
    double pval=0;
    pval = ks2sample_cpp(m,n,kind,M.begin(),lengthM,q,w_vec.begin(),lengthw,tol);
    return(pval);
}
