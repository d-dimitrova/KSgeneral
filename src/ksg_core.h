#ifndef KSG_CORE_H
#define KSG_CORE_H

int gcd(int m, int n);
double ks2sample_c_cpp(int nx, int ny, int kind, const int M[], int lengthM,
                       double q, const double w_vec[], int lengthw,
                       double tol);
double ks2sample_cpp(int nx, int ny, int kind, const int M[], int lengthM,
                     double q, const double w_vec[], int lengthw, double tol);
double kuiper2sample_c_cpp(int nx, int ny, const int M[], int lengthM,
                           double q);
double kuiper2sample_cpp(int nx, int ny, const int M[], int lengthM, double q);

double cont_ks_distribution(long n);

#endif
