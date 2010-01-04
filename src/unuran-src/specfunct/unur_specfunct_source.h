/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_SPECFUNCT_SOURCE_H_SEEN
#define UNUR_SPECFUNCT_SOURCE_H_SEEN
double _unur_cephes_incbet(double a, double b, double x);
#define _unur_sf_incomplete_beta(x,a,b)   _unur_cephes_incbet((a),(b),(x))
double _unur_cephes_lgam(double x);
#define _unur_sf_ln_gamma(x)   _unur_cephes_lgam(x)
#define _unur_sf_ln_factorial(x)   _unur_sf_ln_gamma((x)+1.)
double _unur_cephes_igam(double a, double x);
#define _unur_sf_incomplete_gamma(x,a)  _unur_cephes_igam((a),(x))
double _unur_cephes_ndtr(double x);
#define _unur_sf_cdfnormal(x)   _unur_cephes_ndtr(x)
double _unur_cephes_ndtri(double x);
#define _unur_sf_inv_cdfnormal(x)   _unur_cephes_ndtri(x)
#ifdef HAVE_LIBRMATH
#define HAVE_BESSEL_K 1
double _unur_sf_bessel_k(double x, double nu);
double _unur_sf_bessel_k_expo(double x, double nu);
#endif
#if !HAVE_DECL_LOG1P
double _unur_log1p(double x);
#endif
#endif  
