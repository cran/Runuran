/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#ifdef HAVE_LIBRMATH
#include <unur_Rmath.h>
double
_unur_sf_bessel_k(double x, double nu)
{
  return bessel_k( x, nu, 1);
} 
double
_unur_sf_bessel_k_expo(double x, double nu)
{
  return bessel_k( x, nu, 2);
} 
#else
double _unur_sf_bessel_k(double x ATTRIBUTE__UNUSED, double nu ATTRIBUTE__UNUSED) 
{ 
  return 0./0.; 
}
#endif 
