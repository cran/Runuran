/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#ifndef INFINITY
#  if defined(HAVE_DECL_HUGE_VAL)
      const double INFINITY = HUGE_VAL;
#  elif defined(HAVE_DIVIDE_BY_ZERO)
      const double INFINITY = 1.0 / 0.0;
#  elif defined(HAVE_DECL_DBL_MAX)
      const double INFINITY = DBL_MAX;
#  else
#     error
#     error +--------------------------------------------+
#     error ! Sorry, Cannot define INFINITY correctly!.. !
#     error ! Please contact <unuran@statmath.wu.ac.at>. !
#     error +--------------------------------------------+
#     error
#  endif
#endif
#define ARCMEAN_HARMONIC   (1.e3)   
#define ARCMEAN_ARITHMETIC (1.e-6)  
double
_unur_arcmean( double x0, double x1 )
{
  double a0,a1;
  double r;
  if (x0>x1) {double tmp = x0; x0=x1; x1=tmp;}
  if (x1 < -ARCMEAN_HARMONIC || x0 > ARCMEAN_HARMONIC)
    return (2./(1./x0 + 1./x1));
  a0 = (x0<=-UNUR_INFINITY) ? -M_PI/2. : atan(x0);
  a1 = (x1>= UNUR_INFINITY) ?  M_PI/2. : atan(x1);
  if (fabs(a0-a1) < ARCMEAN_ARITHMETIC)
    r = 0.5*x0 + 0.5*x1;
  else
    r = tan((a0 + a1)/2.);
  return r;
} 
