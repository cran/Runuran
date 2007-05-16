/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#ifndef INFINITY
const double INFINITY = 1.0 / 0.0;
#endif
#define ARCMEAN_HARMONIC 1.e3  
double
_unur_arcmean( double x0, double x1 )
{
  double x;
  if (x0>x1) {double tmp = x0; x0=x1; x1=tmp;}
  if (x1 < -ARCMEAN_HARMONIC || x0 > ARCMEAN_HARMONIC)
    x = 2./(1./x0 + 1./x1);
  else
    x = tan( (((x0<=-INFINITY) ? -M_PI/2. : atan(x0)) + ((x1>=INFINITY) ? M_PI/2. : atan(x1))) / 2. );
  return x;
} 
