/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <config.h>
#include <math.h>
#include "unur_specfunct_source.h"
#if !HAVE_DECL_LOG1P
double _unur_log1p (double x)
{
  volatile double y;
  y = 1 + x;
  return log(y) - ((y-1)-x)/y ;  
}
#endif
