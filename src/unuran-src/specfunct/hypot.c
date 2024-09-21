/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#if !HAVE_DECL_HYPOT
double _unur_hypot (const double x, const double y)
{
  double xabs = fabs(x) ;
  double yabs = fabs(y) ;
  double min, max;
  if (xabs < yabs) {
    min = xabs ;
    max = yabs ;
  } else {
    min = yabs ;
    max = xabs ;
  }
  if (min == 0) 
    {
      return max ;
    }
  {
    double u = min / max ;
    return max * sqrt (1 + u * u) ;
  }
} 
#endif
