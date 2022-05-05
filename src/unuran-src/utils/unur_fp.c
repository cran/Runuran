/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
int
_unur_FP_cmp( double x1, double x2, double eps)
{
  double fx1 = (x1>=0.) ? x1 : -x1;
  double fx2 = (x2>=0.) ? x2 : -x2;
  double delta = eps * _unur_min(fx1,fx2);
  double difference = x1 - x2;
  if (_unur_isinf(delta)) {
    delta = eps * DBL_MAX;
  }
  if (fx1 <= 2.*DBL_MIN && fx2 <= 2.*DBL_MIN)
    return 0;
  if (difference > delta)       
    return +1;
  else if (difference < -delta) 
    return -1;
  else                          
    return 0;                   
} 
#ifndef _unur_iszero
int _unur_iszero (const double x)
{
  return (x==0.);
} 
#endif
#ifndef _unur_isone
int _unur_isone (const double x)
{
  return (x==1.);
} 
#endif
#ifndef _unur_isfsame
int _unur_isfsame (const double x, const double y)
{
  return (x==y);
} 
#endif
int
_unur_isfinite (const double x)
{
#if HAVE_DECL_ISFINITE
  return (isfinite(x) ? TRUE : FALSE);
#elif defined(_MSC_VER) 
  return (_finite(x) ? TRUE : FALSE);
#elif HAVE_IEEE_COMPARISONS
  if (x < UNUR_INFINITY && x > -UNUR_INFINITY)
    return TRUE;
  else
    return FALSE;
#else
# error
# error +--------------------------------------------+
# error ! Sorry, Cannot handle INFINITY correctly! . !
# error ! Please contact <unuran@statmath.wu.ac.at>. !
# error +--------------------------------------------+
# error
#endif
} 
int
_unur_isnan (const double x)
{
#if HAVE_DECL_ISNAN
  return (isnan(x) ? TRUE : FALSE);
#elif defined(_MSC_VER) 
  return (_isnan(x) ? TRUE : FALSE);
#elif HAVE_IEEE_COMPARISONS
  return ((x!=x) ? TRUE : FALSE);
#else
# error
# error +--------------------------------------------+
# error ! Sorry, Cannot handle NaN correctly! ...... !
# error ! Please contact <unuran@statmath.wu.ac.at>. !
# error +--------------------------------------------+
# error
#endif
} 
int
_unur_isinf (const double x)
{
#if HAVE_DECL_ISINF
  return isinf(x);
#elif defined(_MSC_VER) 
  int fpc = _fpclass(x);
  if (fpc == _FPCLASS_PINF)
    return +1;
  else if (fpc == _FPCLASS_NINF)
    return -1;
  else 
    return 0;
#elif HAVE_IEEE_COMPARISONS
  if (x>=UNUR_INFINITY)
    return 1;
  else if (x<=-UNUR_INFINITY)
    return -1;
  else
    return 0;
#else
# error
# error +--------------------------------------------+
# error ! Sorry, Cannot handle INFINITY correctly! . !
# error ! Please contact <unuran@statmath.wu.ac.at>. !
# error +--------------------------------------------+
# error
#endif
} 
