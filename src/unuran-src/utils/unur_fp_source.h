/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_FP_SOURCE_H_SEEN
#define UNUR_FP_SOURCE_H_SEEN
int _unur_FP_cmp( double x1, double x2, double eps);
#define _unur_FP_cmp_same(a,b) (_unur_FP_cmp((a),(b),DBL_EPSILON))
#define _unur_FP_cmp_equal(a,b) (_unur_FP_cmp((a),(b),UNUR_EPSILON))
#define _unur_FP_cmp_approx(a,b) (_unur_FP_cmp((a),(b),UNUR_SQRT_DBL_EPSILON))
#define _unur_FP_same(a,b) (_unur_FP_cmp((a),(b),DBL_EPSILON)==0)
#define _unur_FP_equal(a,b) (_unur_FP_cmp((a),(b),UNUR_EPSILON)==0)
#define _unur_FP_approx(a,b) (_unur_FP_cmp((a),(b),UNUR_SQRT_DBL_EPSILON)==0)
#define _unur_FP_less(a,b) ((_unur_FP_cmp((a),(b),UNUR_EPSILON)<0) ? TRUE : FALSE)
#define _unur_FP_greater(a,b) ((_unur_FP_cmp((a),(b),UNUR_EPSILON)>0) ? TRUE : FALSE)
#define _unur_iszero(x)     ((x)==0.0)
#define _unur_isone(x)      ((x)==1.0)
#define _unur_isfsame(x,y)  ((x)==(y))
#ifndef _unur_iszero
int _unur_iszero (const double x);
#endif
#ifndef _unur_isone
int _unur_isone (const double x);
#endif
#ifndef _unur_isfsame
int _unur_isfsame (const double x, const double y);
#endif
int _unur_isfinite (const double x);
int _unur_isnan (const double x);
int _unur_isinf (const double x);
#define _unur_FP_is_infinity(a)  ((a) >= INFINITY)
#define _unur_FP_is_minus_infinity(a)  ((a) <= -INFINITY)
#endif  
