/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <float.h>
#include <math.h>
#include <utils/unur_math_source.h>
#include <utils/unur_fp_source.h>
#include <utils/unur_fp_const_source.h>
#define PI      M_PI              
#define SQRTH   M_SQRTH           
#define MAXNUM  DBL_MAX           
double _unur_cephes_gamma( double x );
double _unur_cephes_lgam( double x );
double _unur_cephes_igamc( double a, double x );
double _unur_cephes_igam( double a, double x );
double _unur_cephes_incbet( double aa, double bb, double xx );
double _unur_cephes_ndtr( double a );
double _unur_cephes_erfc( double a );
double _unur_cephes_erf( double x );
double _unur_cephes_ndtri( double y0 );
double _unur_cephes_polevl( double x, double coef[], int N );
double _unur_cephes_p1evl( double x, double coef[], int N );
