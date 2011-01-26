/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include "cephes_source.h"
double _unur_cephes_polevl( double x, double coef[], int N )
{
double ans;
int i;
double *p;
p = coef;
ans = *p++;
i = N;
do
	ans = ans * x  +  *p++;
while( --i );
return( ans );
}
double _unur_cephes_p1evl( double x, double coef[], int N )
{
double ans;
double *p;
int i;
p = coef;
ans = x + *p++;
i = N-1;
do
	ans = ans * x  + *p++;
while( --i );
return( ans );
}
