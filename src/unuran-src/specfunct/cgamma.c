/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
double
_unur_Relcgamma (double x, double y)
{
  double t, x0, x1;
  double q1, logq1, q2, th ;
  double gr1,  sr, si;
  int j, k, na;
  static const double a[] = {
    8.333333333333333e-02, -2.777777777777778e-03,
    7.936507936507937e-04, -5.952380952380952e-04,
    8.417508417508418e-04, -1.917526917526918e-03,
    6.410256410256410e-03, -2.955065359477124e-02,
    1.796443723688307e-01, -1.39243221690590 };
  double gr; 
  x1 = 0.0;
  na = 0;
  if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
    return UNUR_INFINITY;
  else if (x < 0.0) {
    x1 = x;
    x = -x;
    y = -y;
  }
  x0 = x;
  if (x <= 7.0) {
    na = (int)(7.0-x);
    x0 = x+na;
  }
  q1 = hypot(x0,y);    
  th = atan(y/x0);
  logq1 = log(q1);
  gr = (x0-0.5)*logq1-th*y-x0+0.5*log(2.0*M_PI);
  for (k=0; k<10; k++){
    t = pow(q1,-1.0-2.0*k);
    gr += (a[k]*t*cos((2.0*k+1.0)*th));
  }
  if (x <= 7.0) {
    gr1 = 0.0;
    for (j=0; j<na; j++) {
      gr1 += (0.5*log((x+j)*(x+j)+y*y));
    }
    gr -= gr1;
  }
  if (x1 < 0.0) {
    q1 = hypot(x,y);    
    sr = -sin(M_PI*x)*cosh(M_PI*y);
    si = -cos(M_PI*x)*sinh(M_PI*y);
    q2 = hypot(sr,si);  
    gr = log(M_PI/(q1*q2))-gr;
  }
  return gr;
} 
