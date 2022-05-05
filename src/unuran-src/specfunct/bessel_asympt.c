/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
double
_unur_bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled)
{
  double z;                   
  double sz, t, t2, eta;      
  double d, u1t,u2t,u3t,u4t;  
  double res;                 
  z = x / nu;
  sz = hypot(1,z);   
  t = 1. / sz;
  t2 = t*t;
  eta = (expon_scaled) ? (1./(z + sz)) : sz;
  eta += log(z) - log1p(sz);                  
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125. 
		 + t2 * (-94121676.
			 + t2 * (349922430. 
				 + t2 * (-446185740. 
					 + t2 * 185910725.)))) / 39813120.;
  d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;
  res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);
  return (islog ? res : exp(res));
} 
