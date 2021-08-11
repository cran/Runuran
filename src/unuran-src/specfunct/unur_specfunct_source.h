/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_SPECFUNCT_SOURCE_H_SEEN
#define UNUR_SPECFUNCT_SOURCE_H_SEEN
#ifdef HAVE_LIBRMATH
#  ifdef R_UNURAN
#  else
#    define MATHLIB_STANDALONE
#  endif
#  include <Rmath.h>
#ifdef trunc
#undef trunc
#endif
#ifdef beta
#undef beta
#endif
#define _unur_SF_incomplete_beta(x,a,b)   pbeta((x),(a),(b),TRUE,FALSE)
#define _unur_SF_ln_gamma(x)              lgammafn(x)
#define _unur_SF_ln_factorial(x)          lgammafn((x)+1.)
#define _unur_SF_incomplete_gamma(x,a)    pgamma(x,a,1.,TRUE,FALSE)
#define _unur_SF_bessel_k(x,nu)           bessel_k((x),(nu),1)
#define _unur_SF_ln_bessel_k(x,nu)        (log(bessel_k((x),(nu),2)) - (x))
#define _unur_SF_cdf_normal(x)            pnorm((x),0.,1.,TRUE,FALSE)
#define _unur_SF_invcdf_normal(u)         qnorm((u),0.,1.,TRUE,FALSE)
#define _unur_SF_invcdf_beta(u,p,q)       qbeta((u),(p),(q),TRUE,FALSE)
#define _unur_SF_cdf_F(x,nua,nub)         pf((x),(nua),(nub),TRUE,FALSE)
#define _unur_SF_invcdf_F(u,nua,nub)      qf((u),(nua),(nub),TRUE,FALSE)
#define _unur_SF_invcdf_gamma(u,shape,scale)  qgamma((u),(shape),(scale),TRUE,FALSE)
#define _unur_SF_cdf_student(x,nu)         pt((x),(nu),TRUE,FALSE)
#define _unur_SF_invcdf_student(u,nu)      qt((u),(nu),TRUE,FALSE)
#define _unur_SF_invcdf_binomial(u,n,p)   qbinom((u),(n),(p),TRUE,FALSE)
#define _unur_SF_cdf_hypergeometric(x,N,M,n)  phyper((x),(M),(N)-(M),(n),TRUE,FALSE)
#define _unur_SF_invcdf_hypergeometric(u,N,M,n)  qhyper((u),(M),(N)-(M),(n),TRUE,FALSE)
#define _unur_SF_cdf_negativebinomial(x,n,p)      pnbinom((x),(n),(p),TRUE,FALSE)
#define _unur_SF_invcdf_negativebinomial(u,n,p)   qnbinom((u),(n),(p),TRUE,FALSE)
#define _unur_SF_invcdf_poisson(u,theta)   qpois((u),(theta),TRUE,FALSE)
#else
#define COMPILE_CEPHES
double _unur_cephes_incbet(double a, double b, double x);
#define _unur_SF_incomplete_beta(x,a,b)   _unur_cephes_incbet((a),(b),(x))
double _unur_cephes_lgam(double x);
#define _unur_SF_ln_gamma(x)              _unur_cephes_lgam(x)
#define _unur_SF_ln_factorial(x)          _unur_cephes_lgam((x)+1.)
double _unur_cephes_igam(double a, double x);
#define _unur_SF_incomplete_gamma(x,a)    _unur_cephes_igam((a),(x))
double _unur_cephes_ndtr(double x);
#define _unur_SF_cdf_normal(x)            _unur_cephes_ndtr(x)
double _unur_cephes_ndtri(double x);
#define _unur_SF_invcdf_normal(x)         _unur_cephes_ndtri(x)
#endif
double _unur_bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled);
#define _unur_SF_bessel_k_nuasympt(x,nu,islog,exponscaled) \
  _unur_bessel_k_nuasympt((x),(nu),(islog),(exponscaled))
double _unur_Relcgamma (double x, double y);
#define _unur_SF_Relcgamma(x,y)  _unur_Relcgamma((x),(y))
#if !HAVE_DECL_LOG1P
double _unur_log1p(double x);
#define log1p _unur_log1p
#endif
#if !HAVE_DECL_HYPOT
double _unur_hypot(const double x, const double y);
#define hypot _unur_hypot
#endif
#endif  
