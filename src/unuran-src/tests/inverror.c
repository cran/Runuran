/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <distr/distr_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include <methods/hinv.h>
#include <methods/ninv.h>
#include <methods/pinv.h>
#include "unuran_tests.h"
static char test_name[] = "InvError";
static double qrng (int i, int samplesize);
double
unur_test_u_error( const UNUR_GEN *gen, 
		   double *max_error, double *MAE, double threshold,
		   int samplesize, int randomized, int testtails, 
		   int verbosity, FILE *out )
{
#define DISTR   gen->distr->data.cont
  double CDFmin, CDFmax;     
  double (*quantile)(const UNUR_GEN *, double);  
  double U, X;               
  double cdfX;               
  double uerror, umax, usum; 
  double penalty = 0;        
  int j;                     
  _unur_check_NULL(test_name,gen,-1.);
  if (verbosity) { _unur_check_NULL(test_name,out,-1.); }
  if (samplesize < 1000) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"samplesize too small --> increased to 1000");
    samplesize = 1000;
  }
  switch (gen->method) {
  case UNUR_METH_HINV:
    quantile = unur_hinv_eval_approxinvcdf;
    break;
  case UNUR_METH_NINV:
    quantile = unur_ninv_eval_approxinvcdf;
    break;
  case UNUR_METH_PINV:
    quantile = unur_pinv_eval_approxinvcdf;
    break;
  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
    return -1.;
  }
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required");
    return -2.;
  }
  CDFmin = (DISTR.trunc[0] > -INFINITY) ? _unur_cont_CDF((DISTR.trunc[0]),(gen->distr)) : 0.;
  CDFmax = (DISTR.trunc[1] < INFINITY)  ? _unur_cont_CDF((DISTR.trunc[1]),(gen->distr)) : 1.;
  umax = 0.;
  usum = 0.;
  for(j=0;j<samplesize;j++) {
    if (randomized)
      U = _unur_call_urng(gen->urng);
    else
      U = (testtails) ? qrng(j,samplesize) : (j+0.5) / ((double) samplesize);
    X = quantile(gen,U);
    cdfX = _unur_cont_CDF(X,gen->distr);
    uerror = fabs( U*(CDFmax - CDFmin) - (cdfX-CDFmin));
    usum += uerror;
    if (uerror > umax) {
      umax = uerror;
    }
    if (_unur_FP_less(threshold,uerror)) {
      penalty += 1. + 10.*(uerror - threshold) / threshold;
      if (verbosity)
	fprintf(out,"\tmax u-error exceeded at %g: %g (>%g)\n",
		X,uerror,threshold);
    }
  }
  *max_error = umax;
  *MAE = usum/samplesize;
  return penalty/samplesize;
#undef DISTR
} 
double 
qrng (int i, int samplesize)
{
  double U;
  int tail;
  tail = (int) (0.05 * samplesize);
  i = i % samplesize;
  if (i < tail) {
    U = (i+0.5) / (1.e5 * tail); 
  }
  else if (i >= samplesize-tail) {
    i -= samplesize-tail;
    U = 1. - (i+0.5) / (1.e5 * tail);
  }
  else {
    i -= tail;
    U = (i+0.5) / (samplesize - 2.*tail);
  }
  return U;
}  
