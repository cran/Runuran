/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <distr/distr_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <methods/dgt.h>
#include <methods/dstd.h>
#include <methods/dstd_struct.h>
#include <methods/hinv.h>
#include <methods/mixt.h>
#include <methods/mixt_struct.h>
#include <methods/ninv.h>
#include <methods/pinv.h>
#include "unuran_tests.h"
static char test_name[] = "InvError";
static double uerror_cont ( const UNUR_GEN *gen, 
			    double *max_error, double *MAE, double threshold,
			    int samplesize, int randomized, int testtails, 
			    int verbosity, FILE *out );
static double uerror_discr( const UNUR_GEN *gen, 
			    double *max_error, double *MAE, double threshold,
			    int samplesize, int randomized, int testtails, 
			    int verbosity, FILE *out );
static double qrng (int i, int samplesize);
double
unur_test_u_error( const UNUR_GEN *gen, 
		   double *max_error, double *MAE, double threshold,
		   int samplesize, int randomized, int testtails, 
		   int verbosity, FILE *out )
{
  _unur_check_NULL(test_name,gen,-1.);
  if (verbosity) { _unur_check_NULL(test_name,out,-1.); }
  if (samplesize < 1000) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"samplesize too small --> increased to 1000");
    samplesize = 1000;
  }
  if ( (gen->method == UNUR_METH_HINV) ||
       (gen->method == UNUR_METH_NINV) ||
       (gen->method == UNUR_METH_PINV) ||
       (gen->method == UNUR_METH_CSTD && ((struct unur_cstd_gen*)gen->datap)->is_inversion) ||
       (gen->method == UNUR_METH_MIXT && ((struct unur_mixt_gen*)gen->datap)->is_inversion)
       ) {
    return uerror_cont(gen,max_error,MAE,threshold,samplesize,
		       randomized,testtails,verbosity,out);
  }
  if ( (gen->method == UNUR_METH_DGT) ||
       (gen->method == UNUR_METH_DSTD && ((struct unur_dstd_gen*)gen->datap)->is_inversion)
       ) {
    return uerror_discr(gen,max_error,MAE,threshold,samplesize,
			randomized,testtails,verbosity,out);
  }
  _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
  return -1.;
} 
double
uerror_cont( const UNUR_GEN *gen, 
	     double *max_error, double *MAE, double threshold,
	     int samplesize, int randomized, int testtails, 
	     int verbosity, FILE *out )
{
#define DISTR   gen->distr->data.cont
  double CDFmin, CDFmax;     
  double (*quantile)(const UNUR_GEN *, double);  
  double U;                  
  double X;                  
  double cdfX;               
  double uerror, umax, usum; 
  double penalty = 0;        
  int j;                     
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
  case UNUR_METH_CSTD:
    if (! (((struct unur_cstd_gen*)gen->datap)->is_inversion))
      return -1.;
    quantile = unur_cstd_eval_invcdf;
    break;
  case UNUR_METH_MIXT:
    if (! (((struct unur_mixt_gen*)gen->datap)->is_inversion))
      return -1.;
    quantile = unur_cstd_eval_invcdf;
    break;
  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
    return -1.;
  }
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required");
    return -2.;
  }
  CDFmin = ((DISTR.trunc[0] > -UNUR_INFINITY) 
	    ? _unur_cont_CDF((DISTR.trunc[0]),(gen->distr)) : 0.);
  CDFmax = ((DISTR.trunc[1] < UNUR_INFINITY)
	    ? _unur_cont_CDF((DISTR.trunc[1]),(gen->distr)) : 1.);
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
uerror_discr( const UNUR_GEN *gen, 
		double *max_error, double *MAE, double threshold,
		int samplesize, int randomized,
		int testtails ATTRIBUTE__UNUSED, 
		int verbosity, FILE *out )
{
#define DISTR   gen->distr->data.discr
  int (*quantile)(const UNUR_GEN *, double);  
  double U;                  
  int K;                     
  double cdfK;               
  double uerror, umax, usum; 
  double penalty = 0;        
  int j;                     
  switch (gen->method) {
  case UNUR_METH_DGT:
    quantile = unur_dgt_eval_invcdf;
    break;
  case UNUR_METH_DSTD:
    if (! (((struct unur_dstd_gen*)gen->datap)->is_inversion))
      return -1.;
    quantile = unur_dstd_eval_invcdf;
    break;
  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
    return -1.;
  }
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required");
    return -2.;
  }
  umax = 0.;
  usum = 0.;
  for(j=0;j<samplesize;j++) {
    if (randomized)
      U = _unur_call_urng(gen->urng);
    else
      U = (j+0.5) / ((double) samplesize);
    K = (int) quantile(gen,U);
    uerror = 0.;
    cdfK = _unur_discr_CDF(K,gen->distr);
    if (cdfK < U) {
      uerror = U - cdfK;
    }
    else {
      cdfK = _unur_discr_CDF(K-1,gen->distr);
      uerror = cdfK - U;
      uerror = _unur_max(0.,uerror);
    }
    usum += uerror;
    if (uerror > umax)
      umax = uerror;
    if (uerror > umax) {
      umax = uerror;
    }
    if (_unur_FP_less(threshold,uerror)) {
      penalty += 1. + 10.*(uerror - threshold) / threshold;
      if (verbosity)
	fprintf(out,"\tmax u-error exceeded at U=%g: %g (>%g)\n",
		U,uerror,threshold);
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
