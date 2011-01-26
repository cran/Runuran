/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "logarithmic";
#define theta  params[0]
#define DISTR distr->data.discr
#define NORMCONSTANT (distr->data.discr.norm_constant)
#undef  HAVE_CDF
static double _unur_pmf_logarithmic( int k, const UNUR_DISTR *distr );
#ifdef HAVE_CDF
static double _unur_cdf_logarithmic( int k, const UNUR_DISTR *distr );      
#endif
static int _unur_upd_mode_logarithmic( UNUR_DISTR *distr );
static int _unur_upd_sum_logarithmic( UNUR_DISTR *distr );
static int _unur_set_params_logarithmic( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_logarithmic(int k, const UNUR_DISTR *distr)
{ 
  return ((k<1) ? 0. : pow( DISTR.theta, (double)k ) / k * NORMCONSTANT);
} 
#ifdef HAVE_CDF
double
_unur_cdf_logarithmic(int k, const UNUR_DISTR *distr)
{ 
  return 0.;
} 
#endif
int
_unur_upd_mode_logarithmic( UNUR_DISTR *distr )
{
  DISTR.mode = 1;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_sum_logarithmic( UNUR_DISTR *distr )
{
  NORMCONSTANT = -1. / log( 1.-DISTR.theta);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
#ifdef HAVE_CDF
  DISTR.sum = ( _unur_cdf_logarithmic( DISTR.domain[1],distr) 
		 - _unur_cdf_logarithmic( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif
} 
int
_unur_set_params_logarithmic( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (theta <= 0. || theta >= 1.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0 || theta >= 1");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.theta = theta;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 1;           
    DISTR.domain[1] = INT_MAX;     
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_logarithmic( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_LOGARITHMIC;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_logarithmic_init;
  DISTR.pmf  = _unur_pmf_logarithmic;   
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_logarithmic;   
#endif           
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE | 
		 UNUR_DISTR_SET_PMFSUM );
  if (_unur_set_params_logarithmic(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = -1. / log( 1.-DISTR.theta);
  DISTR.mode = 1;
  DISTR.sum = 1.;
  DISTR.set_params = _unur_set_params_logarithmic;
  DISTR.upd_mode = _unur_upd_mode_logarithmic; 
  DISTR.upd_sum  = _unur_upd_sum_logarithmic; 
  return distr;
} 
#undef theta
#undef DISTR
