/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "poisson";
#define theta  params[0]
#define DISTR distr->data.discr
static double _unur_pmf_poisson( int k, const UNUR_DISTR *distr );
static double _unur_cdf_poisson( int k, const UNUR_DISTR *distr );      
#ifdef _unur_SF_invcdf_binomial
static int    _unur_invcdf_poisson( double u, const UNUR_DISTR *distr ); 
#endif
static int _unur_upd_mode_poisson( UNUR_DISTR *distr );
static int _unur_upd_sum_poisson( UNUR_DISTR *distr );
static int _unur_set_params_poisson( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_poisson(int k, const UNUR_DISTR *distr)
{ 
  if (k>=0)
    return exp( -DISTR.theta + k * log(DISTR.theta) - _unur_SF_ln_factorial(k) );
  else
    return 0.;
} 
double
_unur_cdf_poisson(int k, const UNUR_DISTR *distr)
{ 
  if (k>=0)
    return (1.-_unur_SF_incomplete_gamma(DISTR.theta,k+1.));
  else
    return 0.;
} 
#ifdef _unur_SF_invcdf_poisson
int
_unur_invcdf_poisson(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;
  x = _unur_SF_invcdf_poisson(u,theta);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} 
#endif
int
_unur_upd_mode_poisson( UNUR_DISTR *distr )
{
  DISTR.mode = (int) DISTR.theta;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_sum_poisson( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.sum = ( _unur_cdf_poisson( DISTR.domain[1],distr) 
		 - _unur_cdf_poisson( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_poisson( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.theta = theta;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;           
    DISTR.domain[1] = INT_MAX;     
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_poisson( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_POISSON;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_poisson_init;
  DISTR.pmf  = _unur_pmf_poisson;   
  DISTR.cdf  = _unur_cdf_poisson;   
#ifdef _unur_SF_invcdf_poisson
  DISTR.invcdf = _unur_invcdf_poisson;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_poisson(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = (int) DISTR.theta;
  DISTR.sum = 1.;
  DISTR.set_params = _unur_set_params_poisson;
  DISTR.upd_mode = _unur_upd_mode_poisson; 
  DISTR.upd_sum  = _unur_upd_sum_poisson;  
  return distr;
} 
#undef theta
#undef DISTR
