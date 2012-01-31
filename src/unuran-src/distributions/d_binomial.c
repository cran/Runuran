/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "binomial";
#define n  params[0]
#define p  params[1]
#define DISTR distr->data.discr
static double _unur_pmf_binomial( int k, const UNUR_DISTR *distr );
static double _unur_cdf_binomial( int k, const UNUR_DISTR *distr ); 
#ifdef _unur_SF_invcdf_binomial
static int    _unur_invcdf_binomial( double u, const UNUR_DISTR *distr ); 
#endif
static int _unur_upd_mode_binomial( UNUR_DISTR *distr );
static int _unur_upd_sum_binomial( UNUR_DISTR *distr );
static int _unur_set_params_binomial( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_binomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  if ( k<0 || k>(n+0.5) )
    return 0.;
  else
    return exp( k * log(p) + (n-k) * log(1.-p) +
		_unur_SF_ln_factorial(n) - _unur_SF_ln_factorial(k) - _unur_SF_ln_factorial(n-k) ) ;
} 
double
_unur_cdf_binomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  if (k<0)
    return 0.;
  if (k==0)
    return exp(n*(log(1.-p)));
  if(k>(n-0.5))
    return 1.;
  return(_unur_SF_incomplete_beta(1.-p, n-k, k+1.));
} 
#ifdef _unur_SF_invcdf_binomial
int
_unur_invcdf_binomial(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;
  x = _unur_SF_invcdf_binomial(u,n,p);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} 
#endif
int
_unur_upd_mode_binomial( UNUR_DISTR *distr )
{
  DISTR.mode = (int) ((DISTR.n + 1) * DISTR.p);
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_sum_binomial( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.sum = ( _unur_cdf_binomial( DISTR.domain[1],distr) 
		 - _unur_cdf_binomial( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_binomial( UNUR_DISTR *distr, const double *params, int n_params )
{
  int nh;
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (p <= 0. || p >= 1. || n <= 0.) { 
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1 || n <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.n = nh;
  DISTR.p = p;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;        
    DISTR.domain[1] = nh;       
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_binomial( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_BINOMIAL;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_binomial_init;
  DISTR.pmf  = _unur_pmf_binomial;   
  DISTR.cdf  = _unur_cdf_binomial;   
#ifdef _unur_SF_invcdf_binomial
  DISTR.invcdf = _unur_invcdf_binomial;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_binomial(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = (int) ((DISTR.n + 1) * DISTR.p);
  DISTR.sum = 1.;
  DISTR.set_params = _unur_set_params_binomial;
  DISTR.upd_mode = _unur_upd_mode_binomial; 
  DISTR.upd_sum  = _unur_upd_sum_binomial;  
  return distr;
} 
#undef p
#undef r
#undef DISTR
