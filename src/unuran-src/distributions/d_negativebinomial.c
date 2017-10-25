/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "negativebinomial";
#define p  params[0]
#define r  params[1]
#define DISTR distr->data.discr
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)
static double _unur_pmf_negativebinomial( int k, const UNUR_DISTR *distr );
#ifdef _unur_SF_cdf_negativebinomial
static double _unur_cdf_negativebinomial( int k, const UNUR_DISTR *distr ); 
#endif
#ifdef _unur_SF_invcdf_negativebinomial
static int    _unur_invcdf_negativebinomial( double u, const UNUR_DISTR *distr ); 
#endif
static int _unur_upd_mode_negativebinomial( UNUR_DISTR *distr );
static int _unur_upd_sum_negativebinomial( UNUR_DISTR *distr );
static int _unur_set_params_negativebinomial( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_negativebinomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  if (k<0)
    return 0.;
  else
    return exp( k*log(1-p) 
		+ _unur_SF_ln_gamma(k+r) - _unur_SF_ln_gamma(k+1.) - LOGNORMCONSTANT ) ;
} 
#ifdef _unur_SF_cdf_negativebinomial
double
_unur_cdf_negativebinomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  if (k<0)
    return 0.;
  else
    return _unur_SF_cdf_negativebinomial(k,r,p);
} 
#endif
#ifdef _unur_SF_invcdf_negativebinomial
int
_unur_invcdf_negativebinomial(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;
  x = _unur_SF_invcdf_negativebinomial(u,r,p);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} 
#endif
int
_unur_upd_mode_negativebinomial( UNUR_DISTR *distr )
{
  if (DISTR.r > 1.) {
    DISTR.mode = (int) ((1.+UNUR_EPSILON) * (DISTR.r - 1.) * (1. - DISTR.p) / DISTR.p);
  }
  else { 
    DISTR.mode = 0;
  }
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_sum_negativebinomial( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = - DISTR.r * log(DISTR.p) + _unur_SF_ln_gamma(DISTR.r);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
#ifdef _unur_SF_cdf_negativebinomial
  DISTR.sum = ( _unur_cdf_negativebinomial( DISTR.domain[1],distr) 
		 - _unur_cdf_negativebinomial( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif
} 
int
_unur_set_params_negativebinomial( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (p <= 0. || p >= 1. || r <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1 || r <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.p = p;
  DISTR.r = r;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;           
    DISTR.domain[1] = INT_MAX;     
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_negativebinomial( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_NEGATIVEBINOMIAL;
  distr->name = distr_name;
  DISTR.pmf  = _unur_pmf_negativebinomial;   
#ifdef _unur_SF_cdf_negativebinomial
  DISTR.cdf  = _unur_cdf_negativebinomial;   
#endif
#ifdef _unur_SF_invcdf_negativebinomial
  DISTR.invcdf = _unur_invcdf_negativebinomial;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_negativebinomial(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  _unur_upd_sum_negativebinomial( distr );
  _unur_upd_mode_negativebinomial(distr);
  DISTR.sum = 1.;
  DISTR.set_params = _unur_set_params_negativebinomial;
  DISTR.upd_mode = _unur_upd_mode_negativebinomial; 
  DISTR.upd_sum  = _unur_upd_sum_negativebinomial;  
  return distr;
} 
#undef p
#undef r
#undef DISTR
