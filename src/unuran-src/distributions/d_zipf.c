/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "zipf";
#define rho  params[0]
#define tau  params[1]
#define DISTR distr->data.discr
#undef HAVE_CDF
#undef HAVE_SUM
static double _unur_pmf_zipf( int k, const UNUR_DISTR *distr );
#ifdef HAVE_CDF
static double _unur_cdf_zipf( int k, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_zipf( UNUR_DISTR *distr );
#ifdef HAVE_SUM
static int _unur_upd_sum_zipf( UNUR_DISTR *distr );
#endif
static int _unur_set_params_zipf( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_zipf(int k, const UNUR_DISTR *distr)
{ 
  return ((k<1) ? 0. : exp( log(k + DISTR.tau) * (-DISTR.rho - 1.) ) );
} 
#ifdef HAVE_CDF
double
_unur_cdf_zipf(int k, const UNUR_DISTR *distr)
{ 
  return ((k<1) ? 0. : 1.);
} 
#endif
int
_unur_upd_mode_zipf( UNUR_DISTR *distr )
{
  DISTR.mode = 1;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
#ifdef HAVE_SUM
int
_unur_upd_sum_zipf( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
#ifdef HAVE_CDF
  DISTR.sum = ( _unur_cdf_zipf( DISTR.domain[1],distr) 
		 - _unur_cdf_zipf( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif
} 
#endif
int
_unur_set_params_zipf( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (rho <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"rho <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params > 1 && tau < 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"tau < 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.rho = rho;
  DISTR.tau = 0.;
  switch (n_params) {
  case 2:
    DISTR.tau = tau;
  default:
    n_params = 2;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 1;           
    DISTR.domain[1] = INT_MAX;     
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_zipf( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_ZIPF;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_zipf_init;
  DISTR.pmf  = _unur_pmf_zipf;   
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_zipf;   
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_SUM
		 UNUR_DISTR_SET_PMFSUM |
#endif
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_zipf(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = 1;
  DISTR.sum  = 1.;
  DISTR.set_params = _unur_set_params_zipf;
  DISTR.upd_mode = _unur_upd_mode_zipf; 
#ifdef HAVE_SUM
  DISTR.upd_sum  = _unur_upd_sum_zipf;  
#endif
  return distr;
} 
#undef rho
#undef tau
#undef DISTR
