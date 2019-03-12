/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "powerexponential";
#define tau  params[0]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_powerexponential( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_powerexponential( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_powerexponential( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_powerexponential( double x, const UNUR_DISTR *distr );
static double _unur_cdf_powerexponential( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_powerexponential( UNUR_DISTR *distr );
static int _unur_upd_area_powerexponential( UNUR_DISTR *distr );
static int _unur_set_params_powerexponential( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_powerexponential( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return exp( - pow( fabs(x), tau ) - LOGNORMCONSTANT);
} 
double
_unur_logpdf_powerexponential( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( - pow( fabs(x), tau ) - LOGNORMCONSTANT);
} 
double
_unur_dpdf_powerexponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  register double tmp;
  if (_unur_iszero(x))    
    return 0.;            
  tmp = exp( -pow(fabs(x),tau) - LOGNORMCONSTANT + (tau-1.)*log(fabs(x)) ) * tau;
  return ( (x<0.) ? tmp : -tmp );
} 
double
_unur_dlogpdf_powerexponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (_unur_iszero(x))    
    return 0.;            
  return (x<0. ? 1. : -1.) * (tau-1.)* pow(fabs(x), tau-1.);
} 
double
_unur_cdf_powerexponential( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double cdf;
  cdf = _unur_SF_incomplete_gamma(pow(fabs(x),tau),1./tau) / 2.;
  return ((x<0.) ? 0.5 - cdf : 0.5 + cdf);
} 
int
_unur_upd_mode_powerexponential( UNUR_DISTR *distr )
{
  DISTR.mode = 0;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_powerexponential( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = _unur_SF_ln_gamma(1. + 1./DISTR.tau) + M_LN2;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_powerexponential( DISTR.domain[1],distr) 
		 - _unur_cdf_powerexponential( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_powerexponential( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (tau <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"tau <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.tau = tau;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -UNUR_INFINITY;       
    DISTR.domain[1] = UNUR_INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_powerexponential( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_POWEREXPONENTIAL;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_powerexponential_init;
  DISTR.pdf     = _unur_pdf_powerexponential;     
  DISTR.logpdf  = _unur_logpdf_powerexponential;  
  DISTR.dpdf    = _unur_dpdf_powerexponential;    
  DISTR.dlogpdf = _unur_dlogpdf_powerexponential; 
  DISTR.cdf     = _unur_cdf_powerexponential;     
  DISTR.cdf     = _unur_cdf_powerexponential;     
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_powerexponential(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_SF_ln_gamma(1. + 1./DISTR.tau) + M_LN2;
  DISTR.mode = 0;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_powerexponential;
  DISTR.upd_mode  = _unur_upd_mode_powerexponential; 
  DISTR.upd_area  = _unur_upd_area_powerexponential; 
  return distr;
} 
#undef tau
#undef DISTR
