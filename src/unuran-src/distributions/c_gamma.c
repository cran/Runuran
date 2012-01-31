/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "gamma";
#define alpha  params[0]   
#define beta   params[1]   
#define gamma  params[2]   
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_cdf_gamma( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_gamma
static double _unur_invcdf_gamma( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_gamma( UNUR_DISTR *distr );
static int _unur_upd_area_gamma( UNUR_DISTR *distr );
static double _unur_lognormconstant_gamma(const double *params, int n_params);
static int _unur_set_params_gamma( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x-gamma) / beta;
  if (_unur_isone(alpha) && x >= 0.)
    return exp( -x - LOGNORMCONSTANT);
  if (x > 0.)
    return exp( (alpha-1.)*log(x) - x - LOGNORMCONSTANT);
  if (_unur_iszero(x))
    return (alpha>1. ? 0. : INFINITY);
  return 0.;
} 
double
_unur_logpdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x-gamma) / beta;
  if (_unur_isone(alpha) && x >= 0.)
    return ( -x - LOGNORMCONSTANT);
  if (x > 0.)
    return ( (alpha-1.)*log(x) - x - LOGNORMCONSTANT);
  if (_unur_iszero(x))
    return (alpha>1. ? -INFINITY : INFINITY);
  return -INFINITY;
} 
double
_unur_dpdf_gamma( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x-gamma) / beta;
  if (_unur_isone(alpha) && x>=0)
    return( -exp(-x - LOGNORMCONSTANT) / beta );
  if (x > 0.)
    return ( exp( log(x) * (alpha-2.) - x - LOGNORMCONSTANT) *  ((alpha-1.) -x) / beta ); 
  if (_unur_iszero(x) && alpha < 2.)
    return (alpha>1. ? INFINITY : -INFINITY);
  return 0.;
} 
double
_unur_dlogpdf_gamma( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x-gamma) / beta;
  if (_unur_isone(alpha) && x >= 0.)
    return -1./beta;
  if (x > 0.)
    return ((alpha-1.)/x - 1)/beta;
  if (_unur_iszero(x))
    return (alpha>1. ? INFINITY : -INFINITY);
  return 0.;
} 
double
_unur_cdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x-gamma) / beta;
  if (x <= 0.)
    return 0.;
  if (_unur_isinf(x)==1)
    return 1.;
  return _unur_SF_incomplete_gamma(x,alpha);
} 
#ifdef _unur_SF_invcdf_gamma
double
_unur_invcdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  if (DISTR.n_params == 1)
    return _unur_SF_invcdf_gamma(x, alpha, 1.);
  else
    return (gamma + _unur_SF_invcdf_gamma(x, alpha, beta));
} 
#endif
int
_unur_upd_mode_gamma( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  DISTR.mode = (alpha >= 1.) ? (alpha - 1.) : 0.;
  if (DISTR.n_params > 1)
    DISTR.mode = DISTR.mode * beta + gamma;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  if (alpha < 1.) {
    double center = alpha * beta + gamma;
    center = _unur_max(center,DISTR.domain[0]);
    center = _unur_min(center,DISTR.domain[1]);
    unur_distr_cont_set_center(distr,center);
  }
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_gamma( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_gamma( DISTR.domain[1],distr) 
		 - _unur_cdf_gamma( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
double
_unur_lognormconstant_gamma( const double *params, int n_params )
{
  if (n_params > 1)
    return ( _unur_SF_ln_gamma(alpha) + log(beta) );
  else
    return (_unur_SF_ln_gamma(alpha));
} 
int
_unur_set_params_gamma( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params > 1 && beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.alpha = alpha;
  DISTR.beta  = 1.;
  DISTR.gamma = 0.;
  switch (n_params) {
  case 3:
    DISTR.gamma = gamma;
  case 2:
    DISTR.beta = beta;
    n_params = 3;           
  default:
    break;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.gamma;  
    DISTR.domain[1] = INFINITY;     
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_gamma( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_GAMMA;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_gamma_init;
  DISTR.pdf     = _unur_pdf_gamma;     
  DISTR.logpdf  = _unur_logpdf_gamma;  
  DISTR.dpdf    = _unur_dpdf_gamma;    
  DISTR.dlogpdf = _unur_dlogpdf_gamma; 
  DISTR.cdf     = _unur_cdf_gamma;     
#ifdef _unur_SF_invcdf_gamma
  DISTR.invcdf  = _unur_invcdf_gamma;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_gamma(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);
  _unur_upd_mode_gamma( distr );
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_gamma;
  DISTR.upd_mode  = _unur_upd_mode_gamma; 
  DISTR.upd_area  = _unur_upd_area_gamma; 
  return distr;
} 
#undef alpha
#undef beta 
#undef gamma
#undef DISTR
