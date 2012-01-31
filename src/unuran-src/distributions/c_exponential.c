/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "exponential";
#define sigma  params[0]
#define theta  params[1]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_cdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_exponential( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_exponential( UNUR_DISTR *distr );
static int _unur_upd_area_exponential( UNUR_DISTR *distr );
static int _unur_set_params_exponential( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / sigma;
  return ( (x<0.) ? 0. : exp(-x - LOGNORMCONSTANT) );
} 
double
_unur_logpdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / sigma;
  return ( (x<0.) ? -INFINITY : (-x - LOGNORMCONSTANT) );
} 
double
_unur_dpdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / sigma;
  return ( (x<0.) ? 0. : -exp(-x - 2.*LOGNORMCONSTANT) );
} 
double
_unur_dlogpdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / sigma;
  return ( (x<0.) ? 0. : -1./sigma );
} 
double
_unur_cdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / sigma;
  return ( (x<0.) ? 0. : 1.-exp(-x) );
} 
double
_unur_invcdf_exponential( double U, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  double X;
  X = - log( 1. - U );
  return ((DISTR.n_params==0) ? X : theta + sigma * X);
} 
int
_unur_upd_mode_exponential( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_exponential( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = log(DISTR.sigma);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_exponential( DISTR.domain[1],distr) 
		 - _unur_cdf_exponential( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_exponential( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params > 0 && sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.sigma = 1.;
  DISTR.theta = 0.;
  switch (n_params) {
  case 2:
    DISTR.theta = theta;
  case 1:
    DISTR.sigma = sigma;
    n_params = 2;           
  default:
    break;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.theta;     
    DISTR.domain[1] = INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_exponential( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_EXPONENTIAL;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_exponential_init;
  DISTR.pdf     = _unur_pdf_exponential;     
  DISTR.logpdf  = _unur_logpdf_exponential;  
  DISTR.dpdf    = _unur_dpdf_exponential;    
  DISTR.dlogpdf = _unur_dlogpdf_exponential; 
  DISTR.cdf     = _unur_cdf_exponential;     
  DISTR.invcdf  = _unur_invcdf_exponential;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_exponential(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = log(DISTR.sigma);
  DISTR.mode = DISTR.theta;   
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_exponential;
  DISTR.upd_mode  = _unur_upd_mode_exponential; 
  DISTR.upd_area  = _unur_upd_area_exponential; 
  return distr;
} 
#undef sigma 
#undef theta 
#undef DISTR
