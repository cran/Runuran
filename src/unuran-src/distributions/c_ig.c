/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "ig";
#define mu    params[0]
#define lambda params[1]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_cdf_ig( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_ig( UNUR_DISTR *distr );
static int _unur_upd_area_ig( UNUR_DISTR *distr );
static int _unur_set_params_ig( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_ig( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<=0.)
    return 0.;
  else
    return
      ( sqrt(lambda/(2*M_PI*x*x*x)) * exp( -lambda*(x-mu)*(x-mu) / (2*mu*mu*x) ) );
} 
double
_unur_logpdf_ig( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<0.)
    return -UNUR_INFINITY;
  else
    return ( 0.5* log ( lambda/(2*M_PI*x*x*x) )
	     -lambda*(x-mu)*(x-mu) / (2*mu*mu*x) );
} 
double
_unur_dpdf_ig( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  double res;
  if (x<=0.)
    return 0.;
  res = -exp( -lambda*(x-mu)*(x-mu) / (2*mu*mu*x) );
  res *= sqrt(lambda / (x*x*x));
  res *= (3*mu*mu*x + lambda*(-mu*mu + x*x));
  res /= 2*mu*mu*sqrt(2*M_PI)*x*x;
  return res;
} 
double
_unur_dlogpdf_ig( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (x<=0.)
    return 0.;
  else 
    return 
      (0.5 * lambda*(-1./(mu*mu) + 1./(x*x)) - 3./x);
} 
double
_unur_cdf_ig( double x, const UNUR_DISTR *distr ) 
{
  register const double *params = DISTR.params;
#define Phi(x)   (_unur_SF_cdf_normal(x))
  if (x<=0.)
    return 0.;
  return 
    ( Phi(sqrt(lambda/x)*(x/mu-1.)) 
      + exp(2*lambda/mu) * Phi(-sqrt(lambda/x)*(x/mu+1.)) );
#undef Phi
} 
int
_unur_upd_mode_ig( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  DISTR.mode = 
    (-3.*mu*mu + mu*sqrt(4.*lambda*lambda + 9*mu*mu))/(2*lambda);
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_ig( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = 0.;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_ig( DISTR.domain[1],distr) 
		 - _unur_cdf_ig( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_ig( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (mu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"mu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.mu     = mu;
  DISTR.lambda = lambda;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;             
    DISTR.domain[1] = UNUR_INFINITY;  
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_ig( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_IG;
  distr->name = distr_name;
  DISTR.init = NULL;         
  DISTR.pdf     = _unur_pdf_ig;     
  DISTR.logpdf  = _unur_logpdf_ig;  
  DISTR.dpdf    = _unur_dpdf_ig;    
  DISTR.dlogpdf = _unur_dlogpdf_ig; 
  DISTR.cdf     = _unur_cdf_ig;     
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_ig(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = 0;
  _unur_upd_mode_ig(distr);
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_ig;
  DISTR.upd_mode  = _unur_upd_mode_ig; 
  DISTR.upd_area  = _unur_upd_area_ig; 
  return distr;
} 
#undef mu
#undef lambda
#undef DISTR
