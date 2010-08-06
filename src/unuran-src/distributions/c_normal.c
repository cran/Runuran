/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "normal";
#define mu    params[0]
#define sigma params[1]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_normal( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_normal( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_normal( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_normal( double x, const UNUR_DISTR *distr );
static double _unur_cdf_normal( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_normal( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_normal( UNUR_DISTR *distr );
static int _unur_upd_area_normal( UNUR_DISTR *distr );
static int _unur_set_params_normal( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_normal( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - mu) / sigma;
  return exp(-x*x/2. + LOGNORMCONSTANT); 
} 
double
_unur_logpdf_normal( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - mu) / sigma;
  return (-x*x/2. + LOGNORMCONSTANT); 
} 
double
_unur_dpdf_normal( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0) {
    x = (x - mu) / sigma;
  }
  return ( -x * exp(-x*x/2. + LOGNORMCONSTANT) / sigma );
} 
double
_unur_dlogpdf_normal( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    return (mu - x)/(sigma*sigma);
  else
    return (-x);
} 
double
_unur_cdf_normal( double x, const UNUR_DISTR *distr ) 
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - mu) / sigma;
  return _unur_SF_cdf_normal(x);
} 
double
_unur_invcdf_normal( double u, const UNUR_DISTR *distr ) 
{
  register const double *params = DISTR.params;
  double X;
  X = _unur_SF_invcdf_normal(u);
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
int
_unur_upd_mode_normal( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.mu;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_normal( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = - log(M_SQRTPI * M_SQRT2 * DISTR.sigma);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_normal( DISTR.domain[1],distr) 
		 - _unur_cdf_normal( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_normal( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params > 1 && sigma <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.mu    = 0.;
  DISTR.sigma = 1.;
  switch (n_params) {
  case 2:
    DISTR.sigma = sigma;
  case 1:
    DISTR.mu = mu;
    n_params = 2;           
  default:
    break;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       
    DISTR.domain[1] = INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_normal( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_NORMAL;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_normal_init;
  DISTR.pdf     = _unur_pdf_normal;     
  DISTR.logpdf  = _unur_logpdf_normal;  
  DISTR.dpdf    = _unur_dpdf_normal;    
  DISTR.dlogpdf = _unur_dlogpdf_normal; 
  DISTR.cdf     = _unur_cdf_normal;     
  DISTR.invcdf  = _unur_invcdf_normal;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_normal(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = - log(M_SQRTPI * M_SQRT2 * DISTR.sigma);
  DISTR.mode = DISTR.mu;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_normal;
  DISTR.upd_mode  = _unur_upd_mode_normal; 
  DISTR.upd_area  = _unur_upd_area_normal; 
  return distr;
} 
#undef mu
#undef sigma
#undef DISTR
