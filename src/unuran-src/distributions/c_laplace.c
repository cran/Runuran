/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "laplace";
#define theta  params[0]
#define phi    params[1]
#define DISTR distr->data.cont
static double _unur_pdf_laplace( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_laplace( double x, const UNUR_DISTR *distr );
static double _unur_cdf_laplace( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_laplace( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_laplace( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_laplace( UNUR_DISTR *distr );
static int _unur_upd_area_laplace( UNUR_DISTR *distr );
static int _unur_set_params_laplace( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_laplace( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<theta) 
    return ( exp( (x-theta)/phi ) / (2.*phi) ); 
  return ( exp( (theta-x)/phi ) / (2.*phi) ); 
} 
double
_unur_logpdf_laplace( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<theta) 
    return ( (x-theta)/phi - log(2.*phi) ); 
  return (  (theta-x)/phi - log(2.*phi) ); 
} 
double
_unur_dpdf_laplace( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double z;
  z = (x>theta) ? (x-theta)/phi : (theta-x)/phi;
  if (_unur_iszero(z))   
    return 0.;           
  return ( ((x>theta) ? -exp(-z)/phi : exp(-z)/phi) / (2.*phi) );
} 
double
_unur_dlogpdf_laplace( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<theta) ? 1./phi : -1./phi );
} 
double
_unur_cdf_laplace( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double z;
  z = (x-theta)/phi;
  return ( (x>theta) ? 1.-0.5 * exp(-z) : 0.5*exp(z) );
} 
int
_unur_upd_mode_laplace( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_laplace( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_laplace( DISTR.domain[1],distr) 
		 - _unur_cdf_laplace( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_laplace( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params == 2 && phi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"phi <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.theta = 0.;
  DISTR.phi   = 1.;
  switch (n_params) {
  case 2:
    DISTR.phi = phi;
  case 1:
    DISTR.theta = theta;
  default:
    n_params = 2;           
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       
    DISTR.domain[1] = INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_laplace( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_LAPLACE;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_laplace_init;
  DISTR.pdf     = _unur_pdf_laplace;     
  DISTR.logpdf  = _unur_logpdf_laplace;  
  DISTR.dpdf    = _unur_dpdf_laplace;    
  DISTR.dlogpdf = _unur_dlogpdf_laplace; 
  DISTR.cdf     = _unur_cdf_laplace;     
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_laplace(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = DISTR.theta;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_laplace;
  DISTR.upd_mode  = _unur_upd_mode_laplace; 
  DISTR.upd_area  = _unur_upd_area_laplace; 
  return distr;
} 
#undef theta
#undef phi  
#undef DISTR
