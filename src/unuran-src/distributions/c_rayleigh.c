/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] =  "rayleigh";
#define sigma  params[0]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_rayleigh( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_rayleigh( double x, const UNUR_DISTR *distr );
static double _unur_cdf_rayleigh( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_rayleigh( UNUR_DISTR *distr );
static int _unur_upd_area_rayleigh( UNUR_DISTR *distr );
static int _unur_set_params_rayleigh( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_rayleigh( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<=0.) 
    return 0.;
  return (x * exp(-x*x/(2.*sigma*sigma) - LOGNORMCONSTANT ) ); 
} 
double
_unur_dpdf_rayleigh( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double z;
  z = x*x/(sigma*sigma);
  return ( (x<=0.) ? 0. : exp(-z/2 - LOGNORMCONSTANT) * (1-z) ); 
} 
double
_unur_cdf_rayleigh( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<=0.) ? 0. : 1. - exp(-x*x/(2.*sigma*sigma)) );
} 
int
_unur_upd_mode_rayleigh( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.sigma;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_rayleigh( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT =   2. * log(DISTR.sigma);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_rayleigh( DISTR.domain[1],distr) 
		 - _unur_cdf_rayleigh( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_rayleigh( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.sigma = sigma;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;              
    DISTR.domain[1] = INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_rayleigh( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_RAYLEIGH;
  distr->name = distr_name;
  DISTR.init = NULL;    
  DISTR.pdf  = _unur_pdf_rayleigh;  
  DISTR.dpdf = _unur_dpdf_rayleigh; 
  DISTR.cdf  = _unur_cdf_rayleigh;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_rayleigh(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT =   2. * log(DISTR.sigma);
  DISTR.mode = DISTR.sigma;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_rayleigh;
  DISTR.upd_mode  = _unur_upd_mode_rayleigh; 
  DISTR.upd_area  = _unur_upd_area_rayleigh; 
  return distr; 
} 
#undef sigma
#undef DISTR
