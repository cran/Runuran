/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "chi";
#define nu  params[0]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_chi( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_chi( double x, const UNUR_DISTR *distr );
static double _unur_cdf_chi( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_chi( UNUR_DISTR *distr );
static int _unur_upd_area_chi( UNUR_DISTR *distr );
static int _unur_set_params_chi( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_chi(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return (exp( log(x)*(nu - 1.) -x*x/2. - LOGNORMCONSTANT ));
} 
double
_unur_dpdf_chi(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return ( exp(log(x)*(nu - 2.)-x*x/2. - LOGNORMCONSTANT) * (nu - 1. - x*x) );
} 
double
_unur_cdf_chi(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return _unur_sf_incomplete_gamma(x*x/2.,nu/2.);
} 
int
_unur_upd_mode_chi( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.nu >= 1.) ? sqrt(DISTR.nu - 1.) : 0.;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_chi( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = _unur_sf_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2. - 1.);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_chi( DISTR.domain[1],distr) 
		 - _unur_cdf_chi( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_chi( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.nu = nu;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          
    DISTR.domain[1] = INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_chi( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_CHI;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_chi_init;
  DISTR.pdf  = _unur_pdf_chi;   
  DISTR.dpdf = _unur_dpdf_chi;  
  DISTR.cdf  = _unur_cdf_chi;   
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_chi(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_sf_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2. - 1.);
  DISTR.mode = (DISTR.nu >= 1.) ? sqrt(DISTR.nu - 1.) : 0.;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_chi;
  DISTR.upd_mode  = _unur_upd_mode_chi; 
  DISTR.upd_area  = _unur_upd_area_chi; 
  return distr;
} 
#undef nu
#undef DISTR
