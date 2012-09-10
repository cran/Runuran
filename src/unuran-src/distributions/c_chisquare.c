/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "chisquare";
#define nu  params[0]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_chisquare( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_chisquare( double x, const UNUR_DISTR *distr );
static double _unur_cdf_chisquare( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_gamma
static double _unur_invcdf_chisquare( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_chisquare( UNUR_DISTR *distr );
static int _unur_upd_area_chisquare( UNUR_DISTR *distr );
static int _unur_set_params_chisquare( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_chisquare(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  if (_unur_isfsame(nu,2.))
    return exp(-x/2. - LOGNORMCONSTANT);
  return exp( log(x) * (nu/2. - 1.) - x/2. - LOGNORMCONSTANT );
} 
double
_unur_dpdf_chisquare(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  if (_unur_isfsame(nu,2.))
    return ( -exp(-x/2. - LOGNORMCONSTANT) / 2. );
  return ( exp( log(x) * (nu/2. - 2.) - x/2. - LOGNORMCONSTANT) * (nu - 2. - x)/2. );
} 
double
_unur_cdf_chisquare(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return _unur_SF_incomplete_gamma(x/2.,nu/2.);
} 
#ifdef _unur_SF_invcdf_gamma
double
_unur_invcdf_chisquare( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  return _unur_SF_invcdf_gamma(x, 0.5*nu, 2.);
} 
#endif
int
_unur_upd_mode_chisquare( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.nu >= 2.) ? (DISTR.nu - 2.) : 0.;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_chisquare( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = _unur_SF_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2.);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_chisquare( DISTR.domain[1],distr) 
		 - _unur_cdf_chisquare( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_chisquare( UNUR_DISTR *distr, const double *params, int n_params )
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
    DISTR.domain[1] = UNUR_INFINITY; 
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_chisquare( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_CHISQUARE;
  distr->name = distr_name;
  DISTR.init = NULL;
  DISTR.pdf  = _unur_pdf_chisquare;   
  DISTR.dpdf = _unur_dpdf_chisquare;  
  DISTR.cdf  = _unur_cdf_chisquare;   
#ifdef _unur_SF_invcdf_gamma
  DISTR.invcdf = _unur_invcdf_chisquare;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_chisquare(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_SF_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2.);
  DISTR.mode = (DISTR.nu >= 2.) ? (DISTR.nu - 2.) : 0.;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_chisquare;
  DISTR.upd_mode  = _unur_upd_mode_chisquare; 
  DISTR.upd_area  = _unur_upd_area_chisquare; 
  return distr;
} 
#undef nu
#undef DISTR
