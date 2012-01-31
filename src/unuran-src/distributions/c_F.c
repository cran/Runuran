/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "F";
#define nua  params[0]
#define nub  params[1]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_F( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_F( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_F( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_F( double x, const UNUR_DISTR *distr );
static double _unur_cdf_F( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_F
static double _unur_invcdf_F( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_F( UNUR_DISTR *distr );
static int _unur_upd_area_F( UNUR_DISTR *distr );
inline static double _unur_lognormconstant_F( const double *params, int n_params );
static int _unur_set_params_F( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x < 0.)
    return 0.;
  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return INFINITY;
    else if (_unur_isfsame(nua,2.))
      return exp(-LOGNORMCONSTANT);
    else
      return 0.;
  }
  else 
    return exp ((nua/2. - 1.)*log(x) - 0.5*(nua + nub)*log(1. + x * nua / nub) - LOGNORMCONSTANT);
} 
double
_unur_logpdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x < 0.)
    return -INFINITY;
  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return INFINITY;
    else if (_unur_isfsame(nub,2.))
      return -LOGNORMCONSTANT;
    else
      return -INFINITY;
  }
  else 
    return ((nua/2. - 1.)*log(x) - 0.5*(nua + nub)*log(1. + x * nua / nub) - LOGNORMCONSTANT);
} 
double
_unur_dpdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x < 0.)
    return 0.;
  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return -INFINITY;
    else if (_unur_isfsame(nub,2.))
      return -(2.+nub)/nub * exp(-LOGNORMCONSTANT);
    else
      return 0.;
  }
  else
    return _unur_pdf_F(x,distr) * _unur_dlogpdf_F(x,distr);
} 
double
_unur_dlogpdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x < 0.)
    return 0.;
  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return -INFINITY;
    else if (_unur_isfsame(nub,2.))
      return -(2.+nub)/nub;
    else
      return INFINITY;
  }
  else
    return ((nua/2.-1.)/x - nua*(nua+nub)/(2.*nub)/(1.+x*nua/nub));
} 
double
_unur_cdf_F(double x, const UNUR_DISTR *distr)
{ 
#ifdef _unur_SF_cdf_F
  return _unur_SF_cdf_F(x,DISTR.nua,DISTR.nub);
#else
  const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  if (nua * x > nub)
    return 1. - _unur_SF_incomplete_beta(nub / (nub + nua * x), nub/2., nua/2.);
  else
    return _unur_SF_incomplete_beta(nua * x / (nub + nua * x), nua/2., nub/2.);
#endif
} 
#ifdef _unur_SF_invcdf_F
double
_unur_invcdf_F(double x, const UNUR_DISTR *distr)
{
  return _unur_SF_invcdf_F(x,DISTR.nua,DISTR.nub);
} 
#endif
int
_unur_upd_mode_F( UNUR_DISTR *distr )
{
  if (DISTR.nua >= 2.)
    DISTR.mode = ((DISTR.nua - 2.) * DISTR.nub) / (DISTR.nua * (DISTR.nub + 2.));
  else
    DISTR.mode = 0.;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_F( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = _unur_lognormconstant_F(DISTR.params,DISTR.n_params);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_F( DISTR.domain[1],distr)
		 - _unur_cdf_F( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
double
_unur_lognormconstant_F(const double *params, int n_params ATTRIBUTE__UNUSED)
{ 
  return ((_unur_SF_ln_gamma(nua/2.) + _unur_SF_ln_gamma(nub/2.) - _unur_SF_ln_gamma((nua+nub)/2.))
	  - 0.5 * nua * log(nua/nub));
} 
int
_unur_set_params_F( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (nua <= 0. || nub <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.nua = nua;
  DISTR.nub = nub;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          
    DISTR.domain[1] = INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_F( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_F;
  distr->name = distr_name;
  DISTR.init = NULL;
  DISTR.pdf     = _unur_pdf_F;           
  DISTR.logpdf  = _unur_logpdf_F;        
  DISTR.dpdf    = _unur_dpdf_F;          
  DISTR.dlogpdf = _unur_dlogpdf_F;       
  DISTR.cdf     = _unur_cdf_F;           
#ifdef _unur_SF_invcdf_student
  DISTR.invcdf  = _unur_invcdf_F;        
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_F(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_lognormconstant_F(DISTR.params,DISTR.n_params);
  _unur_upd_mode_F( distr );
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_F;
  DISTR.upd_mode  = _unur_upd_mode_F;   
  DISTR.upd_area  = _unur_upd_area_F;   
  return distr;
} 
#undef nu1
#undef nu2
#undef DISTR
