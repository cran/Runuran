/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "beta";
#define p  params[0]
#define q  params[1]
#define a  params[2]
#define b  params[3]
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_cdf_beta( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_beta
static double _unur_invcdf_beta( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_beta( UNUR_DISTR *distr );
static int _unur_upd_area_beta( UNUR_DISTR *distr );
inline static double _unur_lognormconstant_beta( const double *params, int n_params );
static int _unur_set_params_beta( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 2)
    x = (x-a) / (b-a);
  if (x > 0. && x < 1.)
    return exp((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);
  if ((_unur_iszero(x) && _unur_isone(p)) || (_unur_isone(x) && _unur_isone(q)))
    return exp(-LOGNORMCONSTANT);
  if ((_unur_iszero(x) && p<1.) || (_unur_isone(x) && q<1.))
    return UNUR_INFINITY;
  return 0.;
} 
double
_unur_logpdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 2)
    x = (x-a) / (b-a);
  if (x > 0. && x < 1.)
    return ((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);
  if ((_unur_iszero(x) && _unur_isone(p)) 
      || (_unur_isone(x) && _unur_isone(q)))
    return (-LOGNORMCONSTANT);
  if ((_unur_iszero(x) && p<1.) || (_unur_isone(x) && q<1.))
    return UNUR_INFINITY;
  return -UNUR_INFINITY;
} 
double
_unur_dpdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 2) {
    x = (x-a) / (b-a);
  }
  if (x > 0. && x < 1.)
    return (exp((p-2.)*log(x) + (q-2.)*log(1.-x) - LOGNORMCONSTANT) * ( (p-1.)*(1.-x) - (q-1.)*x ) / (b-a) );
  if (_unur_iszero(x) && _unur_isone(p))
    return (1.-q)*exp(-LOGNORMCONSTANT)/(b-a);
  if (_unur_iszero(x) && _unur_isfsame(p,2.))
    return exp(-LOGNORMCONSTANT)/(b-a);
  if (_unur_iszero(x) && p<2.)
    return (p>1. ? UNUR_INFINITY : -UNUR_INFINITY);
  if (_unur_isone(x) && _unur_isone(q))
    return (p-1.)*exp(-LOGNORMCONSTANT)/(b-a);
  if (_unur_isone(x) && _unur_isfsame(q,2.))
    return -exp(-LOGNORMCONSTANT)/(b-a);
  if (_unur_isone(x) && q<2.)
    return (q>1. ? -UNUR_INFINITY : UNUR_INFINITY);
  return 0.;
} 
double
_unur_dlogpdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 2) {
    x = (x-a) / (b-a);
  }
  if (x > 0. && x < 1.)
    return (((p-1.)/x - (q-1.)/(1.-x)) / (b-a));
  if (_unur_iszero(x) && p<1.)
    return -UNUR_INFINITY;
  if (_unur_iszero(x) && _unur_isone(p))
    return (-(q-1.)/((1.-x)*(b-a)));
  if (_unur_iszero(x) && p>1.)
    return UNUR_INFINITY;
  if (_unur_isone(x) && q<1.)
    return UNUR_INFINITY;
  if (_unur_isone(x) && _unur_isone(q))
    return ((p-1.)/(b-a));
  if (_unur_isone(x) && q>1.)
    return -UNUR_INFINITY;
  return 0.;
} 
double
_unur_cdf_beta(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 2)
    x = (x-a) / (b-a);
  if (x <= 0.) return 0.;
  if (x >= 1.) return 1.;
  return _unur_SF_incomplete_beta(x,p,q);
} 
#ifdef _unur_SF_invcdf_beta
double
_unur_invcdf_beta(double x, const UNUR_DISTR *distr)
{
  const double *params = DISTR.params;
  if (DISTR.n_params == 2)
    return _unur_SF_invcdf_beta(x,p,q);
  else
    return (a + _unur_SF_invcdf_beta(x,p,q))*(b-a);
} 
#endif
int
_unur_upd_mode_beta( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  if (p <= 1. && q > 1.)
    DISTR.mode = 0.;              
  else if (p > 1. && q <= 1.)
    DISTR.mode = 1.;              
  else if (p > 1. && q > 1.)
    DISTR.mode = (p - 1.) / (p + q - 2.);
  else {
    DISTR.mode = UNUR_INFINITY;
    return UNUR_ERR_DISTR_PROP;
  }
  if (DISTR.n_params > 2)
    DISTR.mode = DISTR.mode * (b - a) + a;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_beta( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_beta( DISTR.domain[1],distr) 
		 - _unur_cdf_beta( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
double
_unur_lognormconstant_beta(const double *params, int n_params)
{ 
  if (n_params > 2)
    return (_unur_SF_ln_gamma(p) + _unur_SF_ln_gamma(q) - _unur_SF_ln_gamma(p+q) + log(b-a) );
  else
    return (_unur_SF_ln_gamma(p) + _unur_SF_ln_gamma(q) - _unur_SF_ln_gamma(p+q));
} 
int
_unur_set_params_beta( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params == 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"");
    n_params = 2; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (p <= 0. || q <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 or q <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params > 2 && a >= b) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.p = p;
  DISTR.q = q;
  if (n_params > 2) {
    DISTR.a = a;
    DISTR.b = b;
  }
  else { 
    DISTR.a = 0.;      
    DISTR.b = 1.;      
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.a; 
    DISTR.domain[1] = DISTR.b; 
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_beta( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_BETA;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_beta_init;
  DISTR.pdf     = _unur_pdf_beta;     
  DISTR.logpdf  = _unur_logpdf_beta;  
  DISTR.dpdf    = _unur_dpdf_beta;    
  DISTR.dlogpdf = _unur_dlogpdf_beta; 
  DISTR.cdf     = _unur_cdf_beta;     
#ifdef _unur_SF_invcdf_beta
  DISTR.invcdf  = _unur_invcdf_beta;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_beta(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);
  _unur_upd_mode_beta( distr );
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_beta;
  DISTR.upd_mode  = _unur_upd_mode_beta; 
  DISTR.upd_area  = _unur_upd_area_beta; 
  return distr;
} 
#undef p
#undef q
#undef a
#undef b
#undef DISTR
