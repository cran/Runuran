/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "hyperbolic";
#define alpha  params[0]    
#define beta   params[1]    
#define delta  params[2]    
#define mu     params[3]    
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_hyperbolic( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_hyperbolic( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_hyperbolic( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_hyperbolic( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_hyperbolic( UNUR_DISTR *distr );
static double _unur_normconstant_hyperbolic( const double *params, int n_params );
static int _unur_set_params_hyperbolic( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  return (NORMCONSTANT * exp(-alpha * sqrt(delta*delta + (x-mu)*(x-mu)) + beta*(x-mu) ) );
} 
double
_unur_logpdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  return (-alpha * sqrt(delta*delta + (x-mu)*(x-mu)) + beta*(x-mu) + log(NORMCONSTANT) );
} 
double
_unur_dpdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  return (NORMCONSTANT * _unur_pdf_hyperbolic(x,distr) * _unur_dlogpdf_hyperbolic(x,distr));
} 
double
_unur_dlogpdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  return (beta - (alpha*(x-mu))/sqrt(delta*delta + (x-mu)*(x-mu)) + log(NORMCONSTANT));
} 
int
_unur_upd_mode_hyperbolic( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  DISTR.mode =
    mu + delta*beta / sqrt(alpha*alpha - beta*beta);
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
double
_unur_normconstant_hyperbolic(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
#ifdef HAVE_BESSEL_K
  double gamm = sqrt(alpha*alpha-beta*beta);
  return ( gamm / ( 2 * alpha * delta * _unur_sf_bessel_k(delta*gamm, 1) ) );
#else
  return 1.;
#endif
} 
int
_unur_set_params_hyperbolic( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 4) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (delta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"delta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (alpha <= fabs(beta)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= |beta|");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.mu = mu;
  DISTR.alpha = alpha;
  DISTR.beta = beta;
  DISTR.delta = delta;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;   
    DISTR.domain[1] = INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_hyperbolic( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_HYPERBOLIC;
  distr->name = distr_name;
  DISTR.pdf     = _unur_pdf_hyperbolic;     
  DISTR.logpdf  = _unur_logpdf_hyperbolic;  
  DISTR.dpdf    = _unur_dpdf_hyperbolic;    
  DISTR.dlogpdf = _unur_dlogpdf_hyperbolic; 
  DISTR.cdf  = NULL;                 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   );
  if (_unur_set_params_hyperbolic(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = _unur_normconstant_hyperbolic(params,n_params);
  _unur_upd_mode_hyperbolic(distr);
  DISTR.set_params = _unur_set_params_hyperbolic;
  DISTR.upd_mode  = _unur_upd_mode_hyperbolic; 
  return distr;
} 
#undef mu
#undef alpha
#undef beta
#undef delta
#undef DISTR
