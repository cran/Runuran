/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "ghyp";
#define lambda  params[0]    
#define alpha   params[1]    
#define beta    params[2]    
#define delta   params[3]    
#define mu      params[4]    
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
#ifdef HAVE_BESSEL_K
static double _unur_pdf_ghyp( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_center_ghyp( UNUR_DISTR *distr );
static double _unur_normconstant_ghyp( const double *params, int n_params );
static int _unur_set_params_ghyp( UNUR_DISTR *distr, const double *params, int n_params );
#ifdef HAVE_BESSEL_K
double
_unur_pdf_ghyp(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  double tmp = delta*delta + (x-mu)*(x-mu);
  return ( NORMCONSTANT 
	   * pow( tmp, 0.5*lambda-0.25 ) 
	   * exp(beta*(x-mu))
	   * _unur_sf_bessel_k( alpha * sqrt(tmp), lambda-0.5 ) );
} 
#endif
int
_unur_upd_center_ghyp( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  DISTR.center = mu;
  if (DISTR.center < DISTR.domain[0]) 
    DISTR.center = DISTR.domain[0];
  else if (DISTR.center > DISTR.domain[1]) 
    DISTR.center = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
double
_unur_normconstant_ghyp(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
#ifdef HAVE_BESSEL_K
  double gamm = sqrt(alpha*alpha-beta*beta);
  return ( pow(gamm/delta, lambda ) 
	   / ( (M_SQRTPI*M_SQRT2) * pow(alpha, lambda-0.5)
	       * _unur_sf_bessel_k( delta*gamm, lambda ) ) );
#else
  return 1.;
#endif
} 
int
_unur_set_params_ghyp( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 5) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 5) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 5; }
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
  DISTR.lambda = lambda;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;   
    DISTR.domain[1] = INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_ghyp( const double *params, int n_params)
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_GHYP;
  distr->name = distr_name;
#ifdef HAVE_BESSEL_K
  DISTR.pdf     = _unur_pdf_ghyp;     
#endif
#ifdef HAVE_BESSEL_K
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_CENTER );
#else
  distr->set = ( UNUR_DISTR_SET_DOMAIN | UNUR_DISTR_SET_STDDOMAIN );
#endif
  if (_unur_set_params_ghyp(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = _unur_normconstant_ghyp(params,n_params);
  if (_unur_upd_center_ghyp(distr)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.set_params = _unur_set_params_ghyp;
  return distr;
} 
#undef mu
#undef alpha
#undef beta
#undef delta
#undef lambda
#undef DISTR
