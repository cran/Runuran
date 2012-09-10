/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
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
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
#ifdef _unur_SF_bessel_k
static double _unur_pdf_ghyp( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_ghyp( double x, const UNUR_DISTR *distr );
static double _unur_normconstant_ghyp( const double *params, int n_params );
#endif
static int _unur_upd_center_ghyp( UNUR_DISTR *distr );
static int _unur_set_params_ghyp( UNUR_DISTR *distr, const double *params, int n_params );
#ifdef _unur_SF_bessel_k
double
_unur_pdf_ghyp(double x, const UNUR_DISTR *distr)
{ 
  return exp(_unur_logpdf_ghyp(x,distr));
} 
double
_unur_logpdf_ghyp(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  double res = 0.;            
  double nu = lambda - 0.5;   
  double y;                   
  y = sqrt(delta*delta + (x-mu)*(x-mu));
  do {
    if (y>0.) {
      double besk;
      if (nu < 100)
        besk = _unur_SF_ln_bessel_k(alpha*y, nu);
      else
        besk = _unur_SF_bessel_k_nuasympt(alpha*y, nu, TRUE, FALSE);
      if (_unur_isfinite(besk) && besk < MAXLOG - 20.0) {
        res = LOGNORMCONSTANT + besk + nu*log(y) + beta*(x-mu);
        break;
      }
    }
    if (y < 1.0) {
      res = LOGNORMCONSTANT + beta*(x-mu);
      res += -M_LN2 + _unur_SF_ln_gamma(nu) + nu*log(2./alpha);
      if (nu > 1.0) {
        double xi = 0.25*(alpha*y)*(alpha*y);
        double sum = 1.0 - xi/(nu-1.0);
        if(nu > 2.0) sum += (xi/(nu-1.0)) * (xi/(nu-2.0));
        res += log(sum);
      }
    }
    else {
      res = -UNUR_INFINITY;
    }
  } while(0);
  return res;
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
#ifdef _unur_SF_bessel_k
double
_unur_normconstant_ghyp(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
  double gamm = sqrt(alpha*alpha-beta*beta);
  double logconst =  -0.5*(M_LNPI+M_LN2) + lambda * log(gamm/delta);
  logconst -= (lambda-0.5) * log(alpha);
  if (lambda < 50)
    logconst -= _unur_SF_ln_bessel_k( delta*gamm, lambda );
  else
    logconst -= _unur_SF_bessel_k_nuasympt( delta*gamm, lambda, TRUE, FALSE );
  return logconst;
} 
#endif
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
  DISTR.lambda = lambda;
  DISTR.alpha = alpha;
  DISTR.beta = beta;
  DISTR.delta = delta;
  DISTR.mu = mu;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -UNUR_INFINITY;   
    DISTR.domain[1] = UNUR_INFINITY;    
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
#ifdef _unur_SF_bessel_k
  DISTR.pdf     = _unur_pdf_ghyp;     
  DISTR.logpdf  = _unur_logpdf_ghyp;  
#endif
#ifdef _unur_SF_bessel_k
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_CENTER |
		 UNUR_DISTR_SET_PDFAREA );
#else
  distr->set = ( UNUR_DISTR_SET_DOMAIN | UNUR_DISTR_SET_STDDOMAIN );
#endif
  if (_unur_set_params_ghyp(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
#ifdef _unur_SF_bessel_k
  LOGNORMCONSTANT = _unur_normconstant_ghyp(DISTR.params,DISTR.n_params);
#else
  LOGNORMCONSTANT = 0.;
#endif
  if (_unur_upd_center_ghyp(distr)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
#ifdef _unur_SF_bessel_k
  DISTR.area = 1;
#endif
  DISTR.set_params = _unur_set_params_ghyp;
  return distr;
} 
#undef mu
#undef alpha
#undef beta
#undef delta
#undef lambda
#undef DISTR
