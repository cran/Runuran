/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "vg";
#define lambda  params[0]    
#define alpha   params[1]    
#define beta    params[2]    
#define mu      params[3]    
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
#ifdef _unur_SF_bessel_k
static double _unur_pdf_vg( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_vg( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_center_vg( UNUR_DISTR *distr );
static double _unur_lognormconstant_vg( const double *params, int n_params );
static int _unur_set_params_vg( UNUR_DISTR *distr, const double *params, int n_params );
#ifdef _unur_SF_bessel_k
double
_unur_pdf_vg(double x, const UNUR_DISTR *distr)
{
  return exp(_unur_logpdf_vg(x,distr));
} 
double
_unur_logpdf_vg(double x, const UNUR_DISTR *distr)
{
  const double *params = DISTR.params;
  double nu = lambda - 0.5;   
  double res;                 
  double y, absy;             
  y = x - mu;
  absy = fabs(y);
  do {
    if (absy>0) {
      double besk;
      if (nu < 100) 
	besk = _unur_SF_ln_bessel_k(alpha*absy, nu);
      else
	besk = _unur_SF_bessel_k_nuasympt(alpha*absy, nu, TRUE, FALSE);
      if (_unur_isfinite(besk) && besk < MAXLOG - 20.0) {
	res = LOGNORMCONSTANT + besk + log(absy)*nu + beta*y;
	break;
      }
    }
    if (absy < 1.0) {
      res = LOGNORMCONSTANT + beta*y;
      res += -M_LN2 + _unur_SF_ln_gamma(nu) + nu*log(2./alpha);
      if (nu > 1.0) {
	double xi = 0.25*(alpha*absy)*(alpha*absy);
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
_unur_upd_center_vg( UNUR_DISTR *distr )
{
  const double *params = DISTR.params;
  double gam = sqrt(alpha*alpha-beta*beta);
  DISTR.center = mu + 2*beta*lambda / (gam*gam);
  if (!_unur_isfinite(DISTR.center))
    DISTR.center = mu;
  if (DISTR.center < DISTR.domain[0])
    DISTR.center = DISTR.domain[0];
  else if (DISTR.center > DISTR.domain[1])
    DISTR.center = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
double
_unur_lognormconstant_vg(const double *params, int n_params ATTRIBUTE__UNUSED)
{
  return (lambda*log(alpha*alpha - beta*beta) - 0.5*M_LNPI 
	  - (lambda-0.5)*log(2*alpha) - _unur_SF_ln_gamma(lambda));
} 
int
_unur_set_params_vg( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 4) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (alpha <= fabs(beta)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= |beta|");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.lambda = lambda;
  DISTR.alpha = alpha;
  DISTR.beta = beta;
  DISTR.mu = mu;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -UNUR_INFINITY;   
    DISTR.domain[1] = UNUR_INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_vg( const double *params, int n_params)
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_VG;
  distr->name = distr_name;
#ifdef _unur_SF_bessel_k
  DISTR.pdf     = _unur_pdf_vg;     
  DISTR.logpdf  = _unur_logpdf_vg;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_CENTER |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_vg(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_lognormconstant_vg(DISTR.params,DISTR.n_params);
  if (_unur_upd_center_vg(distr)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.area = 1;
  DISTR.set_params = _unur_set_params_vg;
  return distr;
} 
#undef mu
#undef alpha
#undef beta
#undef lambda
#undef DISTR
