/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "meixner";
#define alpha   params[0]    
#define beta    params[1]    
#define delta   params[2]    
#define mu      params[3]    
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_meixner( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_meixner( double x, const UNUR_DISTR *distr );
static int _unur_upd_center_meixner( UNUR_DISTR *distr );
static double _unur_lognormconstant_meixner( const double *params, int n_params );
static int _unur_set_params_meixner( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_meixner(double x, const UNUR_DISTR *distr)
{
  return exp(_unur_logpdf_meixner(x,distr));
} 
double
_unur_logpdf_meixner(double x, const UNUR_DISTR *distr)
{
  const double *params = DISTR.params;
  double res;           
  double y;             
  y = (x-mu) / alpha;
  res = LOGNORMCONSTANT + beta*y + 2*_unur_SF_Relcgamma(delta, y);
  return res;
} 
int
_unur_upd_center_meixner( UNUR_DISTR *distr )
{
  const double *params = DISTR.params;
  DISTR.center = mu;
  if (DISTR.center < DISTR.domain[0])
    DISTR.center = DISTR.domain[0];
  else if (DISTR.center > DISTR.domain[1])
    DISTR.center = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
double
_unur_lognormconstant_meixner(const double *params, int n_params ATTRIBUTE__UNUSED)
{
  return ( 2.*delta*log(2.*cos(beta/2.))
	   - (log(2.*alpha*M_PI) + _unur_SF_ln_gamma(2.*delta)));
} 
int
_unur_set_params_meixner( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 4) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (alpha <= 0. || delta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha or delta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (fabs(beta) >= M_PI) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta not in (-PI,PI)");
    return UNUR_ERR_DISTR_DOMAIN;
  }
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
unur_distr_meixner( const double *params, int n_params)
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_MEIXNER;
  distr->name = distr_name;
  DISTR.pdf     = _unur_pdf_meixner;     
  DISTR.logpdf  = _unur_logpdf_meixner;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_CENTER |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_meixner(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = _unur_lognormconstant_meixner(DISTR.params,DISTR.n_params);
  if (_unur_upd_center_meixner(distr)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.area = 1;
  DISTR.set_params = _unur_set_params_meixner;
  return distr;
} 
#undef alpha
#undef beta
#undef delta
#undef mu
#undef DISTR
