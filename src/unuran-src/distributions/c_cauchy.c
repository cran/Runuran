/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "cauchy";
#define theta  params[0]
#define lambda params[1]
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_cdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_cauchy( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_cauchy( UNUR_DISTR *distr );
static int _unur_upd_area_cauchy( UNUR_DISTR *distr );
static int _unur_set_params_cauchy( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_cauchy(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / lambda; 
  return (1./((1+x*x)*NORMCONSTANT));
} 
double
_unur_logpdf_cauchy(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / lambda; 
  return (-log1p(x*x)-log(NORMCONSTANT)); 
} 
double
_unur_dpdf_cauchy(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / lambda; 
  return ( -2.*x/(lambda*(1.+x*x)*(1.+x*x)*NORMCONSTANT) );
} 
double
_unur_dlogpdf_cauchy(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - theta) / lambda; 
  return -2.*x/(lambda*(1.+x*x));
} 
double
_unur_cdf_cauchy(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  register double Fx;
  if (DISTR.n_params > 0)
    x = (x - theta) / lambda; 
  Fx = 0.5 + atan(x)/M_PI;
  if (Fx<0.)  Fx = 0.;
  if (Fx>1.)  Fx = 1.;
  return Fx;
} 
double
_unur_invcdf_cauchy(double u, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  double X;
  X = tan( M_PI * (u - 0.5) );
  return ((DISTR.n_params==0) ? X : theta + lambda * X );
} 
int
_unur_upd_mode_cauchy( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta; 
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_cauchy( UNUR_DISTR *distr )
{
  NORMCONSTANT = M_PI * DISTR.lambda;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_cauchy( DISTR.domain[1],distr) 
		 - _unur_cdf_cauchy( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_cauchy( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params == 2 && lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.theta  = 0.;
  DISTR.lambda = 1.;
  switch (n_params) {
  case 2:
    DISTR.lambda = lambda;
    /* FALLTHROUGH */
  case 1:
    DISTR.theta  = theta;
    n_params = 2;           
    /* FALLTHROUGH */
  default:
    break;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -UNUR_INFINITY;   
    DISTR.domain[1] = UNUR_INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_cauchy( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_CAUCHY;
  distr->name = distr_name;
  DISTR.pdf     = _unur_pdf_cauchy;     
  DISTR.logpdf  = _unur_logpdf_cauchy;  
  DISTR.dpdf    = _unur_dpdf_cauchy;    
  DISTR.dlogpdf = _unur_dlogpdf_cauchy; 
  DISTR.cdf     = _unur_cdf_cauchy;     
  DISTR.invcdf  = _unur_invcdf_cauchy;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_cauchy(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = M_PI * DISTR.lambda;
  DISTR.mode = DISTR.theta; 
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_cauchy;
  DISTR.upd_mode  = _unur_upd_mode_cauchy; 
  DISTR.upd_area  = _unur_upd_area_cauchy; 
  return distr;
} 
#undef theta 
#undef lambda
#undef DISTR
