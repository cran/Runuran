/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "weibull";
#define c      params[0]
#define alpha  params[1]
#define zeta   params[2]
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_weibull( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_weibull( double x, const UNUR_DISTR *distr );
static double _unur_cdf_weibull( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_weibull( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_weibull( UNUR_DISTR *distr );
static int _unur_upd_area_weibull( UNUR_DISTR *distr );
static int _unur_set_params_weibull( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_weibull( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x - zeta) / alpha;
  if (x < 0.)
    return 0.;
  if (_unur_iszero(x) && _unur_isone(c))
    return NORMCONSTANT;
  if (_unur_iszero(x) && !_unur_isone(c))
    return 0.;
  return (exp (-pow (x, c) + (c-1.) * log (x)) * NORMCONSTANT);
} 
double
_unur_dpdf_weibull( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double factor = 1.;
  register double xc;
  if (DISTR.n_params > 1) {
    factor = 1. / alpha;
    x = (x - zeta) / alpha;
  }
  if (x < 0.)
    return 0.;
  if (_unur_iszero(x) && _unur_isone(c))
    return 0.; 
  xc = -pow (x, c);
  return (exp (xc + (c-2.) * log (x)) * (-1. - c * (-xc-1.)) * NORMCONSTANT * factor);
} 
double
_unur_cdf_weibull( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x - zeta) / alpha;
  if (x <= 0.)
    return 0.;
  return (1. - exp(-pow (x, c)));
} 
double
_unur_invcdf_weibull( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = pow( -log(1.-U), 1./c );
  return ((DISTR.n_params==1) ? X : zeta + alpha * X );
} 
int
_unur_upd_mode_weibull( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.c<=1.) ? 0. : DISTR.alpha * pow((DISTR.c - 1.)/DISTR.c, 1./DISTR.c) + DISTR.zeta;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_weibull( UNUR_DISTR *distr )
{
  NORMCONSTANT = DISTR.c / DISTR.alpha;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_weibull( DISTR.domain[1],distr) 
		 - _unur_cdf_weibull( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_weibull( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (c <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"c <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params > 1 && alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.c = c;
  DISTR.alpha = 1.;
  DISTR.zeta  = 0.;
  switch (n_params) {
  case 3:
    DISTR.zeta = zeta;
  case 2:
    DISTR.alpha = alpha;
    n_params = 3;           
  default:
    break;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.zeta;      
    DISTR.domain[1] = UNUR_INFINITY;   
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_weibull( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_WEIBULL;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_weibull;    
  DISTR.dpdf   = _unur_dpdf_weibull;   
  DISTR.cdf    = _unur_cdf_weibull;    
  DISTR.invcdf = _unur_invcdf_weibull; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_weibull(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = DISTR.c / DISTR.alpha;
  DISTR.mode = (DISTR.c<=1.) ? 0. : DISTR.alpha * pow((DISTR.c - 1.)/DISTR.c, 1./DISTR.c) + DISTR.zeta;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_weibull;
  DISTR.upd_mode  = _unur_upd_mode_weibull; 
  DISTR.upd_area  = _unur_upd_area_weibull; 
  return distr;
} 
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
