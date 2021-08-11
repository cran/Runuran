/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "extremeII";
#define k      params[0]    
#define zeta   params[1]    
#define theta  params[2]    
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_extremeII( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_extremeII( double x, const UNUR_DISTR *distr );
static double _unur_cdf_extremeII( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_extremeII( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_extremeII( UNUR_DISTR *distr );
static int _unur_upd_area_extremeII( UNUR_DISTR *distr );
static int _unur_set_params_extremeII( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_extremeII( double x, const UNUR_DISTR *distr )
{ 
  register double xk;
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x - zeta) / theta;
  if (x<=0.)
    return 0.;
  xk = pow( x, -k - 1.);
  return ( exp( -xk * x - LOGNORMCONSTANT) * xk * k );
} 
double
_unur_dpdf_extremeII( double x, const UNUR_DISTR *distr )
{ 
  register double factor = 1.;
  register double xk;
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1) {
    factor = 1. / (theta * theta);
    x = (x - zeta) / theta;
  }
  if (x<=0.)
    return 0.;
  xk = pow(x, k);
  return (- factor * exp(-1./xk) * k * (xk + k*(xk - 1.)) / pow(x,2. + 2.*k)) ;
} 
double
_unur_cdf_extremeII( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 1)
    x = (x - zeta) / theta;
  if (x<=0.)
    return 0.;
  return ( exp( -pow( x, -k ) ) );
} 
double
_unur_invcdf_extremeII( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = exp( -log( -log(U) )/k );
  return ((DISTR.n_params==1) ? X : zeta + theta * X );
} 
int
_unur_upd_mode_extremeII( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.zeta + pow( DISTR.k/(DISTR.k+1.), 1/DISTR.k ) * DISTR.theta;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_extremeII( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = log(DISTR.theta);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_extremeII( DISTR.domain[1],distr) 
		 - _unur_cdf_extremeII( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_extremeII( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (k <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params > 2 && theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.k = k;
  DISTR.zeta  = 0.;
  DISTR.theta = 1.;
  switch (n_params) {
  case 3:
    DISTR.theta = theta;
    /* FALLTHROUGH */
  case 2:
    DISTR.zeta = zeta;
    n_params = 3;           
    /* FALLTHROUGH */
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
unur_distr_extremeII( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_EXTREME_II;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_extremeII;    
  DISTR.dpdf   = _unur_dpdf_extremeII;   
  DISTR.cdf    = _unur_cdf_extremeII;    
  DISTR.invcdf = _unur_invcdf_extremeII; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_extremeII(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = log(DISTR.theta);
  DISTR.mode = DISTR.zeta + pow( DISTR.k/(DISTR.k+1.), 1/DISTR.k ) * DISTR.theta;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_extremeII;
  DISTR.upd_mode  = _unur_upd_mode_extremeII; 
  DISTR.upd_area  = _unur_upd_area_extremeII; 
  return distr;
} 
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
