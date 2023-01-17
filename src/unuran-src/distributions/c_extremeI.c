/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "extremeI";
#define zeta   params[0]    
#define theta  params[1]    
#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_extremeI( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_extremeI( double x, const UNUR_DISTR *distr );
static double _unur_cdf_extremeI( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_extremeI( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_extremeI( UNUR_DISTR *distr );
static int _unur_upd_area_extremeI( UNUR_DISTR *distr );
static int _unur_set_params_extremeI( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_extremeI( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - zeta) / theta;
  return ( exp( -exp(-x) - x - LOGNORMCONSTANT) );
} 
double
_unur_dpdf_extremeI( double x, const UNUR_DISTR *distr )
{ 
  register double factor = 1.;
  register double expx;
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0) {
    factor = 1. / (theta * theta);
    x = (x - zeta) / theta;
  }
  expx = exp(-x);
  return ( exp( -expx - x ) * (expx - 1.) * factor );
} 
double
_unur_cdf_extremeI( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - zeta) / theta;
  return exp( -exp( -x) );
} 
double
_unur_invcdf_extremeI( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = -log( -log(U) );
  return ((DISTR.n_params==0) ? X : zeta + theta * X );
} 
int
_unur_upd_mode_extremeI( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.zeta;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_extremeI( UNUR_DISTR *distr )
{
  LOGNORMCONSTANT = log(DISTR.theta);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_extremeI( DISTR.domain[1],distr) 
		 - _unur_cdf_extremeI( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_extremeI( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params == 2 && theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.zeta  = 0.;
  DISTR.theta = 1.;
  switch (n_params) {
  case 2:
    DISTR.theta = theta;
    /* FALLTHROUGH */
  case 1:
    DISTR.zeta = zeta;
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
unur_distr_extremeI( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_EXTREME_I;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_extremeI;    
  DISTR.dpdf   = _unur_dpdf_extremeI;   
  DISTR.cdf    = _unur_cdf_extremeI;    
  DISTR.invcdf = _unur_invcdf_extremeI; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_extremeI(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  LOGNORMCONSTANT = log(DISTR.theta);
  DISTR.domain[0] = -UNUR_INFINITY;       
  DISTR.domain[1] = UNUR_INFINITY;        
  DISTR.mode = DISTR.zeta;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_extremeI;
  DISTR.upd_mode  = _unur_upd_mode_extremeI; 
  DISTR.upd_area  = _unur_upd_area_extremeI; 
  return distr;
} 
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
