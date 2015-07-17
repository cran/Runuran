/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "lomax";
#define a params[0]
#define C params[1]
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_lomax( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_lomax( double x, const UNUR_DISTR *distr );
static double _unur_cdf_lomax( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_lomax( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_lomax( UNUR_DISTR *distr );
static int _unur_upd_area_lomax( UNUR_DISTR *distr );
static int _unur_set_params_lomax( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_lomax( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<0.) 
    return 0.;
  return ( pow(x+C,-(a+1.)) * NORMCONSTANT );
} 
double
_unur_dpdf_lomax( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<0.) ? 0. : -(a+1.) * pow(x+C,-(a+2.)) * NORMCONSTANT );
} 
double
_unur_cdf_lomax( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<0.) ? 0. : 1. - pow((C/(x+C)),a) );
} 
double
_unur_invcdf_lomax( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = C * ( pow(1-U, -1/a) - 1. );
  return X;
} 
int
_unur_upd_mode_lomax( UNUR_DISTR *distr )
{
  DISTR.mode = 0.;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_lomax( UNUR_DISTR *distr )
{
  NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_lomax( DISTR.domain[1],distr) 
		 - _unur_cdf_lomax( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_lomax( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (a <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params > 1 && C <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"C <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.a = a;
  DISTR.C = 1.; 
  if (n_params == 2)
    DISTR.C = C;
  DISTR.n_params = 2;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;               
    DISTR.domain[1] = UNUR_INFINITY;   
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_lomax( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_LOMAX;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_lomax;    
  DISTR.dpdf   = _unur_dpdf_lomax;   
  DISTR.cdf    = _unur_cdf_lomax;    
  DISTR.invcdf = _unur_invcdf_lomax; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_lomax(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);
  DISTR.mode = 0.;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_lomax;
  DISTR.upd_mode  = _unur_upd_mode_lomax; 
  DISTR.upd_area  = _unur_upd_area_lomax; 
  return distr;
} 
#undef a
#undef C
#undef DISTR
