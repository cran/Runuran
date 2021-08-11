/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "uniform";
#define a  params[0]
#define b  params[1]
#define DISTR distr->data.cont
static double _unur_pdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_cdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_uniform( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_uniform( UNUR_DISTR *distr );
static int _unur_upd_area_uniform( UNUR_DISTR *distr );
static int _unur_set_params_uniform( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_uniform( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x < a || x > b)
    return 0.;
  return 1./(b-a);
} 
double
_unur_logpdf_uniform( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x < a || x > b)
    return -UNUR_INFINITY;
  return -log(b-a);
} 
double
_unur_dpdf_uniform( double x ATTRIBUTE__UNUSED, const UNUR_DISTR *distr ATTRIBUTE__UNUSED )
{ 
  return 0.;
} 
double
_unur_dlogpdf_uniform( double x ATTRIBUTE__UNUSED, const UNUR_DISTR *distr ATTRIBUTE__UNUSED )
{ 
  return 0.;
} 
double
_unur_cdf_uniform( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  x = (x-a) / (b-a);
  if (x<=0.)
    return 0.;
  if (x>=1.)
    return 1.;
  return x;
} 
double
_unur_invcdf_uniform( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = (DISTR.n_params==0) ? U : a + U * (b - a);
  return X;
} 
int
_unur_upd_mode_uniform( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_uniform( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_uniform( DISTR.domain[1],distr) 
		 - _unur_cdf_uniform( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_uniform( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params == 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params == 2 && (a >= b)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.a = 0.;
  DISTR.b = 1.;
  if (n_params == 2) {
    DISTR.a = a;
    DISTR.b = b;
  }
  DISTR.n_params = 2;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.a;      
    DISTR.domain[1] = DISTR.b;      
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_uniform( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_UNIFORM;
  distr->name = distr_name;
  DISTR.pdf     = _unur_pdf_uniform;     
  DISTR.logpdf  = _unur_logpdf_uniform;  
  DISTR.dpdf    = _unur_dpdf_uniform;    
  DISTR.dlogpdf = _unur_dlogpdf_uniform; 
  DISTR.cdf     = _unur_cdf_uniform;     
  DISTR.invcdf  = _unur_invcdf_uniform;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_uniform(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_uniform;
  DISTR.upd_mode  = _unur_upd_mode_uniform; 
  DISTR.upd_area  = _unur_upd_area_uniform; 
  return distr;
} 
#undef a
#undef b
#undef DISTR
