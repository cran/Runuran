/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "pareto";
#define k  params[0]
#define a  params[1]
#define DISTR distr->data.cont
static double _unur_pdf_pareto( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_pareto( double x, const UNUR_DISTR *distr );
static double _unur_cdf_pareto( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_pareto( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_pareto( UNUR_DISTR *distr );
static int _unur_upd_area_pareto( UNUR_DISTR *distr );
static int _unur_set_params_pareto( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_pareto( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<k)
    return 0.;
  return ((a/k) / pow(x/k, a + 1.) );
} 
double
_unur_dpdf_pareto( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<k) ? 0. : a * (-a-1.) / (k * k) * pow(x/k,-a-2.) );
} 
double
_unur_cdf_pareto( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<k) ? 0. : (1. - pow(k/x,a)) );
} 
double
_unur_invcdf_pareto( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = pow(1-U, -1/a);
  X *= k;
  return X;
} 
int
_unur_upd_mode_pareto( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.k;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_pareto( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_pareto( DISTR.domain[1],distr) 
		 - _unur_cdf_pareto( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_pareto( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (k <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (a <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.k = k;
  DISTR.a = a;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.k;         
    DISTR.domain[1] = UNUR_INFINITY;   
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_pareto( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_PARETO;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_pareto;    
  DISTR.dpdf   = _unur_dpdf_pareto;   
  DISTR.cdf    = _unur_cdf_pareto;    
  DISTR.invcdf = _unur_invcdf_pareto; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_pareto(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = DISTR.k;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_pareto;
  DISTR.upd_mode  = _unur_upd_mode_pareto; 
  DISTR.upd_area  = _unur_upd_area_pareto; 
  return distr;
} 
#undef k
#undef a
#undef DISTR
