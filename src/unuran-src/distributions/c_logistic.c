/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "logistic";
#define alpha  params[0]
#define beta   params[1]
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_logistic( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_logistic( double x, const UNUR_DISTR *distr );
static double _unur_cdf_logistic( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_logistic( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_logistic( UNUR_DISTR *distr );
static int _unur_upd_area_logistic( UNUR_DISTR *distr );
static int _unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_logistic( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double ex;
  if (DISTR.n_params > 0)
    x = (x - alpha) / beta;
  ex = exp( -fabs(x) );
  return (NORMCONSTANT * ex / ((1. + ex) * (1. + ex)));
} 
double
_unur_dpdf_logistic( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double factor = 1.;
  register double ex;
  if (DISTR.n_params > 0) {
    factor = 1. / beta;
    x = (x - alpha) / beta;
  }
  ex = exp(-fabs(x));
  if (x<0)
    factor = -factor;
  return (factor * NORMCONSTANT * ex * (ex - 1.) / ((1.+ex)*(1.+ex)*(1.+ex)));
} 
double
_unur_cdf_logistic( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (DISTR.n_params > 0)
    x = (x - alpha) / beta;
  return ( 1. / (1. + exp(-x)) );
} 
double
_unur_invcdf_logistic( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;
  X = -log(1./U - 1.);
  return ((DISTR.n_params==0) ? X : alpha + beta * X );
} 
int
_unur_upd_mode_logistic( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.alpha;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_logistic( UNUR_DISTR *distr )
{
  NORMCONSTANT = 1. / DISTR.beta;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_logistic( DISTR.domain[1],distr) 
		 - _unur_cdf_logistic( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params > 1 && beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.alpha = 0.;
  DISTR.beta  = 1.;
  switch (n_params) {
  case 2:
    DISTR.beta = beta;
    /* FALLTHROUGH */
  case 1:
    DISTR.alpha = alpha;
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
unur_distr_logistic( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_LOGISTIC;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_logistic;    
  DISTR.dpdf   = _unur_dpdf_logistic;   
  DISTR.cdf    = _unur_cdf_logistic;    
  DISTR.invcdf = _unur_invcdf_logistic; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_logistic(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.set_params = _unur_set_params_logistic;
  NORMCONSTANT = 1. / DISTR.beta;
  DISTR.mode = DISTR.alpha;
  DISTR.area = 1.;
  DISTR.upd_mode  = _unur_upd_mode_logistic; 
  DISTR.upd_area  = _unur_upd_area_logistic; 
  return distr;
} 
#undef alpha
#undef beta 
#undef DISTR
