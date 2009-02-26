/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "lognormal";
#define zeta   params[0]
#define sigma  params[1]
#define theta  params[2]
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_lognormal( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_lognormal( double x, const UNUR_DISTR *distr );
static int _unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_lognormal( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double z;
  if (x <= theta)
    return 0.;
  z = log(x-theta)-zeta;
  return ( 1./(x-theta) * exp( -z*z/(2.*sigma*sigma) ) / NORMCONSTANT );
} 
double
_unur_dpdf_lognormal( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double z, sigmasqu;
  if (x <= theta)
    return 0.;
  z = log(x-theta)-zeta;
  sigmasqu = sigma * sigma;
  return ( 1/((x-theta)*(x-theta)) * exp( -z*z/(2*sigmasqu) ) * (1.+z/sigmasqu) / NORMCONSTANT );
} 
int
_unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.zeta = zeta;
  DISTR.sigma = sigma;
  DISTR.theta = 0.;        
  switch (n_params) {
  case 3:
    DISTR.theta = theta;
  default:
    n_params = 3;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.theta;     
    DISTR.domain[1] = INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_lognormal( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_LOGNORMAL;
  distr->name = distr_name;
  DISTR.init = NULL;         
  DISTR.pdf  = _unur_pdf_lognormal;  
  DISTR.dpdf = _unur_dpdf_lognormal; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_logistic(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = DISTR.sigma * sqrt(2.*M_PI);
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_logistic;
  return distr;
} 
#undef zeta 
#undef sigma
#undef theta
#undef DISTR
