/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
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
static double _unur_cdf_lognormal( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_lognormal( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_lognormal( UNUR_DISTR *distr );
static int _unur_set_params_lognormal( UNUR_DISTR *distr, const double *params, int n_params );
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
double
_unur_cdf_lognormal( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  double z;
  if (x <= theta)
    return 0.;
  z = (log(x-theta)-zeta) / sigma;
  return _unur_SF_cdf_normal(z);
} 
double
_unur_invcdf_lognormal( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  return (theta + exp( _unur_SF_invcdf_normal(x) * sigma + zeta));
} 
int
_unur_upd_mode_lognormal( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  DISTR.mode = 
    exp(-sigma*sigma) * ( exp(zeta) + theta*exp(sigma*sigma) );
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_set_params_lognormal( UNUR_DISTR *distr, const double *params, int n_params )
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
  if (n_params == 3)
    DISTR.theta = theta;
  DISTR.n_params = 3;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.theta;     
    DISTR.domain[1] = UNUR_INFINITY;   
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
  DISTR.pdf    = _unur_pdf_lognormal;     
  DISTR.dpdf   = _unur_dpdf_lognormal;    
  DISTR.cdf    = _unur_cdf_lognormal;     
  DISTR.invcdf = _unur_invcdf_lognormal;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_lognormal(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = DISTR.sigma * sqrt(2.*M_PI);
  _unur_upd_mode_lognormal(distr);
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_lognormal;
  DISTR.upd_mode  = _unur_upd_mode_lognormal;   
  return distr;
} 
#undef zeta 
#undef sigma
#undef theta
#undef DISTR
