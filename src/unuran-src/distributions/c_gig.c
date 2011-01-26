/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "gig";
#define theta  params[0]    
#define omega  params[1]    
#define eta    params[2]    
#define DISTR distr->data.cont
static double _unur_pdf_gig( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_gig( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_gig( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_gig( double x, const UNUR_DISTR *distr );
static int _unur_upd_mode_gig( UNUR_DISTR *distr );
static int _unur_set_params_gig( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return (exp( (theta-1.) * log(x) - 0.5 * omega * (x/eta + eta/x) ));
} 
double
_unur_logpdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return -INFINITY;
  return ( (theta-1.) * log(x) - 0.5 * omega * (x/eta + eta/x) );
} 
double
_unur_dpdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return ( exp( (theta-3.) * log(x) - 0.5 * omega * (x/eta + eta/x) )
	   * (eta*eta*omega + 2.*eta*(theta-1.)*x - omega*x*x) / (2*eta) );
} 
double
_unur_dlogpdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return ( -0.5*(omega*(1/eta - eta/(x*x))) + (theta-1.)/x );
} 
int
_unur_upd_mode_gig( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  DISTR.mode =
    (eta*(-1. + sqrt(omega*omega + (theta-1.)*(theta-1.)) + theta))/omega;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_set_params_gig( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (omega <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"omega <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (n_params == 3 && eta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"eta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.theta = theta;
  DISTR.omega = omega;
  DISTR.eta  = 1.;
  if (n_params > 2)
    DISTR.eta = eta;
  n_params = 3;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          
    DISTR.domain[1] = INFINITY;    
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_gig( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_GIG;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_gig_init;
  DISTR.pdf     = _unur_pdf_gig;     
  DISTR.logpdf  = _unur_logpdf_gig;  
  DISTR.dpdf    = _unur_dpdf_gig;    
  DISTR.dlogpdf = _unur_dlogpdf_gig; 
  DISTR.cdf  = NULL;                 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   );
  if (_unur_set_params_gig(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  _unur_upd_mode_gig(distr);
  DISTR.set_params = _unur_set_params_gig;
  DISTR.upd_mode  = _unur_upd_mode_gig; 
  return distr;
} 
#undef theta
#undef omega
#undef eta
#undef DISTR
