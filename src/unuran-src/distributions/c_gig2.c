/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
#ifdef USE_EXPERIMENTAL_CODE
#include <gsl/gsl_integration.h>
#endif
static const char distr_name[] = "gig2";
#define theta  params[0]    
#define psi    params[1]    
#define chi    params[2]    
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_gig2( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_gig2( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_gig2( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_gig2( double x, const UNUR_DISTR *distr );
#ifdef USE_EXPERIMENTAL_CODE
static double _unur_cdf_gig2( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_gig2( UNUR_DISTR *distr );
static double _unur_normconstant_gig2( const double *params, int n_params );
static int _unur_set_params_gig2( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return (NORMCONSTANT * exp( (theta-1.) * log(x) - 0.5 * (chi/x + psi*x) ));
} 
double
_unur_logpdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return -UNUR_INFINITY;
  return ( (theta-1.) * log(x) - 0.5 * (chi/x + psi*x) + log(NORMCONSTANT) );
} 
double
_unur_dpdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return ( NORMCONSTANT * 0.5 * exp( (theta-3.) * log(x) - (chi + psi*x*x)/(2*x) )
	   * (chi - x*(2 - 2*theta + psi*x)) );
} 
double
_unur_dlogpdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  return ( -0.5*(psi - chi/(x*x)) + (theta-1.)/x  + log(NORMCONSTANT) ) ;
} 
#ifdef USE_EXPERIMENTAL_CODE
#ifndef _unur_SF_bessel_k
#error run ./configure with flag --with-Rmath
#endif 
double
_unur_cdf_gig2(double x, const UNUR_DISTR *distr)
{
  double epsabs = 1.e-13;
  double epsrel = 1.e-13;
  double result, abserr;
  size_t limit = 1000;
  gsl_integration_workspace *work;
  if (x <= 0.)
    return 0.;
  work = gsl_integration_workspace_alloc (limit);
  gsl_function F;
  F.function = _unur_pdf_gig2;
  F.params = distr;
  gsl_integration_qag (&F, 0, x, epsabs, epsrel, limit,  GSL_INTEG_GAUSS61,
		       work, &result, &abserr);
  gsl_integration_workspace_free (work);
  return result;
} 
#endif
int
_unur_upd_mode_gig2( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  DISTR.mode =
    ((theta-1.)+sqrt((theta-1.)*(theta-1.) + psi*chi)) / psi;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
double
_unur_normconstant_gig2(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
#ifdef _unur_SF_bessel_k
  return ( pow(psi/chi, theta/2.) / (2. * _unur_SF_bessel_k(sqrt(psi*chi),theta)) );
#else
  return 1.;
#endif
} 
int
_unur_set_params_gig2( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); 
    return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (psi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"psi <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (chi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"chi <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.theta = theta;
  DISTR.psi = psi;
  DISTR.chi = chi;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;             
    DISTR.domain[1] = UNUR_INFINITY;  
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_gig2( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_GIG2;
  distr->name = distr_name;
  DISTR.pdf     = _unur_pdf_gig2;     
  DISTR.logpdf  = _unur_logpdf_gig2;  
  DISTR.dpdf    = _unur_dpdf_gig2;    
  DISTR.dlogpdf = _unur_dlogpdf_gig2; 
#ifdef USE_EXPERIMENTAL_CODE
  DISTR.cdf     = _unur_cdf_gig2;     
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   );
  if (_unur_set_params_gig2(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = _unur_normconstant_gig2(DISTR.params,DISTR.n_params);
  _unur_upd_mode_gig2(distr);
  DISTR.set_params = _unur_set_params_gig2;
  DISTR.upd_mode  = _unur_upd_mode_gig2; 
  return distr;
} 
#undef theta
#undef psi
#undef chi
#undef DISTR
