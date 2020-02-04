/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "triangular";
#define H   params[0]    
#define DISTR distr->data.cont
static double _unur_pdf_triangular( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_triangular( double x, const UNUR_DISTR *distr );
static double _unur_cdf_triangular( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_triangular( double u, const UNUR_DISTR *distr );
static int _unur_upd_mode_triangular( UNUR_DISTR *distr );
static int _unur_upd_area_triangular( UNUR_DISTR *distr );
static int _unur_set_params_triangular( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_triangular( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x <= 0.)
    return 0.;
  if (x <= H)
    return (2.*x/H);
  if (x < 1.)
    return (2.*(1.-x)/(1.-H));
  return 0.;
} 
double
_unur_dpdf_triangular( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x < 0.)
    return 0.;
  if (x <= H && H > 0.)
    return (2./H);
  if (x <= 1.&& H < 1.)
    return (-2./(1.-H));
  return 0.;
} 
double
_unur_cdf_triangular( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  double Fx;
  if (x <= 0.)
    return 0.;
  if (x <= H)
    return (x*x/H);
  if (x < 1.) {
    if ((Fx = ((H + x * (x-2.))/(H-1.))) < 1.)
      return Fx;
  }
  return 1.;
} 
double
_unur_invcdf_triangular( double U, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  double tmp,X;
  if (U<=H) {
    X = sqrt(H*U);
  }
  else {
    tmp = (1.-H)*(1.-U);
    X = (tmp>0.) ? (1.-sqrt(tmp)) : 1.;
  }
  return X;
} 
int
_unur_upd_mode_triangular( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.H;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_triangular( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_triangular( DISTR.domain[1],distr) 
		 - _unur_cdf_triangular( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_triangular( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 0) n_params = 0;
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);
  if (n_params > 0 && (H < 0. || H > 1.)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"H < 0 || H > 1");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.H = 0.5;   
  if (n_params == 1)
    DISTR.H = H;
  DISTR.n_params = 1;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;        
    DISTR.domain[1] = 1.;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_triangular( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_TRIANGULAR;
  distr->name = distr_name;
  DISTR.pdf    = _unur_pdf_triangular;    
  DISTR.dpdf   = _unur_dpdf_triangular;   
  DISTR.cdf    = _unur_cdf_triangular;    
  DISTR.invcdf = _unur_invcdf_triangular; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_triangular(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = DISTR.H;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_triangular;
  DISTR.upd_mode  = _unur_upd_mode_triangular; 
  DISTR.upd_area  = _unur_upd_area_triangular; 
  return distr;
} 
#undef theta
#undef phi  
#undef DISTR
