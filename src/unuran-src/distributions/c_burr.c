/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "burr";
#define burr_type  params[0]
#define k          params[1]
#define c          params[2]
#define DISTR distr->data.cont
static double _unur_cdf_burr(double x, const UNUR_DISTR *distr);
static int _unur_set_params_burr( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_cdf_burr( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  switch ((int) (burr_type + 0.5)) {
  case  1: 
    if (x<=0.)
      return 0.;
    if (x>=1.)
      return 1.;
    return x;
  case  2: 
    return pow(exp(-x) + 1., -k);
  case  3: 
    if (x<=0.)
      return 0.;
    return pow( pow(x, -c) + 1., -k);
  case  4: 
    if (x<=0.)
      return 0.;
    if (x>=c)
      return 1.;
    return pow( pow( (c-x)/x, 1/c ) + 1., -k);
  case  5: 
    if (x<=-M_PI/2.)
      return 0.;
    if (x>= M_PI/2.)
      return 1.;
    return pow( c * exp(-tan(x)) + 1., -k );
  case  6: 
    return pow( c * exp(-k*sinh(x)) + 1., -k );
  case  7: 
    return pow( (1. + tanh(x))/2, k );
  case  8: 
    return pow( 2./M_PI * atan(exp(x)), k );
  case  9: 
    return (1. - 2. / (2. + c * (pow(1.+exp(x), k) - 1.)));
  case 10: 
    if (x<=0.) return 0.;
    return pow( 1. - exp(-x*x), k );
  case 11: 
    if (x<=0.)
      return 0.;
    if (x>=1.)
      return 1.;
    return pow( x - M_PI/2. * sin( 2. * M_PI * x), k );
  case 12: 
    if (x<=0.)
      return 0.;
    return (1. - pow( 1 + pow(x, c), -k ) );
  default:
    _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }
} 
int
_unur_set_params_burr( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  switch (distr->id) {
  case UNUR_DISTR_BURR_III:
  case UNUR_DISTR_BURR_IV:
  case UNUR_DISTR_BURR_V:
  case UNUR_DISTR_BURR_VI:
  case UNUR_DISTR_BURR_IX:
  case UNUR_DISTR_BURR_XII:
    if (n_params < 3) {
      _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few");
      free( distr ); return UNUR_ERR_DISTR_NPARAMS;
    }
  default: 
    if (n_params == 3) {
      _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
      n_params = 2;
    }
  }
  if (k <= 0. || (c <= 0. && n_params == 3) ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0 || c <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.burr_type = burr_type;
  switch (n_params) {
  case 3:
    DISTR.c = c;
  case 2:
    DISTR.k = k;
  default:
    break;
  }
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;  
    DISTR.domain[1] = INFINITY;   
    switch (distr->id) {
    case UNUR_DISTR_BURR_I:
      DISTR.domain[0] = 0.;       
      DISTR.domain[1] = 1.;       
      break;
    case UNUR_DISTR_BURR_III:
      DISTR.domain[0] = 0.;       
      break;
    case UNUR_DISTR_BURR_IV:
      DISTR.domain[0] = 0.;       
      DISTR.domain[1] = DISTR.c;  
      break;
    case UNUR_DISTR_BURR_V:
      DISTR.domain[0] = -M_PI/2.; 
      DISTR.domain[1] = M_PI/2.;  
      break;
    case UNUR_DISTR_BURR_X:
      DISTR.domain[0] = 0.;       
      break;
    case UNUR_DISTR_BURR_XI:
      DISTR.domain[0] = 0.;       
      DISTR.domain[1] = 1.;       
      break;
    case UNUR_DISTR_BURR_XII:
      DISTR.domain[0] = 0.;       
      break;
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_burr( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  switch ((int) (burr_type + 0.5)) {
  case  1:  distr->id = UNUR_DISTR_BURR_I;    break;
  case  2:  distr->id = UNUR_DISTR_BURR_II;   break;
  case  3:  distr->id = UNUR_DISTR_BURR_III;  break;
  case  4:  distr->id = UNUR_DISTR_BURR_IV;   break;
  case  5:  distr->id = UNUR_DISTR_BURR_V;    break;
  case  6:  distr->id = UNUR_DISTR_BURR_VI;   break;
  case  7:  distr->id = UNUR_DISTR_BURR_VII;  break;
  case  8:  distr->id = UNUR_DISTR_BURR_VIII; break;
  case  9:  distr->id = UNUR_DISTR_BURR_IX;   break;
  case 10:  distr->id = UNUR_DISTR_BURR_X;    break;
  case 11:  distr->id = UNUR_DISTR_BURR_XI;   break;
  case 12:  distr->id = UNUR_DISTR_BURR_XII;  break;
  default:
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"type < 1 || type > 12");
    free( distr ); return NULL;
  }
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_burr_init;
  DISTR.cdf  = _unur_cdf_burr;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN );
  if (_unur_set_params_burr(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.set_params = _unur_set_params_burr;
  return distr;
} 
#undef burr_type
#undef k
#undef c
#undef DISTR
