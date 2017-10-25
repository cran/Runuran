/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "geometric";
#define p  params[0]
#define DISTR distr->data.discr
static double _unur_pmf_geometric( int k, const UNUR_DISTR *distr );
static double _unur_cdf_geometric( int k, const UNUR_DISTR *distr ); 
static int    _unur_invcdf_geometric( double u, const UNUR_DISTR *distr ); 
static int _unur_upd_mode_geometric( UNUR_DISTR *distr );
static int _unur_upd_sum_geometric( UNUR_DISTR *distr );
static int _unur_set_params_geometric( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_geometric(int k, const UNUR_DISTR *distr)
{ 
  return ((k<0) ? 0. : DISTR.p * pow( 1. - DISTR.p, (double)k ));
} 
double
_unur_cdf_geometric(int k, const UNUR_DISTR *distr)
{ 
  return ((k<0) ? 0. : (1. - pow(1. - DISTR.p, k+1.)) );
} 
int
_unur_invcdf_geometric(double u, const UNUR_DISTR *distr)
{ 
  double x;
  if (_unur_isone(DISTR.p))
    return 0;
  x = ceil(log1p(-u) / log1p(-DISTR.p) - 1.);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} 
int
_unur_upd_mode_geometric( UNUR_DISTR *distr )
{
  DISTR.mode = 0;
  if (DISTR.mode < DISTR.domain[0] || DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = (DISTR.domain[0]<0) ? 0 : DISTR.domain[0];
  return UNUR_SUCCESS;
} 
int
_unur_upd_sum_geometric( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.sum = ( _unur_cdf_geometric( DISTR.domain[1],distr) 
		 - _unur_cdf_geometric( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
} 
int
_unur_set_params_geometric( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (p <= 0. || p >= 1.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.p = p;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;           
    DISTR.domain[1] = INT_MAX;     
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_geometric( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_GEOMETRIC;
  distr->name = distr_name;
  DISTR.pmf     = _unur_pmf_geometric;    
  DISTR.cdf     = _unur_cdf_geometric;    
  DISTR.invcdf  = _unur_invcdf_geometric; 
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE | 
		 UNUR_DISTR_SET_PMFSUM );
  if (_unur_set_params_geometric(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  DISTR.mode = 0;
  DISTR.sum = 1.;
  DISTR.set_params = _unur_set_params_geometric;
  DISTR.upd_mode = _unur_upd_mode_geometric; 
  DISTR.upd_sum  = _unur_upd_sum_geometric; 
  return distr;
} 
#undef p
#undef DISTR
