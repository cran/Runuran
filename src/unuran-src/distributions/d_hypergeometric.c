/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "hypergeometric";
#define N  params[0]
#define M  params[1]
#define n  params[2]
#define DISTR distr->data.discr
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)
static double _unur_pmf_hypergeometric( int k, const UNUR_DISTR *distr );
#ifdef _unur_SF_cdf_hypergeometric
static double _unur_cdf_hypergeometric( int k, const UNUR_DISTR *distr ); 
#endif
#ifdef _unur_SF_invcdf_hypergeometric
static int    _unur_invcdf_hypergeometric( double u, const UNUR_DISTR *distr ); 
#endif
static int _unur_upd_mode_hypergeometric( UNUR_DISTR *distr );
static int _unur_upd_sum_hypergeometric( UNUR_DISTR *distr );
static int _unur_set_params_hypergeometric( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pmf_hypergeometric(int k, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  if ( k < _unur_max(0,(n-N+M-0.5)) || k > _unur_min(n,M)+0.5 ) 
    return 0.;
  else
    return exp( LOGNORMCONSTANT - _unur_SF_ln_factorial(k) - _unur_SF_ln_factorial(M-k) -
                _unur_SF_ln_factorial(n-k) - _unur_SF_ln_factorial(N-M-n+k) );
} 
#ifdef _unur_SF_cdf_hypergeometric
double
_unur_cdf_hypergeometric(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  return _unur_SF_cdf_hypergeometric(k,N,M,n);
} 
#endif
#ifdef _unur_SF_invcdf_hypergeometric
int
_unur_invcdf_hypergeometric(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;
  x = _unur_SF_invcdf_hypergeometric(u,N,M,n);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} 
#endif
int
_unur_upd_mode_hypergeometric( UNUR_DISTR *distr )
{
  DISTR.mode = (int) ( (DISTR.n + 1) * (DISTR.M + 1.) / (DISTR.N + 2.) );
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_sum_hypergeometric( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  LOGNORMCONSTANT = _unur_SF_ln_factorial(M) + _unur_SF_ln_factorial(N-M) + _unur_SF_ln_factorial(n) +
    _unur_SF_ln_factorial(N-n) - _unur_SF_ln_factorial(N);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
#ifdef _unur_SF_cdf_hypergeometric
  DISTR.sum = ( _unur_cdf_hypergeometric( DISTR.domain[1],distr) 
		 - _unur_cdf_hypergeometric( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif
} 
int
_unur_set_params_hypergeometric( UNUR_DISTR *distr, const double *params, int n_params )
{
  int nh;
  if (n_params < 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (M <= 0. || N <=0. || n <= 0. || n >= N || M >= N ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"M, N, n must be > 0 and n<N M<N");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  nh = (int)(N+0.5);
  if(fabs(nh-N)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.N = nh; 
  nh = (int)(M+0.5);
  if(fabs(nh-M)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.M = nh; 
  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.n = nh; 
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = (int) (_unur_max(0,(DISTR.n - DISTR.N + DISTR.M + 0.5)));  
    DISTR.domain[1] = (int) (_unur_min(DISTR.n, DISTR.M) + 0.5);                 
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_hypergeometric( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_discr_new();
  distr->id = UNUR_DISTR_HYPERGEOMETRIC;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_hypergeometric_init;
  DISTR.pmf  = _unur_pmf_hypergeometric;   
#ifdef _unur_SF_cdf_hypergeometric
  DISTR.cdf  = _unur_cdf_hypergeometric;   
#endif
#ifdef _unur_SF_invcdf_hypergeometric
  DISTR.invcdf = _unur_invcdf_hypergeometric;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_hypergeometric(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  _unur_upd_sum_hypergeometric( distr );
  _unur_upd_mode_hypergeometric(distr);
  DISTR.sum = 1.;
  DISTR.set_params = _unur_set_params_hypergeometric;
  DISTR.upd_mode = _unur_upd_mode_hypergeometric; 
  DISTR.upd_sum  = _unur_upd_sum_hypergeometric;  
  return distr;
} 
#undef N
#undef M
#undef n
#undef DISTR
