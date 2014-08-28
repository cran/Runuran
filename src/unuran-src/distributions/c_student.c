/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "student";
#define nu  params[0]
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_student( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_student( double x, const UNUR_DISTR *distr );
static double _unur_cdf_student( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_student
static double _unur_invcdf_student( double x, const UNUR_DISTR *distr );
#endif
static int _unur_upd_mode_student( UNUR_DISTR *distr );
static int _unur_upd_area_student( UNUR_DISTR *distr );
static double _unur_normconstant_student(const double *params, int n_params);
static int _unur_set_params_student( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_student( double x, const UNUR_DISTR *distr )
{
  const double *params = DISTR.params;
  return pow( (1. + x*x/nu), (-nu-1.)*0.5 ) / NORMCONSTANT;
}  
double
_unur_dpdf_student( double x, const UNUR_DISTR *distr )
{
  const double *params = DISTR.params;
  return ( (-nu-1.)*x/nu * pow( (1. + x*x/nu), (-nu-3.)*0.5 ) / NORMCONSTANT );
} 
double
_unur_cdf_student(double x, const UNUR_DISTR *distr)
{
#ifdef _unur_SF_cdf_student
  return _unur_SF_cdf_student(x,DISTR.nu);
#else
  const double *params = DISTR.params;
  double xx;
  if (_unur_iszero(nu))
    return 0.; 
  xx=1./(1.+x*x/nu);
  if (x>0)
    return 1-0.5*_unur_SF_incomplete_beta(xx,0.5*nu,0.5)/_unur_SF_incomplete_beta(1.,0.5*nu,0.5);
  else
    return   0.5*_unur_SF_incomplete_beta(xx,0.5*nu,0.5)/_unur_SF_incomplete_beta(1.,0.5*nu,0.5);
#endif
} 
#ifdef _unur_SF_invcdf_student
double
_unur_invcdf_student(double x, const UNUR_DISTR *distr)
{
  return _unur_SF_invcdf_student(x,DISTR.nu);
} 
#endif
int
_unur_upd_mode_student( UNUR_DISTR *distr )
{
  DISTR.mode = 0.;
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
_unur_upd_area_student( UNUR_DISTR *distr )
{
  NORMCONSTANT = _unur_normconstant_student(DISTR.params,DISTR.n_params);
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  DISTR.area = ( _unur_cdf_student( DISTR.domain[1],distr) 
		 - _unur_cdf_student( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} 
double
_unur_normconstant_student( const double *params, int n_params ATTRIBUTE__UNUSED )
{
  return( sqrt(M_PI * nu) * exp(_unur_SF_ln_gamma(0.5*nu) - _unur_SF_ln_gamma(0.5*(nu+1.))) );
} 
int
_unur_set_params_student( UNUR_DISTR *distr, const double *params, int n_params )
{
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  DISTR.nu = nu;
  DISTR.n_params = n_params;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -UNUR_INFINITY;       
    DISTR.domain[1] = UNUR_INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_student( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_STUDENT;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_student_init;
  DISTR.pdf  = _unur_pdf_student;  
  DISTR.dpdf = _unur_dpdf_student; 
  DISTR.cdf  = _unur_cdf_student;  
#ifdef _unur_SF_invcdf_student
  DISTR.invcdf  = _unur_invcdf_student;  
#endif
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
  if (_unur_set_params_student(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = _unur_normconstant_student(DISTR.params,DISTR.n_params);
  DISTR.mode = 0.;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_student;
  DISTR.upd_mode  = _unur_upd_mode_student; 
  DISTR.upd_area  = _unur_upd_area_student; 
  return distr;
} 
#undef nu
#undef DISTR
