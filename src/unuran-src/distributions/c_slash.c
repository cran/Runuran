/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "slash";
#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)
static double _unur_pdf_slash( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_slash( double x, const UNUR_DISTR *distr );
static int _unur_set_params_slash( UNUR_DISTR *distr, const double *params, int n_params );
double
_unur_pdf_slash(double x, const UNUR_DISTR *distr)
{
  if (_unur_iszero(x))
    return (0.5 * NORMCONSTANT);
  else
    return ((1. - exp(-x*x/2.)) / (x*x) * NORMCONSTANT);
} 
double
_unur_dpdf_slash(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
{ 
  register double xsq = x * x;
  if (_unur_iszero(x))
    return 0.;
  else
    return ((-2. + exp(-xsq/2.) * (2. + xsq)) / (xsq * x));
} 
int
_unur_set_params_slash( UNUR_DISTR *distr, const double *params ATTRIBUTE__UNUSED, int n_params )
{
  if (n_params > 0)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
  DISTR.n_params = 0;
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       
    DISTR.domain[1] = INFINITY;        
  }
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_slash( const double *params, int n_params )
{
  register struct unur_distr *distr;
  distr = unur_distr_cont_new();
  distr->id = UNUR_DISTR_SLASH;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_slash_init;
  DISTR.pdf  = _unur_pdf_slash;   
  DISTR.dpdf = _unur_dpdf_slash;  
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );
  if (_unur_set_params_slash(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  NORMCONSTANT = 1. / (M_SQRT2 * M_SQRTPI);
  DISTR.mode = 0.;
  DISTR.area = 1.;
  DISTR.set_params = _unur_set_params_slash;
  return distr;
} 
#undef nu
#undef DISTR
