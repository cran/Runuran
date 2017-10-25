/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/cvec.h>
#include <utils/matrix_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
#ifdef USE_DEPRECATED_CODE
#  include <distr/deprecated_distr.h>
#endif
static const char distr_name[] = "multinormal";
#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)
static double _unur_pdf_multinormal( const double *x, UNUR_DISTR *distr );
static double _unur_logpdf_multinormal( const double *x, UNUR_DISTR *distr );
static int _unur_dlogpdf_multinormal( double *result, const double *x, UNUR_DISTR *distr );
static double _unur_pdlogpdf_multinormal( const double *x, int coord, UNUR_DISTR *distr );
static int _unur_upd_mode_multinormal( UNUR_DISTR *distr );
static int _unur_upd_volume_multinormal( UNUR_DISTR *distr );
double
_unur_pdf_multinormal( const double *x, UNUR_DISTR *distr )
{ 
  return exp( _unur_logpdf_multinormal( x, distr ) );
} 
double
_unur_logpdf_multinormal( const double *x, UNUR_DISTR *distr )
{ 
#define idx(a,b) ((a)*dim+(b))
  int i,j, dim;
  double *mean;
  const double *covar_inv; 
  double xx; 
  double cx; 
  dim = distr->dim;
  if (DISTR.mean == NULL) {
    if (DISTR.covar != NULL) {
      _unur_warning(distr->name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }
    xx=0.;
    for (i=0; i<dim; i++) { xx += x[i]*x[i]; }
    return (-xx/2. + LOGNORMCONSTANT);  
  }
  mean = DISTR.mean;
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    return UNUR_INFINITY;
  xx=0.; 
  for (i=0; i<dim; i++) {
    cx=0.; 
    for (j=0; j<dim; j++) {
      cx += covar_inv[idx(i,j)] * (x[j]-mean[j]);
    }
    xx += (x[i]-mean[i])*cx;
  }
  return (-xx/2. + LOGNORMCONSTANT);
#undef idx
} 
int
_unur_dlogpdf_multinormal( double *result, const double *x, UNUR_DISTR *distr )
{
#define idx(a,b) ((a)*dim+(b))
  int i,j, dim;
  double *mean;
  const double *covar_inv;
  dim = distr->dim;
  mean = DISTR.mean;
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    return UNUR_FAILURE;
  for (i=0; i<dim; i++) {
    result[i] = 0.;
    for (j=0; j<dim; j++) 
      result[i] += -0.5 * (x[j]-mean[j]) * (covar_inv[idx(i,j)]+covar_inv[idx(j,i)]);
  }
  return UNUR_SUCCESS; 
#undef idx
} 
double
_unur_pdlogpdf_multinormal( const double *x, int coord, UNUR_DISTR *distr )
{
#define idx(a,b) ((a)*dim+(b))
  int j, dim;
  double *mean;
  const double *covar_inv;
  double result;
  dim = distr->dim;
  mean = DISTR.mean;
  if (coord < 0 || coord >= dim) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DOMAIN,"invalid coordinate");
    return UNUR_INFINITY;
  }
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    return UNUR_INFINITY;
  result = 0.;
  for (j=0; j<dim; j++)
    result += -0.5 * (x[j]-mean[j]) * (covar_inv[idx(coord,j)]+covar_inv[idx(j,coord)]);
  return result;
#undef idx
} 
int
_unur_upd_mode_multinormal( UNUR_DISTR *distr )
{
  if (DISTR.mode == NULL) 
    DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  memcpy( DISTR.mode, DISTR.mean, distr->dim * sizeof(double) );
  return UNUR_SUCCESS;
} 
int
_unur_upd_volume_multinormal( UNUR_DISTR *distr )
{
  double det_covar;
  det_covar = (DISTR.covar == NULL) 
    ? 1. : _unur_matrix_determinant(distr->dim, DISTR.covar);
  LOGNORMCONSTANT = - ( distr->dim * log(2 * M_PI) + log(det_covar) ) / 2.;
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_multinormal( int dim, const double *mean, const double *covar )
{
  struct unur_distr *distr;
  double det_covar; 
  distr = unur_distr_cvec_new(dim);
  if (distr == NULL) {
    return NULL;
  }
  distr->id = UNUR_DISTR_MNORMAL;
  distr->name = distr_name;
  DISTR.init = _unur_stdgen_multinormal_init;
  if ((unur_distr_cvec_set_mean(distr,mean)!=UNUR_SUCCESS) ||
      (unur_distr_cvec_set_covar(distr,covar)!=UNUR_SUCCESS) ) {
    unur_distr_free( distr );
    return NULL;
  }
  DISTR.pdf      = _unur_pdf_multinormal;       
  DISTR.logpdf   = _unur_logpdf_multinormal;    
  DISTR.dpdf     = _unur_distr_cvec_eval_dpdf_from_dlogpdf;  
  DISTR.dlogpdf  = _unur_dlogpdf_multinormal;    
  DISTR.pdpdf    = _unur_distr_cvec_eval_pdpdf_from_pdlogpdf;  
  DISTR.pdlogpdf = _unur_pdlogpdf_multinormal;  
#ifdef USE_DEPRECATED_CODE
  {
    struct unur_distr *stdmarginal = unur_distr_normal(NULL,0);
    unur_distr_cvec_set_stdmarginals(distr,stdmarginal);
    unur_distr_free(stdmarginal);
  }
#endif
  det_covar = (DISTR.covar == NULL) ? 1. : _unur_matrix_determinant(dim, DISTR.covar);
  LOGNORMCONSTANT = - ( distr->dim * log(2 * M_PI) + log(det_covar) ) / 2.;
  DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  memcpy( DISTR.mode, DISTR.mean, distr->dim * sizeof(double) );
  DISTR.volume = 1.; 
  DISTR.upd_mode   = _unur_upd_mode_multinormal;   
  DISTR.upd_volume = _unur_upd_volume_multinormal; 
  distr->set |= ( UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_PDFVOLUME |
		  UNUR_DISTR_SET_MODE );
  return distr;
} 
#undef DISTR
