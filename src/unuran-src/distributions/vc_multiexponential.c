/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/cvec.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "multiexponential";
#define INDEX_SIGMA 0 
#define INDEX_THETA 1
#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)
static double _unur_pdf_multiexponential( const double *x, UNUR_DISTR *distr );
static double _unur_logpdf_multiexponential( const double *x, UNUR_DISTR *distr );
static int _unur_dlogpdf_multiexponential( double *result, const double *x, UNUR_DISTR *distr );
static int _unur_set_params_multiexponential( UNUR_DISTR *distr, const double *sigma, const double *theta );
static int _unur_upd_mode_multiexponential( UNUR_DISTR *distr );
static int _unur_upd_volume_multiexponential( UNUR_DISTR *distr );
double
_unur_pdf_multiexponential( const double *x, UNUR_DISTR *distr )
{ 
  double flog;
  flog=_unur_logpdf_multiexponential( x, distr );
  if (_unur_isfinite(flog))   
    return exp( flog );
  return 0.;
} 
double
_unur_logpdf_multiexponential( const double *x, UNUR_DISTR *distr )
{ 
  int i, dim;
  double dx, sum; 
  double *sigma, *theta;
  dim = distr->dim ;
  dx=0.;
  sum=0.;
  sigma = DISTR.param_vecs[INDEX_SIGMA];
  theta = DISTR.param_vecs[INDEX_THETA];
  if ( sigma==NULL || theta==NULL ) {
    for (i=0; i<dim; i++) { 
      dx = ((i==0)
	    ? ((x[i]<0)? UNUR_INFINITY: x[i]) 
	    : ((x[i]<x[i-1])? UNUR_INFINITY: x[i]-x[i-1]));
      sum -= (dim-i) * dx;  
    }
  }
  else {
    for (i=0; i<dim; i++) {
      dx = ((i==0) 
	    ? ( (x[i]-theta[i]<0) ? UNUR_INFINITY : x[i]-theta[i])
	    : (((x[i]-theta[i]) < (x[i-1]-theta[i-1])) ? UNUR_INFINITY : x[i]-x[i-1]-theta[i]+theta[i-1])); 
      dx /= sigma[i];
      sum -= (dim-i) * dx;  
    }
  }
  return ( sum + LOGNORMCONSTANT);
} 
int
_unur_dlogpdf_multiexponential( double *result, const double *x, UNUR_DISTR *distr )
{
  int i, dim;
  double dx, fx1, fx2;
  double *xx;
  dim = distr->dim;
  xx=malloc(dim*sizeof(double)); 
  dx=1e7*UNUR_EPSILON;   
  for (i=0; i<dim; i++) {
    memcpy(xx,x,dim*sizeof(double));
    xx[i]=x[i]+dx;
    fx1=_unur_logpdf_multiexponential( x, distr );
    fx2=_unur_logpdf_multiexponential(xx, distr );
    result[i] = (fx2-fx1)/dx;
  }
  if (xx) free(xx);
  return UNUR_SUCCESS; 
} 
int
_unur_set_params_multiexponential( UNUR_DISTR *distr, const double *sigma, const double *theta )
{
  int i;
  double *default_sigma=NULL;
  double *default_theta=NULL;
  if(sigma==NULL) {
    default_sigma = _unur_xmalloc( distr->dim * sizeof(double));
    for (i=0; i<distr->dim; i++) default_sigma[i]=1.;
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_SIGMA, default_sigma, distr->dim );
    if (default_sigma) free(default_sigma);
  }
  else {
    for (i=0; i<distr->dim; i++) {
      if ( sigma[i] <= UNUR_EPSILON ) {
        _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma is too low");
        return UNUR_ERR_DISTR_DOMAIN;
      }
    }
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_SIGMA, sigma, distr->dim );
  }
  if(theta==NULL) {
    default_theta = _unur_xmalloc(distr->dim * sizeof(double) );
    for (i=0; i<distr->dim; i++) default_theta[i]=0.;
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_THETA, default_theta, distr->dim );
    if (default_theta) free(default_theta);  
  }
  else {
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_THETA, theta, distr->dim ); 
  }
  DISTR.n_params = 0; 
  return UNUR_SUCCESS;
} 
int
_unur_upd_mode_multiexponential( UNUR_DISTR *distr )
{
  int i;
  if (DISTR.mode == NULL) 
    DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  for (i=0; i<distr->dim; i++)  DISTR.mode[i]=0.;
  return UNUR_SUCCESS;
} 
int
_unur_upd_volume_multiexponential( UNUR_DISTR *distr )
{
  int i;
  double sumsigma = 0.;
  for (i=0; i<distr->dim; i++) {
    sumsigma += DISTR.param_vecs[INDEX_SIGMA][i];
  }
  LOGNORMCONSTANT = - 1. / sumsigma;   
  return UNUR_SUCCESS;
} 
struct unur_distr *
unur_distr_multiexponential( int dim, const double *sigma, const double *theta )
{
  struct unur_distr *distr;
  struct unur_distr **marginal;
  int i;
  double sumsigma; 
  double alpha;  
  distr = unur_distr_cvec_new(dim);
  if (distr == NULL) {
    return NULL;
  }
  distr->id = UNUR_DISTR_MEXPONENTIAL;
  distr->name = distr_name;
  DISTR.init = NULL;
  DISTR.pdf     = _unur_pdf_multiexponential;       
  DISTR.logpdf  = _unur_logpdf_multiexponential;    
  DISTR.dpdf    = _unur_distr_cvec_eval_dpdf_from_dlogpdf;  
  DISTR.dlogpdf = _unur_dlogpdf_multiexponential;    
  DISTR.pdpdf   = _unur_distr_cvec_eval_pdpdf_from_pdlogpdf;  
  marginal = malloc(distr->dim * sizeof(struct unur_distr*));
  for (i=0; i<distr->dim; i++) {  
    alpha = i+1.; 
    marginal[i] = unur_distr_gamma(&alpha, 1);
  }
  unur_distr_cvec_set_marginal_array(distr, marginal);
  for (i=0; i<distr->dim; i++)
    if (marginal[i]) _unur_distr_free(marginal[i]);
  if (marginal) free(marginal);
  if (_unur_set_params_multiexponential(distr, sigma, theta)!=UNUR_SUCCESS) {
    _unur_distr_free(distr); return NULL;
  }
  sumsigma = 0.;
  for (i=0; i<distr->dim; i++) {
    sumsigma += DISTR.param_vecs[INDEX_SIGMA][i];
  }
  LOGNORMCONSTANT = - 1. / sumsigma;   
  DISTR.mode = _unur_xmalloc(distr->dim * sizeof(double) );
  for (i=0; i<distr->dim; i++)  DISTR.mode[i]=0.;
  DISTR.volume = 1.; 
  DISTR.upd_mode   = _unur_upd_mode_multiexponential;   
  DISTR.upd_volume = _unur_upd_volume_multiexponential; 
  distr->set |= ( UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_PDFVOLUME |
		  UNUR_DISTR_SET_MODE );
  return distr;
} 
#undef DISTR
