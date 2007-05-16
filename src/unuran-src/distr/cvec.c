/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "cvec.h"
#include <utils/matrix_source.h>
#include <stdarg.h>
static struct unur_distr **_unur_distr_cvec_marginals_clone ( struct unur_distr **marginals, int dim );
static void _unur_distr_cvec_marginals_free ( struct unur_distr **marginals, int dim );
static void _unur_distr_cvec_free( struct unur_distr *distr );
#define DISTR distr->data.cvec
struct unur_distr *
unur_distr_cvec_new( int dim )
{
  register struct unur_distr *distr;
  int i;
  if (dim < 1) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"dimension < 1");
    return NULL;
  }
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;
  COOKIE_SET(distr,CK_DISTR_CVEC);
  distr->type = UNUR_DISTR_CVEC;
  distr->id = UNUR_DISTR_GENERIC;
  distr->dim = dim;   
  distr->base = NULL;
  distr->destroy = _unur_distr_cvec_free;
  distr->clone = _unur_distr_cvec_clone;
  DISTR.pdf       = NULL;   
  DISTR.dpdf      = NULL;   
  DISTR.pdpdf     = NULL;   
  DISTR.logpdf    = NULL;   
  DISTR.dlogpdf   = NULL;   
  DISTR.pdlogpdf  = NULL;   
  DISTR.domainrect = NULL;  
  DISTR.init      = NULL;   
  DISTR.mean      = NULL;   
  DISTR.covar     = NULL;   
  DISTR.cholesky  = NULL;   
  DISTR.covar_inv = NULL;   
  DISTR.rankcorr  = NULL;   
  DISTR.rk_cholesky = NULL; 
  DISTR.marginals = NULL;   
  DISTR.upd_mode  = NULL;   
  DISTR.upd_volume = NULL;  
#ifdef USE_DEPRECATED_CODE
  DISTR.stdmarginals = NULL;  
#endif
  DISTR.n_params  = 0;             
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    DISTR.n_param_vec[i] = 0;
    DISTR.param_vecs[i] = NULL;
  }
  DISTR.norm_constant = 1.;        
  DISTR.mode       = NULL;         
  DISTR.center     = NULL;         
  DISTR.volume     = INFINITY;     
  return distr;
} 
struct unur_distr *
_unur_distr_cvec_clone( const struct unur_distr *distr )
{
#define CLONE clone->data.cvec
  struct unur_distr *clone;
  int i;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  memcpy( clone, distr, sizeof( struct unur_distr ) );
  if (DISTR.domainrect) {
    CLONE.domainrect = _unur_xmalloc( 2 * distr->dim * sizeof(double) );
    memcpy( CLONE.domainrect, DISTR.domainrect, 2 * distr->dim * sizeof(double) );
  }
  if (DISTR.mean) {
    CLONE.mean = _unur_xmalloc( distr->dim * sizeof(double) );
    memcpy( CLONE.mean, DISTR.mean, distr->dim * sizeof(double) );
  }
  if (DISTR.covar) {
    CLONE.covar = _unur_xmalloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.covar, DISTR.covar, distr->dim * distr->dim * sizeof(double) );
  }
  if (DISTR.cholesky) {
    CLONE.cholesky = _unur_xmalloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.cholesky, DISTR.cholesky, distr->dim * distr->dim * sizeof(double) );
  }
  if (DISTR.covar_inv) {
    CLONE.covar_inv = _unur_xmalloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.covar_inv, DISTR.covar_inv, distr->dim * distr->dim * sizeof(double) );
  }
  if (DISTR.rankcorr) {
    CLONE.rankcorr = _unur_xmalloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.rankcorr, DISTR.rankcorr, distr->dim * distr->dim * sizeof(double) );
  }
  if (DISTR.rk_cholesky) {
    CLONE.rk_cholesky = _unur_xmalloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.rk_cholesky, DISTR.rk_cholesky, distr->dim * distr->dim * sizeof(double) );
  }
  if (DISTR.mode) {
    CLONE.mode = _unur_xmalloc( distr->dim * sizeof(double) );
    memcpy( CLONE.mode, DISTR.mode, distr->dim * sizeof(double) );
  }
  if (DISTR.center) {
    CLONE.center = _unur_xmalloc( distr->dim * sizeof(double) );
    memcpy( CLONE.center, DISTR.center, distr->dim * sizeof(double) );
  }
  if (DISTR.marginals)
    CLONE.marginals = _unur_distr_cvec_marginals_clone( DISTR.marginals, distr->dim );
#ifdef USE_DEPRECATED_CODE
  if (DISTR.stdmarginals)
    CLONE.stdmarginals = _unur_distr_cvec_marginals_clone( DISTR.stdmarginals, distr->dim );
#endif
  CLONE.n_params = DISTR.n_params;  
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    CLONE.params[i] = DISTR.params[i];
  }
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    CLONE.n_param_vec[i] = DISTR.n_param_vec[i];
    if (DISTR.param_vecs[i]) {
      CLONE.param_vecs[i] = _unur_xmalloc( DISTR.n_param_vec[i] * sizeof(double) );
      memcpy( CLONE.param_vecs[i], DISTR.param_vecs[i], DISTR.n_param_vec[i] * sizeof(double) );
    }
  }
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }
  return clone;
#undef CLONE
} 
void
_unur_distr_cvec_free( struct unur_distr *distr )
{
  int i;
  if( distr == NULL ) 
    return;
  COOKIE_CHECK(distr,CK_DISTR_CVEC,RETURN_VOID);
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    if (DISTR.param_vecs[i]) free( DISTR.param_vecs[i] );
  if (DISTR.domainrect)  free(DISTR.domainrect); 
  if (DISTR.mean)        free(DISTR.mean); 
  if (DISTR.covar)       free(DISTR.covar);
  if (DISTR.covar_inv)   free(DISTR.covar_inv);
  if (DISTR.cholesky)    free(DISTR.cholesky);
  if (DISTR.rankcorr)    free(DISTR.rankcorr);
  if (DISTR.rk_cholesky) free(DISTR.rk_cholesky);
  if (DISTR.mode)        free(DISTR.mode);
  if (DISTR.center)      free(DISTR.center);
  if (DISTR.marginals)
    _unur_distr_cvec_marginals_free(DISTR.marginals, distr->dim);
#ifdef USE_DEPRECATED_CODE
  if (DISTR.stdmarginals)
    _unur_distr_cvec_marginals_free(DISTR.stdmarginals, distr->dim);
#endif
  if (distr->name_str) free(distr->name_str);
  COOKIE_CLEAR(distr);
  free( distr );
} 
int
unur_distr_cvec_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CVEC *pdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdf, UNUR_ERR_NULL);
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.pdf = pdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_set_dpdf( struct unur_distr *distr, UNUR_VFUNCT_CVEC *dpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.dpdf != NULL || DISTR.dlogpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.dpdf = dpdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_set_pdpdf( struct unur_distr *distr, UNUR_FUNCTD_CVEC *pdpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdpdf, UNUR_ERR_NULL);
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdpdf != NULL || DISTR.pdlogpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of pdPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.pdpdf = pdpdf;
  return UNUR_SUCCESS;
} 
UNUR_FUNCT_CVEC *
unur_distr_cvec_get_pdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  return DISTR.pdf;
} 
UNUR_VFUNCT_CVEC *
unur_distr_cvec_get_dpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  return DISTR.dpdf;
} 
UNUR_FUNCTD_CVEC *
unur_distr_cvec_get_pdpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  return DISTR.pdpdf;
} 
double
unur_distr_cvec_eval_pdf( const double *x, struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );
  if (DISTR.pdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cvec_PDF(x,distr);
} 
int
unur_distr_cvec_eval_dpdf( double *result, const double *x, struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.dpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  return _unur_cvec_dPDF(result,x,distr);
} 
double
unur_distr_cvec_eval_pdpdf( const double *x, int coord, struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );
  if (DISTR.pdpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  if (coord < 0 || coord >= distr->dim) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DOMAIN,"invalid coordinate");
    return INFINITY;
  }
  return _unur_cvec_pdPDF(x,coord,distr);
} 
int
unur_distr_cvec_set_logpdf( struct unur_distr *distr, UNUR_FUNCT_CVEC *logpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, logpdf, UNUR_ERR_NULL);
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.logpdf = logpdf;
  DISTR.pdf = _unur_distr_cvec_eval_pdf_from_logpdf;
  return UNUR_SUCCESS;
} 
double
_unur_distr_cvec_eval_pdf_from_logpdf( const double *x, struct unur_distr *distr )
{
  if (DISTR.logpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return exp(_unur_cvec_logPDF(x,distr));
} 
int
unur_distr_cvec_set_dlogpdf( struct unur_distr *distr, UNUR_VFUNCT_CVEC *dlogpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dlogpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.dpdf != NULL || DISTR.dlogpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dlogPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.dlogpdf = dlogpdf;
  DISTR.dpdf = _unur_distr_cvec_eval_dpdf_from_dlogpdf;
  return UNUR_SUCCESS;
} 
int
_unur_distr_cvec_eval_dpdf_from_dlogpdf( double *result, const double *x, struct unur_distr *distr )
{
  int ret, i;
  double fx;
  if (DISTR.logpdf == NULL || DISTR.dlogpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  fx = exp(unur_distr_cvec_eval_logpdf( x, distr ));
  if (!_unur_isfinite(fx)) return UNUR_ERR_DISTR_DATA;
  ret = _unur_cvec_dlogPDF(result,x,distr);
  for (i=0; i<distr->dim; i++)
    result[i] *= fx;
  return ret;
} 
int
unur_distr_cvec_set_pdlogpdf( struct unur_distr *distr, UNUR_FUNCTD_CVEC *pdlogpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdlogpdf, UNUR_ERR_NULL);
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdpdf != NULL || DISTR.pdlogpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of pdlogPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.pdlogpdf = pdlogpdf;
  DISTR.pdpdf = _unur_distr_cvec_eval_pdpdf_from_pdlogpdf;
  return UNUR_SUCCESS;
} 
double
_unur_distr_cvec_eval_pdpdf_from_pdlogpdf( const double *x, int coord, struct unur_distr *distr )
{
  double fx;
  if (DISTR.logpdf == NULL || DISTR.pdlogpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  if (coord < 0 || coord >= distr->dim) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DOMAIN,"invalid coordinate");
    return INFINITY;
  }
  fx = exp(unur_distr_cvec_eval_logpdf( x, distr ));
  if (!_unur_isfinite(fx)) return INFINITY;
  return fx * _unur_cvec_pdlogPDF(x,coord,distr);
} 
UNUR_FUNCT_CVEC *
unur_distr_cvec_get_logpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  return DISTR.logpdf;
} 
UNUR_VFUNCT_CVEC *
unur_distr_cvec_get_dlogpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  return DISTR.dlogpdf;
} 
UNUR_FUNCTD_CVEC *
unur_distr_cvec_get_pdlogpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  return DISTR.pdlogpdf;
} 
double
unur_distr_cvec_eval_logpdf( const double *x, struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );
  if (DISTR.logpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cvec_logPDF(x,distr);
} 
int
unur_distr_cvec_eval_dlogpdf( double *result, const double *x, struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.dlogpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  return _unur_cvec_dlogPDF(result,x,distr);
} 
double
unur_distr_cvec_eval_pdlogpdf( const double *x, int coord, struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );
  if (DISTR.pdlogpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  if (coord < 0 || coord >= distr->dim) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DOMAIN,"invalid coordinate");
    return INFINITY;
  }
  return _unur_cvec_pdlogPDF(x,coord,distr);
} 
int
unur_distr_cvec_set_domain_rect( struct unur_distr *distr, const double *lowerleft, const double *upperright )
{
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, lowerleft, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, upperright, UNUR_ERR_NULL );
  for (i=0; i<distr->dim; i++) {
    if (!(lowerleft[i] < upperright[i] * (1.-UNUR_SQRT_DBL_EPSILON))) {
      _unur_error(distr->name,UNUR_ERR_DISTR_SET,"domain, left >= right");
      return UNUR_ERR_DISTR_SET;
    }
  }
  DISTR.domainrect = _unur_xrealloc(DISTR.domainrect, 2 * distr->dim * sizeof(double));
  for (i=0; i<distr->dim; i++) {
    DISTR.domainrect[2*i] = lowerleft[i];
    DISTR.domainrect[2*i+1] = upperright[i];
  }
  distr->set |= UNUR_DISTR_SET_DOMAIN | UNUR_DISTR_SET_DOMAINBOUNDED;
  distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_MASK_DERIVED );
  if (distr->base) {
    distr->base->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
			  UNUR_DISTR_SET_MASK_DERIVED );
    if ( distr->base->type == UNUR_DISTR_CVEC ) {
      if (unur_distr_cvec_set_domain_rect(distr->base, lowerleft, upperright)!=UNUR_SUCCESS)
	return UNUR_ERR_DISTR_SET;
    }
  }
  return UNUR_SUCCESS;
} 
int
_unur_distr_cvec_has_boundeddomain( const struct unur_distr *distr )
{
  int i;
  double *domain;
  CHECK_NULL( distr, FALSE );
  COOKIE_CHECK(distr,CK_DISTR_CVEC,FALSE);
  if (! (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED && 
	 DISTR.domainrect))
    return FALSE;
  domain = DISTR.domainrect;
  for (i=0; i < 2*distr->dim; i++) 
    if (!_unur_isfinite(domain[i]))
      return FALSE;
  return TRUE;
} 
int
_unur_distr_cvec_is_indomain( const double *x, const struct unur_distr *distr )
{
  int i;
  double *domain;
  CHECK_NULL( distr, FALSE );
  COOKIE_CHECK(distr,CK_DISTR_CVEC,FALSE);
  domain = DISTR.domainrect;
  if (domain==NULL) 
    return TRUE;
  for (i=0; i<distr->dim; i++) {
    if (x[i] < domain[2*i] || x[i] > domain[2*i+1]) 
      return FALSE;
  }
  return TRUE;
} 
int
unur_distr_cvec_is_indomain( const double *x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  _unur_check_distr_object( distr, CVEC, FALSE );
  return _unur_distr_cvec_is_indomain(x, distr);
}  
int
unur_distr_cvec_set_mean( struct unur_distr *distr, const double *mean )
{
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.mean == NULL)
    DISTR.mean = _unur_xmalloc( distr->dim * sizeof(double) );
  if (mean)
    memcpy( DISTR.mean, mean, distr->dim * sizeof(double) );
  else  
    for (i=0; i<distr->dim; i++)
      DISTR.mean[i] = 0.;
  distr->set |= UNUR_DISTR_SET_MEAN;
  return UNUR_SUCCESS;
} 
const double *
unur_distr_cvec_get_mean( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( !(distr->set & UNUR_DISTR_SET_MEAN) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"mean");
    return NULL;
  }
  return DISTR.mean;
} 
int
unur_distr_cvec_set_covar( struct unur_distr *distr, const double *covar )
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;
  int dim;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  dim = distr->dim;
  distr->set &= ~( UNUR_DISTR_SET_COVAR 
		   | UNUR_DISTR_SET_COVAR_IDENT
		   | UNUR_DISTR_SET_CHOLESKY
		   | UNUR_DISTR_SET_COVAR_INV );
  if (DISTR.covar == NULL)
    DISTR.covar = _unur_xmalloc( dim * dim * sizeof(double) );
  if (DISTR.cholesky == NULL)
    DISTR.cholesky = _unur_xmalloc( dim * dim * sizeof(double) );   
  if (covar==NULL) { 
    for (i=0; i<dim; i++) { 
      for (j=0; j<dim; j++) {
         DISTR.covar[idx(i,j)] = (i==j) ? 1. : 0.;
         DISTR.cholesky[idx(i,j)] = (i==j) ? 1. : 0.;
      } 
    } 
    distr->set |= UNUR_DISTR_SET_COVAR_IDENT;
  } 
  else {
    for (i=0; i<dim*dim; i+= dim+1)
      if (covar[i] <= 0.) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"variance <= 0");
	return UNUR_ERR_DISTR_DOMAIN;
      }
    for (i=0; i<dim; i++)
      for (j=i+1; j<dim; j++)
	if (!_unur_FP_same(covar[i*dim+j],covar[j*dim+i])) {
	  _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,
	              "covariance matrix not symmetric");
	  return UNUR_ERR_DISTR_DOMAIN;
	}
    memcpy( DISTR.covar, covar, dim * dim * sizeof(double) );
    if (_unur_matrix_cholesky_decomposition(dim, covar, DISTR.cholesky) != UNUR_SUCCESS) {
      _unur_error(distr->name, UNUR_ERR_DISTR_DOMAIN, 
		  "covariance matrix not positive definite");
      return UNUR_ERR_DISTR_DOMAIN;      
    }
  }
  distr->set |= UNUR_DISTR_SET_COVAR | UNUR_DISTR_SET_CHOLESKY;
  return UNUR_SUCCESS;
#undef idx
} 
int
unur_distr_cvec_set_covar_inv( struct unur_distr *distr, const double *covar_inv )
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;
  int dim;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  dim = distr->dim;
  distr->set &= ~(UNUR_DISTR_SET_COVAR_INV);
  if (DISTR.covar_inv == NULL)
    DISTR.covar_inv = _unur_xmalloc( dim * dim * sizeof(double) );
  if (covar_inv==NULL)
    for (i=0; i<dim; i++)
      for (j=0; j<dim; j++)
         DISTR.covar_inv[idx(i,j)] = (i==j) ? 1. : 0.;
  else {
    for (i=0; i<dim*dim; i+= dim+1)
      if (covar_inv[i] <= 0.) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"diagonals <= 0");
	return UNUR_ERR_DISTR_DOMAIN;
      }
    for (i=0; i<dim; i++)
      for (j=i+1; j<dim; j++)
	if (!_unur_FP_same(covar_inv[i*dim+j],covar_inv[j*dim+i])) {
	  _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,
	              "inverse of covariance matrix not symmetric");
	  return UNUR_ERR_DISTR_DOMAIN;
	}
    memcpy( DISTR.covar_inv, covar_inv, dim * dim * sizeof(double) );
  }
  distr->set |= UNUR_DISTR_SET_COVAR_INV;
  return UNUR_SUCCESS;
#undef idx
} 
const double *
unur_distr_cvec_get_covar( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( !(distr->set & UNUR_DISTR_SET_COVAR) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"covariance matrix");
    return NULL;
  }
  return DISTR.covar;
} 
const double *
unur_distr_cvec_get_cholesky( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( !(distr->set & UNUR_DISTR_SET_CHOLESKY) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"covariance matrix");
    return NULL;
  }
  return DISTR.cholesky;
} 
const double *
unur_distr_cvec_get_covar_inv ( struct unur_distr *distr )
{
  double det; 
  int dim;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  dim = distr->dim;
  if ( !(distr->set & UNUR_DISTR_SET_COVAR) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"covariance matrix not known");
    return NULL;
  }
  if (DISTR.covar_inv == NULL)
    DISTR.covar_inv = _unur_xmalloc( dim * dim * sizeof(double) );   
  if ( !(distr->set & UNUR_DISTR_SET_COVAR_INV) ) {       
      if (_unur_matrix_invert_matrix(dim, DISTR.covar, DISTR.covar_inv, &det) != UNUR_SUCCESS) {
        _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"cannot compute inverse of covariance");
        return NULL;
      }
  }
  distr->set |= UNUR_DISTR_SET_COVAR_INV;
  return DISTR.covar_inv;
} 
int
unur_distr_cvec_set_rankcorr( struct unur_distr *distr, const double *rankcorr )
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;
  int dim;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  dim = distr->dim;
  distr->set &= ~(UNUR_DISTR_SET_RANKCORR | UNUR_DISTR_SET_RK_CHOLESKY);
  if (DISTR.rankcorr == NULL)
    DISTR.rankcorr = _unur_xmalloc( dim * dim * sizeof(double) );
  if (DISTR.rk_cholesky == NULL)
    DISTR.rk_cholesky = _unur_xmalloc( dim * dim * sizeof(double) );   
  if (rankcorr==NULL) { 
    for (i=0; i<dim; i++)
      for (j=0; j<dim; j++) {
	DISTR.rankcorr[idx(i,j)] = (i==j) ? 1. : 0.;
	DISTR.rk_cholesky[idx(i,j)] = (i==j) ? 1. : 0.;
      }
  } 
  else {
    for (i=0; i<dim*dim; i+= dim+1) {
      if (!_unur_FP_same(rankcorr[i],1.)) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"diagonals != 1");
	return UNUR_ERR_DISTR_DOMAIN;
      }
    }
    for (i=0; i<dim; i++)
      for (j=i+1; j<dim; j++)
	if (!_unur_FP_same(rankcorr[i*dim+j],rankcorr[j*dim+i])) {
	  _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,
	              "rank-correlation matrix not symmetric");
	  return UNUR_ERR_DISTR_DOMAIN;
	}
    memcpy( DISTR.rankcorr, rankcorr, dim * dim * sizeof(double) );
    if (_unur_matrix_cholesky_decomposition(dim, rankcorr, DISTR.rk_cholesky) != UNUR_SUCCESS) {
      _unur_error(distr->name, UNUR_ERR_DISTR_DOMAIN, 
		  "rankcorriance matrix not positive definite");
      return UNUR_ERR_DISTR_DOMAIN;      
    }
  }
  distr->set |= UNUR_DISTR_SET_RANKCORR | UNUR_DISTR_SET_RK_CHOLESKY;
  return UNUR_SUCCESS;
#undef idx
} 
const double *
unur_distr_cvec_get_rankcorr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( !(distr->set & UNUR_DISTR_SET_RANKCORR) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"rank-correlation matrix");
    return NULL;
  }
  return DISTR.rankcorr;
} 
const double *
unur_distr_cvec_get_rk_cholesky( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( !(distr->set & UNUR_DISTR_SET_RK_CHOLESKY) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"rank correlation matrix");
    return NULL;
  }
  return DISTR.rk_cholesky;
} 
int
unur_distr_cvec_set_marginals ( struct unur_distr *distr, struct unur_distr *marginal)
{
  struct unur_distr *clone;
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, marginal, UNUR_ERR_NULL );
  _unur_check_distr_object( marginal, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.marginals)
    _unur_distr_cvec_marginals_free(DISTR.marginals, distr->dim);
  clone = _unur_distr_clone( marginal );
  DISTR.marginals = _unur_xmalloc (distr->dim * sizeof(struct unur_distr *));
  for (i=0; i<distr->dim; i++)
    DISTR.marginals[i] = clone;
  distr->set |= UNUR_DISTR_SET_MARGINAL;
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_set_marginal_array ( struct unur_distr *distr, struct unur_distr **marginals)
{
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, marginals, UNUR_ERR_NULL );
  for (i=0; i<distr->dim; i++) {
    _unur_check_NULL( distr->name, *(marginals+i), UNUR_ERR_NULL );
    _unur_check_distr_object( *(marginals+i), CONT, UNUR_ERR_DISTR_INVALID );
  }
  if (DISTR.marginals)
    _unur_distr_cvec_marginals_free(DISTR.marginals, distr->dim);
  DISTR.marginals = _unur_xmalloc (distr->dim * sizeof(struct unur_distr *));
  for (i=0; i<distr->dim; i++) 
    DISTR.marginals[i] = _unur_distr_clone( *(marginals+i) );
  distr->set |= UNUR_DISTR_SET_MARGINAL;
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_set_marginal_list ( struct unur_distr *distr, ... )
{
  int i;
  int failed = FALSE;
  struct unur_distr *marginal;
  struct unur_distr **marginal_list;
  va_list vargs;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  marginal_list = _unur_xmalloc (distr->dim * sizeof(struct unur_distr *));
  for (i=0; i<distr->dim; i++) marginal_list[i] = NULL;
  va_start(vargs, distr);
  for (i=0; i<distr->dim; i++) {
    marginal = (struct unur_distr *) va_arg(vargs, struct unur_distr *);
    if (marginal) {
      marginal_list[i] = _unur_distr_clone( marginal );
      _unur_distr_free(marginal);
    }
    else {
      failed = TRUE;
    }
  }
  va_end(vargs);
  if (failed) {
    _unur_distr_cvec_marginals_free(marginal_list, distr->dim);
    _unur_error(distr->name ,UNUR_ERR_DISTR_SET,"marginals == NULL");
    return UNUR_ERR_DISTR_SET;
  }
  if (DISTR.marginals)
    _unur_distr_cvec_marginals_free(DISTR.marginals, distr->dim);
  DISTR.marginals = marginal_list;
  distr->set |= UNUR_DISTR_SET_MARGINAL;
  return UNUR_SUCCESS;
} 
const struct unur_distr *
unur_distr_cvec_get_marginal( const struct unur_distr *distr, int n )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if (n<=0 || n > distr->dim) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"n not in 1 .. dim");
    return NULL;
  }
  if ( !(distr->set & UNUR_DISTR_SET_MARGINAL) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"marginals");
    return NULL;
  }
  _unur_check_NULL( distr->name, DISTR.marginals, NULL );
  return (DISTR.marginals[n-1]);
} 
struct unur_distr **
_unur_distr_cvec_marginals_clone ( struct unur_distr **marginals, int dim )
{
  struct unur_distr **clone;
  int i;
  _unur_check_NULL( NULL, marginals, NULL );
  if (dim < 1) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"dimension < 1");
    return NULL;
  }
  clone = _unur_xmalloc (dim * sizeof(struct unur_distr *));
  if (_unur_distr_cvec_marginals_are_equal(marginals, dim)) {
      clone[0] = _unur_distr_clone( marginals[0] );
      for (i=1; i<dim; i++)
	clone[i] = clone[0];
  }
  else {
    for (i=0; i<dim; i++) 
      clone[i] = _unur_distr_clone( marginals[i] );
  }
  return clone;
} 
void
_unur_distr_cvec_marginals_free ( struct unur_distr **marginals, int dim )
{
  int i;
  if (_unur_distr_cvec_marginals_are_equal(marginals,dim)) {
    _unur_distr_free(marginals[0]);
  }
  else {
    for (i=0; i<dim; i++) 
      if (marginals[i]) _unur_distr_free(marginals[i]);
  }
  free (marginals);
} 
int 
_unur_distr_cvec_marginals_are_equal( struct unur_distr **marginals, int dim )
{
  return (dim <= 1 || marginals[0] == marginals[1]) ? TRUE : FALSE;
} 
int
_unur_distr_cvec_duplicate_firstmarginal( struct unur_distr *distr )
{
  struct unur_distr *marginal;
  int i;
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  marginal = DISTR.marginals[0];
  if ( !(distr->set & UNUR_DISTR_SET_MARGINAL) || marginal==NULL ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"marginals");
    return UNUR_ERR_DISTR_DATA;
  }
  if (!_unur_distr_cvec_marginals_are_equal(DISTR.marginals,distr->dim)) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"marginals not equal");
    return UNUR_ERR_DISTR_DATA;
  }
  for (i=1; i<distr->dim; i++) 
    DISTR.marginals[i] = _unur_distr_clone( marginal );
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_set_pdfparams( struct unur_distr *distr, const double *params, int n_params )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( NULL, params, UNUR_ERR_NULL );  
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name,params,UNUR_ERR_NULL);
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.n_params = n_params;
  if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_get_pdfparams( const struct unur_distr *distr, const double **params )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );
  *params = (DISTR.n_params) ? DISTR.params : NULL;
  return DISTR.n_params;
} 
int
unur_distr_cvec_set_pdfparams_vec( struct unur_distr *distr, int par, const double *param_vec, int n_param_vec )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }
  if (param_vec != NULL) {
    DISTR.param_vecs[par] = _unur_xrealloc( DISTR.param_vecs[par], n_param_vec * sizeof(double) );
    memcpy( DISTR.param_vecs[par], param_vec, n_param_vec*sizeof(double) );
    DISTR.n_param_vec[par] = n_param_vec;
  }
  else {
    if (DISTR.param_vecs[par]) free(DISTR.param_vecs[par]);
    DISTR.n_param_vec[par] = 0;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  return UNUR_SUCCESS;
} 
int
unur_distr_cvec_get_pdfparams_vec( const struct unur_distr *distr, int par, const double **param_vecs )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_NPARAMS,"");
    *param_vecs = NULL;
    return 0;
  }
  *param_vecs = DISTR.param_vecs[par];
  return (*param_vecs) ? DISTR.n_param_vec[par] : 0;
} 
int
unur_distr_cvec_set_mode( struct unur_distr *distr, const double *mode )
{
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.mode == NULL)
    DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  if (mode)
    memcpy( DISTR.mode, mode, distr->dim * sizeof(double) );
  else  
    for (i=0; i<distr->dim; i++)
      DISTR.mode[i] = 0.;
  distr->set |= UNUR_DISTR_SET_MODE;
  return UNUR_SUCCESS;
} 
int 
unur_distr_cvec_upd_mode( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.upd_mode == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  if ((DISTR.upd_mode)(distr)==UNUR_SUCCESS) {
    distr->set |= UNUR_DISTR_SET_MODE;
    return UNUR_SUCCESS;
  }
  else {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
} 
const double *
unur_distr_cvec_get_mode( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    if (DISTR.upd_mode == NULL) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return NULL;
    }
    else {
      if (unur_distr_cvec_upd_mode(distr)!=UNUR_SUCCESS) {
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return NULL;
      }
    }
  }
  return DISTR.mode;
} 
int
unur_distr_cvec_set_center( struct unur_distr *distr, const double *center )
{
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.center == NULL)
    DISTR.center = _unur_xmalloc( distr->dim * sizeof(double) );
  if (center)
    memcpy( DISTR.center, center, distr->dim * sizeof(double) );
  else  
    for (i=0; i<distr->dim; i++)
      DISTR.center[i] = 0.;
  distr->set |= UNUR_DISTR_SET_CENTER;
  return UNUR_SUCCESS;
} 
const double *
unur_distr_cvec_get_center( struct unur_distr *distr )
{
  int i;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );
  if ( distr->set & UNUR_DISTR_SET_CENTER )
    return DISTR.center;
  if ( distr->set & UNUR_DISTR_SET_MODE ) 
    return DISTR.mode;
  if ( distr->set & UNUR_DISTR_SET_MEAN ) 
    return DISTR.mean;
  if ( DISTR.center == NULL )
    DISTR.center = _unur_xmalloc( distr->dim * sizeof(double) );
  for (i=0; i<distr->dim; i++) 
    DISTR.center[i] = 0.;
  return DISTR.center;
} 
int
unur_distr_cvec_set_pdfvol( struct unur_distr *distr, double volume )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (volume <= 0.) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PDF volume <= 0");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.volume = volume;
  distr->set |= UNUR_DISTR_SET_PDFVOLUME;
  return UNUR_SUCCESS;
} 
int 
unur_distr_cvec_upd_pdfvol( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  if (DISTR.upd_volume == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  if (((DISTR.upd_volume)(distr)!=UNUR_SUCCESS) || DISTR.volume <= 0.) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"upd volume <= 0");
    DISTR.volume = 1.;   
    distr->set &= ~UNUR_DISTR_SET_PDFVOLUME;
    return UNUR_ERR_DISTR_SET;
  }
  distr->set |= UNUR_DISTR_SET_PDFVOLUME;
  return UNUR_SUCCESS;
} 
double
unur_distr_cvec_get_pdfvol( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );
  if ( !(distr->set & UNUR_DISTR_SET_PDFVOLUME) ) {
    if (DISTR.upd_volume == NULL) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"volume");
      return INFINITY;
    }
    else {
      unur_distr_cvec_upd_pdfvol( distr );
    }
  }
  return DISTR.volume;
} 
double
_unur_cvec_PDF(const double *x, struct unur_distr *distr)
{
  if ( (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
       !_unur_distr_cvec_is_indomain(x, distr) )
    return 0.;
  return (*(distr->data.cvec.pdf)) (x,distr);
}
int
_unur_cvec_dPDF(double *result, const double *x, struct unur_distr *distr)
{
  if ( (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
       !_unur_distr_cvec_is_indomain(x, distr) )
    return 0.;
  return (*(distr->data.cvec.dpdf)) (result,x,distr);
}
double
_unur_cvec_pdPDF(const double *x, int coord, struct unur_distr *distr)
{
  if ( (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
       !_unur_distr_cvec_is_indomain(x, distr) )
    return 0.;
  return (*(distr->data.cvec.pdpdf)) (x,coord,distr);
}
double
_unur_cvec_logPDF(const double *x, struct unur_distr *distr)
{
  if ( (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
       !_unur_distr_cvec_is_indomain(x, distr) )
    return -INFINITY;
  return (*(distr->data.cvec.logpdf)) (x,distr);
}
int
_unur_cvec_dlogPDF(double *result, const double *x, struct unur_distr *distr)
{
  if ( (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
       !_unur_distr_cvec_is_indomain(x, distr) )
    return 0.;
  return (*(distr->data.cvec.dlogpdf)) (result,x,distr);
}
double
_unur_cvec_pdlogPDF(const double *x, int coord, struct unur_distr *distr)
{
  if ( (distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
       !_unur_distr_cvec_is_indomain(x, distr) )
    return 0.;
  return (*(distr->data.cvec.pdlogpdf)) (x,coord,distr);
}
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_cvec_debug( const struct unur_distr *distr, const char *genid )
{
  FILE *log;
  double *mat;
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CVEC,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous multivariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);
  fprintf(log,"%s:\tdimension = %d\n",genid,distr->dim);
  fprintf(log,"%s:\tfunctions: ",genid);
  if (DISTR.pdf) fprintf(log,"PDF ");
  if (DISTR.logpdf) fprintf(log,"logPDF ");
  if (DISTR.dpdf) fprintf(log,"dPDF ");
  if (DISTR.dlogpdf) fprintf(log,"dlogPDF ");
  if (DISTR.pdpdf) fprintf(log,"pdPDF ");
  if (DISTR.pdlogpdf) fprintf(log,"pdlogPDF ");
  fprintf(log,"\n%s:\n",genid);
  fprintf(log,"%s:\tdomain = ",genid);
  if (!(distr->set & UNUR_DISTR_SET_DOMAINBOUNDED)) {
    fprintf(log,"unbounded\n");
  }
  else {
    if (DISTR.domainrect) {
      double *domain = DISTR.domainrect;
      int i;
      fprintf(log,"rectangular\n");
      for (i=0; i<distr->dim; i++)
	fprintf(log,"%s:\t %c ( %g, %g)\n",genid, i?'x':' ', 
		domain[2*i], domain[2*i+1]);
    }
  }
  fprintf(log,"%s:\n",genid);
  mat = ((distr->set & UNUR_DISTR_SET_MODE) && DISTR.mode) ? DISTR.mode : NULL;
  _unur_matrix_print_vector( distr->dim, mat, "\tmode =", log, genid, "\t   ");
  mat = ((distr->set & UNUR_DISTR_SET_MEAN) && DISTR.mean) ? DISTR.mean : NULL;
  _unur_matrix_print_vector( distr->dim, mat, "\tmean vector =", log, genid, "\t   ");
  if ((distr->set & UNUR_DISTR_SET_CENTER) && DISTR.center)
    _unur_matrix_print_vector( distr->dim, DISTR.center, "\tcenter vector =", log, genid, "\t   ");
  else {
    fprintf(log,"%s:\tcenter = mode [not given explicitly]\n",genid);
    fprintf(log,"%s:\n",genid);
  }
  mat = ((distr->set & UNUR_DISTR_SET_COVAR) && DISTR.covar) ? DISTR.covar : NULL;
  _unur_matrix_print_matrix( distr->dim, mat, "\tcovariance matrix =", log, genid, "\t   ");
  mat = ((distr->set & UNUR_DISTR_SET_CHOLESKY) && DISTR.cholesky) ? DISTR.cholesky : NULL;
  _unur_matrix_print_matrix( distr->dim, mat, "\tcholesky factor of covariance matrix =", log, genid, "\t   ");
  mat = ((distr->set & UNUR_DISTR_SET_RANKCORR) && DISTR.rankcorr) ? DISTR.rankcorr : NULL;
  _unur_matrix_print_matrix( distr->dim, mat, "\trank correlation matrix =", log, genid, "\t   ");
  mat = ((distr->set & UNUR_DISTR_SET_RK_CHOLESKY) && DISTR.rk_cholesky) ? DISTR.rk_cholesky : NULL;
  _unur_matrix_print_matrix( distr->dim, mat, "\tcholesky factor of rank correlation matrix =", log, genid, "\t   ");
  fprintf(log,"%s:\tmarginal distributions:\n",genid);
  if (distr->set & UNUR_DISTR_SET_MARGINAL) {
    if (_unur_distr_cvec_marginals_are_equal(DISTR.marginals, distr->dim)) {
      fprintf(log,"%s: all mariginals [1-%d]:\n",genid,distr->dim);
      _unur_distr_cont_debug( DISTR.marginals[0], genid );
    }
    else {
      int i;
      for (i=0; i<distr->dim; i++) {
	fprintf(log,"%s: mariginal [%d]:\n",genid,i+1);
	_unur_distr_cont_debug( DISTR.marginals[i], genid );
      }
    }
  }
  else {
    fprintf(log,"%s:\t   [unknown]\n",genid);
  }
  fprintf(log,"%s:\n",genid);
#ifdef USE_DEPRECATED_CODE
  fprintf(log,"%s:\tstandardized marginal distributions:   [see also marginal distributions]\n",genid);
  if (distr->set & UNUR_DISTR_SET_STDMARGINAL) {
    if (_unur_distr_cvec_marginals_are_equal(DISTR.stdmarginals, distr->dim)) {
      fprintf(log,"%s: all standardized mariginals [1-%d]:\n",genid,distr->dim);
      _unur_distr_cont_debug( DISTR.stdmarginals[0], genid );
    }
    else {
      int i;
      for (i=0; i<distr->dim; i++) {
	fprintf(log,"%s: mariginal [%d]:\n",genid,i+1);
	_unur_distr_cont_debug( DISTR.stdmarginals[i], genid );
      }
    }
  }
  else {
    fprintf(log,"%s:\t   [unknown]\n",genid);
  }
  fprintf(log,"%s:\n",genid);
#endif
} 
#endif    
#undef DISTR
