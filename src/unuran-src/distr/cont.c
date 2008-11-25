/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include <parser/functparser_source.h>
#include <utils/fmax_source.h>
#include "distr_source.h"
#include "distr.h"
#include "cont.h"
#define DISTR distr->data.cont
#define BASE  distr->base->data.cont
static double _unur_distr_cont_eval_pdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_logpdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_dpdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_dlogpdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_cdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_logcdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_hr_tree( double x, const struct unur_distr *distr );
static void _unur_distr_cont_free( struct unur_distr *distr );
static int _unur_distr_cont_find_mode( struct unur_distr *distr );
static double _unur_aux_pdf(double x, void *p);
struct unur_distr *
unur_distr_cont_new( void )
{
  register struct unur_distr *distr;
  int i;
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;
  COOKIE_SET(distr,CK_DISTR_CONT);
  distr->type = UNUR_DISTR_CONT;
  distr->id = UNUR_DISTR_GENERIC;
  distr->dim = 1;   
  distr->destroy = _unur_distr_cont_free;
  distr->clone = _unur_distr_cont_clone;
  DISTR.pdf       = NULL;          
  DISTR.dpdf      = NULL;          
  DISTR.logpdf    = NULL;          
  DISTR.dlogpdf   = NULL;          
  DISTR.cdf       = NULL;          
  DISTR.logcdf    = NULL;          
  DISTR.hr        = NULL;          
  DISTR.init      = NULL;          
  DISTR.n_params  = 0;               
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    DISTR.n_param_vec[i] = 0;
    DISTR.param_vecs[i] = NULL;
  }  
  DISTR.norm_constant = 1.;        
  DISTR.mode       = INFINITY;     
  DISTR.center     = 0.;           
  DISTR.area       = 1.;           
  DISTR.trunc[0] = DISTR.domain[0] = -INFINITY; 
  DISTR.trunc[1] = DISTR.domain[1] = INFINITY;  
  DISTR.set_params = NULL;         
  DISTR.upd_mode   = _unur_distr_cont_find_mode; 
  DISTR.upd_area   = NULL;         
  DISTR.pdftree    = NULL;         
  DISTR.dpdftree   = NULL;         
  DISTR.logpdftree = NULL;         
  DISTR.dlogpdftree = NULL;        
  DISTR.cdftree    = NULL;         
  DISTR.logcdftree = NULL;         
  DISTR.hrtree     = NULL;         
  return distr;
} 
struct unur_distr *
_unur_distr_cont_clone( const struct unur_distr *distr )
{
#define CLONE clone->data.cont
  struct unur_distr *clone;
  int i;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  memcpy( clone, distr, sizeof( struct unur_distr ) );
  CLONE.pdftree  = (DISTR.pdftree)  ? _unur_fstr_dup_tree(DISTR.pdftree)  : NULL;
  CLONE.dpdftree = (DISTR.dpdftree) ? _unur_fstr_dup_tree(DISTR.dpdftree) : NULL;
  CLONE.logpdftree  = (DISTR.logpdftree)  ? _unur_fstr_dup_tree(DISTR.logpdftree)  : NULL;
  CLONE.dlogpdftree = (DISTR.dlogpdftree) ? _unur_fstr_dup_tree(DISTR.dlogpdftree) : NULL;
  CLONE.cdftree  = (DISTR.cdftree)  ? _unur_fstr_dup_tree(DISTR.cdftree)  : NULL;
  CLONE.logcdftree  = (DISTR.logcdftree)  ? _unur_fstr_dup_tree(DISTR.logcdftree)  : NULL;
  CLONE.hrtree   = (DISTR.hrtree)   ? _unur_fstr_dup_tree(DISTR.hrtree)   : NULL;
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
  if (distr->base != NULL) {
    clone->base = _unur_distr_clone(distr->base);
  }
  return clone;
#undef CLONE
} 
void
_unur_distr_cont_free( struct unur_distr *distr )
{
  int i;
  if( distr == NULL ) 
    return;
  _unur_check_distr_object( distr, CONT, RETURN_VOID );
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    if (DISTR.param_vecs[i]) free( DISTR.param_vecs[i] );
  if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
  if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
  if (DISTR.logpdftree)  _unur_fstr_free(DISTR.logpdftree);
  if (DISTR.dlogpdftree) _unur_fstr_free(DISTR.dlogpdftree);
  if (DISTR.cdftree)  _unur_fstr_free(DISTR.cdftree);
  if (DISTR.logcdftree)  _unur_fstr_free(DISTR.logcdftree);
  if (DISTR.hrtree)   _unur_fstr_free(DISTR.hrtree);
  if (distr->base) _unur_distr_free(distr->base);
  if (distr->name_str) free(distr->name_str);
  COOKIE_CLEAR(distr);
  free( distr );
} 
int
unur_distr_cont_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CONT *pdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.pdf = pdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_dpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *dpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.dpdf != NULL || DISTR.dlogpdf != NULL ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.dpdf = dpdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_logpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *logpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, logpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.logpdf = logpdf;
  DISTR.pdf = _unur_distr_cont_eval_pdf_from_logpdf;
  return UNUR_SUCCESS;
} 
double
_unur_distr_cont_eval_pdf_from_logpdf( double x, const struct unur_distr *distr )
{
  if (DISTR.logpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return exp(_unur_cont_logPDF(x,distr));
} 
int
unur_distr_cont_set_dlogpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *dlogpdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dlogpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.dpdf != NULL || DISTR.dlogpdf != NULL ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dlogPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.dlogpdf = dlogpdf;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_from_dlogpdf;
  return UNUR_SUCCESS;
} 
double
_unur_distr_cont_eval_dpdf_from_dlogpdf( double x, const struct unur_distr *distr )
{
  if (DISTR.logpdf == NULL || DISTR.dlogpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return exp(_unur_cont_logPDF(x,distr)) * _unur_cont_dlogPDF(x,distr);
} 
int
unur_distr_cont_set_cdf( struct unur_distr *distr, UNUR_FUNCT_CONT *cdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, cdf,UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.cdf != NULL || DISTR.logcdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.cdf = cdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_logcdf( struct unur_distr *distr, UNUR_FUNCT_CONT *logcdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, logcdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.cdf != NULL || DISTR.logcdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logCDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.logcdf = logcdf;
  DISTR.cdf = _unur_distr_cont_eval_cdf_from_logcdf;
  return UNUR_SUCCESS;
} 
double
_unur_distr_cont_eval_cdf_from_logcdf( double x, const struct unur_distr *distr )
{
  if (DISTR.logcdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return exp(_unur_cont_logCDF(x,distr));
} 
int
unur_distr_cont_set_hr( struct unur_distr *distr, UNUR_FUNCT_CONT *hr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, hr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.hr != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of HR not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.hr = hr;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_pdfstr( struct unur_distr *distr, const char *pdfstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, pdfstr, UNUR_ERR_NULL );
  if ( DISTR.pdftree || DISTR.logpdftree ) {
    if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
    if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
    if (DISTR.logpdftree)  _unur_fstr_free(DISTR.logpdftree);
    if (DISTR.dlogpdftree) _unur_fstr_free(DISTR.dlogpdftree);
    DISTR.pdf = NULL;
    DISTR.dpdf = NULL;
    DISTR.logpdf = NULL;
    DISTR.dlogpdf = NULL;
  }
  if (DISTR.pdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.pdftree = _unur_fstr2tree(pdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.pdf  = _unur_distr_cont_eval_pdf_tree;
  if ( (DISTR.dpdftree = _unur_fstr_make_derivative(DISTR.pdftree)) == NULL )
    return UNUR_ERR_DISTR_DATA;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_tree;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_logpdfstr( struct unur_distr *distr, const char *logpdfstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, logpdfstr, UNUR_ERR_NULL );
  if ( DISTR.pdftree || DISTR.logpdftree ) {
    if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
    if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
    if (DISTR.logpdftree)  _unur_fstr_free(DISTR.logpdftree);
    if (DISTR.dlogpdftree) _unur_fstr_free(DISTR.dlogpdftree);
    DISTR.pdf = NULL;
    DISTR.dpdf = NULL;
    DISTR.logpdf = NULL;
    DISTR.dlogpdf = NULL;
  }
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.logpdftree = _unur_fstr2tree(logpdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.logpdf  = _unur_distr_cont_eval_logpdf_tree;
  DISTR.pdf = _unur_distr_cont_eval_pdf_from_logpdf;
  if ( (DISTR.dlogpdftree = _unur_fstr_make_derivative(DISTR.logpdftree)) == NULL )
    return UNUR_ERR_DISTR_DATA;
  DISTR.dlogpdf = _unur_distr_cont_eval_dlogpdf_tree;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_from_dlogpdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_cdfstr( struct unur_distr *distr, const char *cdfstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, cdfstr, UNUR_ERR_NULL );
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.cdftree = _unur_fstr2tree(cdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.cdf  = _unur_distr_cont_eval_cdf_tree;
  if (DISTR.pdftree == NULL)
    if ( (DISTR.pdftree = _unur_fstr_make_derivative(DISTR.cdftree)) != NULL )
      DISTR.pdf = _unur_distr_cont_eval_pdf_tree;
  if (DISTR.dpdftree == NULL)
    if ( (DISTR.dpdftree = _unur_fstr_make_derivative(DISTR.pdftree)) != NULL )
      DISTR.dpdf = _unur_distr_cont_eval_dpdf_tree;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_logcdfstr( struct unur_distr *distr, const char *logcdfstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, logcdfstr, UNUR_ERR_NULL );
  if (DISTR.cdf != NULL || DISTR.logcdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logCDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.logcdftree = _unur_fstr2tree(logcdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.logcdf  = _unur_distr_cont_eval_logcdf_tree;
  DISTR.cdf = _unur_distr_cont_eval_cdf_from_logcdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_hrstr( struct unur_distr *distr, const char *hrstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, hrstr, UNUR_ERR_NULL );
  if (DISTR.hr != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.hrtree = _unur_fstr2tree(hrstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.hr  = _unur_distr_cont_eval_hr_tree;
  return UNUR_SUCCESS;
} 
double
_unur_distr_cont_eval_pdf_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.pdftree) ? _unur_fstr_eval_tree(DISTR.pdftree,x) : INFINITY);
} 
double
_unur_distr_cont_eval_logpdf_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.logpdftree) ? _unur_fstr_eval_tree(DISTR.logpdftree,x) : INFINITY);
} 
double
_unur_distr_cont_eval_dpdf_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.dpdftree) ? _unur_fstr_eval_tree(DISTR.dpdftree,x) : INFINITY);
} 
double
_unur_distr_cont_eval_dlogpdf_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.dlogpdftree) ? _unur_fstr_eval_tree(DISTR.dlogpdftree,x) : INFINITY);
} 
double
_unur_distr_cont_eval_cdf_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.cdftree) ? _unur_fstr_eval_tree(DISTR.cdftree,x) : INFINITY);
} 
double
_unur_distr_cont_eval_logcdf_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.logcdftree) ? _unur_fstr_eval_tree(DISTR.logcdftree,x) : INFINITY);
} 
double
_unur_distr_cont_eval_hr_tree( double x, const struct unur_distr *distr )
{
  return ((DISTR.hrtree) ? _unur_fstr_eval_tree(DISTR.hrtree,x) : INFINITY);
} 
char *
unur_distr_cont_get_pdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.pdftree, NULL );
  return _unur_fstr_tree2string(DISTR.pdftree,"x","PDF",TRUE);
} 
char *
unur_distr_cont_get_dpdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.dpdftree, NULL );
  return _unur_fstr_tree2string(DISTR.dpdftree,"x","dPDF",TRUE);
} 
char *
unur_distr_cont_get_logpdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.logpdftree, NULL );
  return _unur_fstr_tree2string(DISTR.logpdftree,"x","logPDF",TRUE);
} 
char *
unur_distr_cont_get_dlogpdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.dlogpdftree, NULL );
  return _unur_fstr_tree2string(DISTR.dlogpdftree,"x","dlogPDF",TRUE);
} 
char *
unur_distr_cont_get_cdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.cdftree, NULL );
  return _unur_fstr_tree2string(DISTR.cdftree,"x","CDF",TRUE);
} 
char *
unur_distr_cont_get_logcdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.logcdftree, NULL );
  return _unur_fstr_tree2string(DISTR.logcdftree,"x","logCDF",TRUE);
} 
char *
unur_distr_cont_get_hrstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.hrtree, NULL );
  return _unur_fstr_tree2string(DISTR.hrtree,"x","HR",TRUE);
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_pdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.pdf;
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_dpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.dpdf;
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_logpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.logpdf;
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_dlogpdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.dlogpdf;
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_cdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.cdf;
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_logcdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.logcdf;
} 
UNUR_FUNCT_CONT *
unur_distr_cont_get_hr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  return DISTR.hr;
} 
double
unur_distr_cont_eval_pdf( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.pdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_PDF(x,distr);
} 
double
unur_distr_cont_eval_dpdf( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.dpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_dPDF(x,distr);
} 
double
unur_distr_cont_eval_logpdf( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.logpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_logPDF(x,distr);
} 
double
unur_distr_cont_eval_dlogpdf( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.dlogpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_dlogPDF(x,distr);
} 
double
unur_distr_cont_eval_cdf( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.cdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_CDF(x,distr);
} 
double
unur_distr_cont_eval_logcdf( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.logcdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_logCDF(x,distr);
} 
double
unur_distr_cont_eval_hr( double x, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if (DISTR.hr == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }
  return _unur_cont_HR(x,distr);
} 
int
unur_distr_cont_set_pdfparams( struct unur_distr *distr, const double *params, int n_params )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name,params,UNUR_ERR_NULL);
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if (distr->base && BASE.set_params) 
    return (BASE.set_params(distr->base,params,n_params));
  if (DISTR.set_params)
    return (DISTR.set_params(distr,params,n_params));
  if (distr->base) {
    BASE.n_params = n_params;
    if (n_params) memcpy( BASE.params, params, n_params*sizeof(double) );
  }
  else {
    DISTR.n_params = n_params;
    if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );
  }
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_get_pdfparams( const struct unur_distr *distr, const double **params )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );
  if (distr->base) {
    *params = (BASE.n_params) ? BASE.params : NULL;
    return BASE.n_params;
  }
  else {
    *params = (DISTR.n_params) ? DISTR.params : NULL;
    return DISTR.n_params;
  }
} 
int
unur_distr_cont_set_pdfparams_vec( struct unur_distr *distr, int par, const double *param_vec, int n_param_vec )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"invalid parameter position");
    return UNUR_ERR_DISTR_NPARAMS;
  }
  if (param_vec != NULL) {
    DISTR.param_vecs[par] = _unur_xrealloc( DISTR.param_vecs[par], n_param_vec * sizeof(double) );
    memcpy( DISTR.param_vecs[par], param_vec, n_param_vec*sizeof(double) );
    DISTR.n_param_vec[par] = n_param_vec;
  }
  else {
    if (DISTR.param_vecs[par]) free(DISTR.param_vecs[par]);
    DISTR.param_vecs[par] = NULL;
    DISTR.n_param_vec[par] = 0;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_get_pdfparams_vec( const struct unur_distr *distr, int par, const double **param_vecs )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"invalid parameter position");
    *param_vecs = NULL;
    return 0;
  }
  *param_vecs = DISTR.param_vecs[par];
  return (*param_vecs) ? DISTR.n_param_vec[par] : 0;
} 
int
unur_distr_cont_set_domain( struct unur_distr *distr, double left, double right )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (left >= right) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  if ( (distr->set & UNUR_DISTR_SET_MODE) &&
       (left  >= DISTR.domain[0]) &&
       (right <= DISTR.domain[1]) ) {
    if ( DISTR.mode < left)       DISTR.mode = left;
    else if ( DISTR.mode > right) DISTR.mode = right;
  }
  DISTR.trunc[0] = DISTR.domain[0] = left;
  DISTR.trunc[1] = DISTR.domain[1] = right;
  distr->set |= UNUR_DISTR_SET_DOMAIN;
  if (distr->set & UNUR_DISTR_SET_MODE) {
    distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		    UNUR_DISTR_SET_TRUNCATED | 
		    UNUR_DISTR_SET_MASK_DERIVED );
    distr->set |= UNUR_DISTR_SET_MODE;
  }
  else {
    distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		    UNUR_DISTR_SET_TRUNCATED | 
		    UNUR_DISTR_SET_MASK_DERIVED );
  }
  if (distr->base) {
    BASE.trunc[0] = BASE.domain[0] = left;
    BASE.trunc[1] = BASE.domain[1] = right;
    distr->base->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
			  UNUR_DISTR_SET_TRUNCATED | 
			  UNUR_DISTR_SET_MASK_DERIVED );
  }
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_get_domain( const struct unur_distr *distr, double *left, double *right )
{
  *left = -INFINITY;
  *right = INFINITY;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_get_truncated( const struct unur_distr *distr, double *left, double *right )
{
  *left = -INFINITY;
  *right = INFINITY;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  *left  = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[0] : DISTR.domain[0];
  *right = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[1] : DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
unur_distr_cont_set_mode( struct unur_distr *distr, double mode )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (mode < DISTR.domain[0] || mode > DISTR.domain[1]) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"mode not in domain");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.mode = mode;
  distr->set |= UNUR_DISTR_SET_MODE;
  return UNUR_SUCCESS;
} 
int 
unur_distr_cont_upd_mode( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
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
double
unur_distr_cont_get_mode( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    if (DISTR.upd_mode == NULL) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INFINITY;
    }
    else {
      if (unur_distr_cont_upd_mode(distr)!=UNUR_SUCCESS) {
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return INFINITY;
      }
    }
  }
  return DISTR.mode;
} 
int
unur_distr_cont_set_center( struct unur_distr *distr, double center )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  DISTR.center = center;
  distr->set |= UNUR_DISTR_SET_CENTER;
  return UNUR_SUCCESS;
} 
double
unur_distr_cont_get_center( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, 0. );
  _unur_check_distr_object( distr, CONT, 0. );
  if ( distr->set & UNUR_DISTR_SET_CENTER )
    return DISTR.center;
  if ( distr->set & UNUR_DISTR_SET_MODE ) 
    return DISTR.mode;
  return 0.;
} 
int
unur_distr_cont_set_pdfarea( struct unur_distr *distr, double area )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (area <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"pdf area <= 0");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.area = area;
  distr->set |= UNUR_DISTR_SET_PDFAREA;
  return UNUR_SUCCESS;
} 
int 
unur_distr_cont_upd_pdfarea( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.upd_area == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  if (((DISTR.upd_area)(distr)!=UNUR_SUCCESS) || DISTR.area <= 0.) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"upd area <= 0");
    DISTR.area = 1.;   
    distr->set &= ~UNUR_DISTR_SET_PDFAREA;
    return UNUR_ERR_DISTR_SET;
  }
  distr->set |= UNUR_DISTR_SET_PDFAREA;
  return UNUR_SUCCESS;
} 
double
unur_distr_cont_get_pdfarea( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );
  if ( !(distr->set & UNUR_DISTR_SET_PDFAREA) ) {
    if ( unur_distr_cont_upd_pdfarea(distr) != UNUR_SUCCESS ) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"area");
      return INFINITY;
    }
  }
  return DISTR.area;
} 
double
_unur_aux_pdf(double x, void *p) 
{
  struct unur_distr *distr = p;
  return (DISTR.pdf(x, distr));
} 
int
_unur_distr_cont_find_mode( struct unur_distr *distr )
{
  struct unur_funct_generic pdf;  
  double mode;
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"PDF required for finding mode numerically"); 
    return UNUR_ERR_DISTR_DATA;
  }
  pdf.f = _unur_aux_pdf;
  pdf.params = distr;
  mode = _unur_util_find_max( pdf, DISTR.domain[0], DISTR.domain[1], DISTR.center );
  if (_unur_isfinite(mode)){
    DISTR.mode = mode;
    distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
    return UNUR_SUCCESS;
  }
  else {
    return UNUR_ERR_DISTR_DATA;
  }
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_cont_debug( const struct unur_distr *distr, const char *genid )
{
  FILE *LOG;
  int i;
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CONT,RETURN_VOID);
  LOG = unur_get_stream();
  if (distr->base) {
    switch (distr->id) {
    case UNUR_DISTR_CORDER:
      _unur_distr_corder_debug(distr,genid);
      return;
    case UNUR_DISTR_CXTRANS:
      _unur_distr_cxtrans_debug(distr,genid);
      return;
    case UNUR_DISTR_CONDI:
      _unur_distr_condi_debug(distr,genid);
      return;
    default:
      fprintf(LOG,"%s: derived distribution.\n",genid);
      fprintf(LOG,"%s:\n",genid);
    }
  }
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous univariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);
  fprintf(LOG,"%s:\tPDF with %d argument(s)\n",genid,DISTR.n_params);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);
  fprintf(LOG,"%s:\tfunctions: ",genid);
  if (DISTR.cdf) fprintf(LOG,"CDF ");
  if (DISTR.logcdf) fprintf(LOG,"logCDF ");
  if (DISTR.pdf) fprintf(LOG,"PDF ");
  if (DISTR.logpdf) fprintf(LOG,"logPDF ");
  if (DISTR.dpdf) fprintf(LOG,"dPDF ");
  if (DISTR.dlogpdf) fprintf(LOG,"dlogPDF ");
  if (DISTR.hr) fprintf(LOG,"HR ");
  fprintf(LOG,"\n");
  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %g\n",genid,DISTR.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);
  if (distr->set & UNUR_DISTR_SET_CENTER)
    fprintf(LOG,"%s:\tcenter = %g\n",genid,DISTR.center);
  else
    fprintf(LOG,"%s:\tcenter = %g [default]\n",genid,
	    unur_distr_cont_get_center(distr));
  fprintf(LOG,"%s:\tdomain = (%g, %g)",genid,DISTR.domain[0],DISTR.domain[1]);
  _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
  fprintf(LOG,"\n%s:\tarea below p.d.f. = %g",genid,DISTR.area);
  _unur_print_if_default(distr,UNUR_DISTR_SET_PDFAREA);
  fprintf(LOG,"\n%s:\n",genid);
} 
#endif    
#undef DISTR
