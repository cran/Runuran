/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include <utils/matrix_source.h>
#include "distr.h"
#include "cont.h"
#include "condi.h"
#include "distr_source.h"
static const char distr_name[] = "conditional";
#define DISTR distr->data.cvec    
#define CONDI condi->data.cont    
#define iK           0            
#define K            (condi->data.cont.params[iK])
#define iPOSITION    0            
#define iDIRECTION   1            
#define iXARG        2            
#define iGRADF       3            
#define POSITION     (condi->data.cont.param_vecs[iPOSITION])
#define DIRECTION    (condi->data.cont.param_vecs[iDIRECTION])
#define XARG         (condi->data.cont.param_vecs[iXARG])
#define GRADF        (condi->data.cont.param_vecs[iGRADF])
static double _unur_pdf_condi( double x, const struct unur_distr *condi );
static double _unur_dpdf_condi( double x, const struct unur_distr *condi );
static double _unur_logpdf_condi( double x, const struct unur_distr *condi );
static double _unur_dlogpdf_condi( double x, const struct unur_distr *condi );
struct unur_distr *
unur_distr_condi_new( const struct unur_distr *distr, const double *pos, const double *dir, int k )
{
  struct unur_distr *condi;
  double *ar;
  _unur_check_NULL( distr_name, distr, NULL );
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);
  _unur_check_NULL( distr_name, pos, NULL );
  if ( dir==NULL && (k<0 || k>=distr->dim)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"k < 0 or k >= dim");
    return NULL;
  }
  condi = unur_distr_cont_new();
  if (!condi) return NULL;
  condi->id = UNUR_DISTR_CONDI;
  condi->name = distr_name;
  condi->base = _unur_distr_cvec_clone( distr );
  if (!condi->base) { _unur_distr_free(condi); return NULL; }
  CONDI.n_params = 1;                 
  if ( unur_distr_condi_set_condition( condi, pos, dir, k ) != UNUR_SUCCESS ) {
    _unur_distr_free(condi); return NULL; 
  }
  ar = _unur_xmalloc( distr->dim * sizeof(double) );
  memset( ar, 0, distr->dim * sizeof(double) );
  if ( (unur_distr_cont_set_pdfparams_vec( condi, iXARG, ar, distr->dim ) != UNUR_SUCCESS) ||
       (unur_distr_cont_set_pdfparams_vec( condi, iGRADF, ar, distr->dim ) != UNUR_SUCCESS) ) {
    _unur_distr_free(condi); free(ar); return NULL; 
  }
  free(ar);
  if (DISTR.pdf) {
    CONDI.pdf = _unur_pdf_condi;      
    if (DISTR.dpdf)
      CONDI.dpdf = _unur_dpdf_condi;  
  }
  if (DISTR.logpdf) {
    CONDI.logpdf = _unur_logpdf_condi;      
    if (DISTR.dlogpdf)
      CONDI.dlogpdf = _unur_dlogpdf_condi;  
  }
  return condi;
} 
int
unur_distr_condi_set_condition( struct unur_distr *condi, const double *pos, const double *dir, int k )
{
  int dim;
  double *domain;
  _unur_check_NULL( distr_name, condi, UNUR_ERR_NULL );
  _unur_check_distr_object( condi, CONT, UNUR_ERR_DISTR_INVALID );
  if (condi->id != UNUR_DISTR_CONDI) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); 
    return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(condi,CK_DISTR_CONT,UNUR_ERR_COOKIE);
  dim = condi->base->dim;
  _unur_check_NULL( condi->name, pos, UNUR_ERR_NULL );
  if ( dir== NULL && (k<0 || k>=dim)) {
    _unur_error(condi->name,UNUR_ERR_DISTR_INVALID,"k < 0 or k >= dim");
    return UNUR_ERR_DISTR_INVALID;
  }
  K = (double) k;
  if ( unur_distr_cont_set_pdfparams_vec( condi, iPOSITION, pos, dim ) != UNUR_SUCCESS ) {
    return UNUR_ERR_DISTR_INVALID;
  }
  if ( unur_distr_cont_set_pdfparams_vec( condi, iDIRECTION, dir, dim ) != UNUR_SUCCESS ) {
    return UNUR_ERR_DISTR_INVALID;
  }
  if ( (domain = condi->base->data.cvec.domainrect) != NULL ) {
    if (dir == NULL) {
      CONDI.trunc[0] = CONDI.domain[0] = domain[2*k];
      CONDI.trunc[1] = CONDI.domain[1] = domain[2*k+1];
    }
    else {
      CONDI.trunc[0] = CONDI.domain[0] = -UNUR_INFINITY;
      CONDI.trunc[1] = CONDI.domain[1] = UNUR_INFINITY;
    }
  }
  condi->set &= ~UNUR_DISTR_SET_MODE; 
  return UNUR_SUCCESS;
} 
int
unur_distr_condi_get_condition( struct unur_distr *condi, const double **pos, const double **dir, int *k )
{
  _unur_check_NULL( distr_name, condi, UNUR_ERR_NULL );
  _unur_check_distr_object( condi, CONT, UNUR_ERR_DISTR_INVALID );
  if (condi->id != UNUR_DISTR_CONDI) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); 
    return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(condi,CK_DISTR_CONT,UNUR_ERR_COOKIE);
  *k = (int) K;
  *pos = POSITION;
  *dir = DIRECTION;
  return UNUR_SUCCESS;
} 
const struct unur_distr *
unur_distr_condi_get_distribution( const struct unur_distr *condi )
{
  _unur_check_NULL( distr_name, condi, NULL );
  _unur_check_distr_object( condi, CONT, NULL );
  if (condi->id != UNUR_DISTR_CONDI) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL; 
  }
  return condi->base;
} 
double
_unur_pdf_condi( double x, const struct unur_distr *condi )
{ 
  int dim = condi->base->dim;   
  int k = (int) K;              
  int i;
  if (DIRECTION==NULL) {  
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
  }
  else {   
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
  }
  return _unur_cvec_PDF(XARG,condi->base);
} 
double
_unur_logpdf_condi( double x, const struct unur_distr *condi )
{ 
  int dim = condi->base->dim;   
  int k = (int) K;              
  int i;
  if (DIRECTION==NULL) {  
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
  }
  else {   
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
  }
  return _unur_cvec_logPDF(XARG,condi->base);
} 
double
_unur_dpdf_condi( double x, const struct unur_distr *condi )
{
  int dim = condi->base->dim;    
  int k = (int) K;               
  int i;
  double df;
  if (DIRECTION==NULL) {  
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
    if (condi->base->data.cvec.pdpdf) {
      df = _unur_cvec_pdPDF(XARG, k, condi->base);
    }
    else {
      _unur_cvec_dPDF(GRADF, XARG, condi->base);
      df = GRADF[k];
    }
  }
  else {   
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
    _unur_cvec_dPDF(GRADF, XARG, condi->base);
    for (df=0.,i=0; i<dim; i++)
      df += GRADF[i]*DIRECTION[i];
  }  
  return df;
} 
double
_unur_dlogpdf_condi( double x, const struct unur_distr *condi )
{
  int dim = condi->base->dim;    
  int k = (int) K;               
  int i;
  double df;
  if (DIRECTION==NULL) {  
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
    if (condi->base->data.cvec.pdlogpdf) {
      df = _unur_cvec_pdlogPDF(XARG, k, condi->base);
    }
    else {
      _unur_cvec_dlogPDF(GRADF, XARG, condi->base);
      df = GRADF[k];
    }
  }
  else {   
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
    _unur_cvec_dlogPDF(GRADF, XARG, condi->base);
    for (df=0.,i=0; i<dim; i++)
      df += GRADF[i]*DIRECTION[i];
  }  
  return df;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_condi_debug( const struct unur_distr *condi, const char *genid )
{
  FILE *LOG;
  CHECK_NULL(condi,RETURN_VOID);
  COOKIE_CHECK(condi,CK_DISTR_CONT,RETURN_VOID);
  CHECK_NULL(condi->base,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = full conditional distribution of continuous multivariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,condi->name);
  fprintf(LOG,"%s:\n",genid);
  fprintf(LOG,"%s:\tcondition (at time of creation):\n",genid);
  if (DIRECTION==NULL) {
    fprintf(LOG,"%s:\tvariable = %d\n",genid,(int)(K));
    _unur_matrix_print_vector( condi->base->dim, POSITION, "\tpoint =", LOG, genid, "\t   ");
  }
  else {
    _unur_matrix_print_vector( condi->base->dim, DIRECTION, "\tdirection =", LOG, genid, "\t   ");
    _unur_matrix_print_vector( condi->base->dim, POSITION, "\tpoint =", LOG, genid, "\t   ");
  }
  fprintf(LOG,"%s:\n",genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)\n",genid,CONDI.domain[0],CONDI.domain[1]);
  fprintf(LOG,"%s:\n",genid);
  fprintf(LOG,"%s: Underlying distribution:\n",genid);
  _unur_distr_cvec_debug(condi->base, genid);
} 
#endif    
