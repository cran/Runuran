/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr.h"
#include "cont.h"
#include "corder.h"
#include "distr_source.h"
static const char distr_name[] = "order statistics";
#define DISTR distr->data.cont    
#define OS    os->data.cont       
#define LOGNORMCONSTANT (os->data.cont.norm_constant)
static double _unur_pdf_corder( double x, const struct unur_distr *os );
static double _unur_dpdf_corder( double x, const struct unur_distr *os );
static double _unur_cdf_corder( double x, const struct unur_distr *os );
static int _unur_upd_area_corder( struct unur_distr *os );
struct unur_distr *
unur_distr_corder_new( const struct unur_distr *distr, int n, int k )
{
  struct unur_distr *os;
  _unur_check_NULL( distr_name,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (distr->id == UNUR_DISTR_CORDER) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"No order statistics of order statistics allowed");
    return NULL; 
  }
  if (n < 2 || k < 1 || k > n) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"n < 2 or k < 1 or k > n");
    return NULL;
  }
  os = unur_distr_cont_new();
  if (!os) return NULL;
  os->id = UNUR_DISTR_CORDER;
  os->name = distr_name;
  os->base = _unur_distr_cont_clone( distr );
  if (!os->base) { free(os); return NULL; }
  OS.n_params = 2;                 
  OS.params[0] = (double) n;
  OS.params[1] = (double) k;
  OS.area = DISTR.area;               
  OS.trunc[0] = OS.domain[0] = DISTR.domain[0];  
  OS.trunc[1] = OS.domain[1] = DISTR.domain[1];  
  if (DISTR.cdf) {
    OS.cdf = _unur_cdf_corder;        
    if (DISTR.pdf) {
      OS.pdf = _unur_pdf_corder;      
      if (DISTR.dpdf)
	OS.dpdf = _unur_dpdf_corder;  
    }
  }
  OS.upd_area  = _unur_upd_area_corder;
  os->set = distr->set & ~UNUR_DISTR_SET_MODE; 
  if (_unur_upd_area_corder(os)==UNUR_SUCCESS)
    os->set |= UNUR_DISTR_SET_PDFAREA;
  return os;
} 
const struct unur_distr *
unur_distr_corder_get_distribution( const struct unur_distr *os )
{
  _unur_check_NULL( distr_name, os, NULL );
  _unur_check_distr_object( os, CONT, NULL );
  if (os->id != UNUR_DISTR_CORDER) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL; 
  }
  return os->base;
} 
int
unur_distr_corder_set_rank( struct unur_distr *os, int n, int k )
{
  _unur_check_NULL( distr_name, os, UNUR_ERR_NULL );
  _unur_check_distr_object( os, CONT, UNUR_ERR_DISTR_INVALID );
  if (os->id != UNUR_DISTR_CORDER) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(os,CK_DISTR_CONT,UNUR_ERR_COOKIE);
  if (n < 2 || k < 1 || k > n) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"n < 2 or k < 1 or k > n");
    return UNUR_ERR_DISTR_SET;
  }
  os->set &= ~UNUR_DISTR_SET_MODE; 
  OS.params[0] = (double) n;
  OS.params[1] = (double) k;
  _unur_upd_area_corder(os);
  return UNUR_SUCCESS;
} 
int 
unur_distr_corder_get_rank( const struct unur_distr *os, int *n, int *k )
{
  _unur_check_NULL( distr_name, os, UNUR_ERR_NULL );
  _unur_check_distr_object( os, CONT, UNUR_ERR_DISTR_INVALID );
  if (os->id != UNUR_DISTR_CORDER) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(os,CK_DISTR_CONT,UNUR_ERR_COOKIE);
  *n = (int)(OS.params[0] + 0.5);
  *k = (int)(OS.params[1] + 0.5);
  return UNUR_SUCCESS;
} 
double
_unur_pdf_corder( double x, const struct unur_distr *os )
{ 
  double Fx;    
  double fx;    
  double p,q;   
  _unur_check_NULL( NULL, os, UNUR_INFINITY );
  _unur_check_distr_object( os, CONT, UNUR_INFINITY );
  CHECK_NULL( os->base, UNUR_INFINITY );
  _unur_check_distr_object( os->base, CONT, UNUR_INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, UNUR_INFINITY );
  CHECK_NULL( os->base->data.cont.pdf, UNUR_INFINITY );
  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  fx = (*(os->base->data.cont.pdf)) (x, os->base);
  p = OS.params[1];                       
  q = OS.params[0] - OS.params[1] + 1.;   
  if (fx <= 0. || Fx <= 0. || Fx >= 1.)
    return 0.;
  else
    return exp(log(fx) + (p-1.)*log(Fx) + (q-1.)*log(1.-Fx) - LOGNORMCONSTANT);
} 
double
_unur_dpdf_corder( double x, const struct unur_distr *os )
{
  double Fx;    
  double fx;    
  double dfx;   
  double p,q;   
  double dpdf;  
  double lFx, lFy;
  _unur_check_NULL( NULL, os, UNUR_INFINITY );
  _unur_check_distr_object( os, CONT, UNUR_INFINITY );
  CHECK_NULL( os->base, UNUR_INFINITY );
  _unur_check_distr_object( os->base, CONT, UNUR_INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, UNUR_INFINITY );
  CHECK_NULL( os->base->data.cont.pdf, UNUR_INFINITY );
  CHECK_NULL( os->base->data.cont.dpdf, UNUR_INFINITY );
  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  fx = (*(os->base->data.cont.pdf)) (x, os->base);
  dfx = (*(os->base->data.cont.dpdf)) (x, os->base);
  p = OS.params[1];                       
  q = OS.params[0] - OS.params[1] + 1.;   
  if (fx <= 0. || Fx <= 0. || Fx >= 1.)
    return 0.;
  lFx = log(Fx);
  lFy = log(1.-Fx);
  dpdf = ( exp(2.*log(fx) + (p-2.)*lFx + (q-2.)*lFy - LOGNORMCONSTANT)
	   * ( (p-1.)*(1.-Fx) - (q-1.)*Fx ));
  dpdf += exp((p-1.)*lFx + (q-1.)*lFy - LOGNORMCONSTANT) * dfx;
  return dpdf;
} 
double
_unur_cdf_corder( double x, const struct unur_distr *os ) 
{
  double Fx;    
  double p,q;   
  _unur_check_NULL( NULL, os, UNUR_INFINITY );
  _unur_check_distr_object( os, CONT, UNUR_INFINITY );
  CHECK_NULL( os->base, UNUR_INFINITY );
  _unur_check_distr_object( os->base, CONT, UNUR_INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, UNUR_INFINITY );
  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  p = OS.params[1];                       
  q = OS.params[0] - OS.params[1] + 1.;   
  return _unur_SF_incomplete_beta(Fx,p,q);
} 
int
_unur_upd_area_corder( UNUR_DISTR *os )
{
  LOGNORMCONSTANT = ( _unur_SF_ln_gamma(OS.params[1]) 
		      + _unur_SF_ln_gamma(OS.params[0] - OS.params[1] + 1.) 
		      - _unur_SF_ln_gamma(OS.params[0] + 1.) );
  if (!(os->set & UNUR_DISTR_SET_STDDOMAIN)) {
    if (OS.cdf == NULL)
      return UNUR_ERR_DISTR_REQUIRED;
    OS.area  = (OS.domain[1] < UNUR_INFINITY)  ? _unur_cdf_corder(OS.domain[1],os) : 1.;
    OS.area -= (OS.domain[0] > -UNUR_INFINITY) ? _unur_cdf_corder(OS.domain[0],os) : 0.;
  }
  return (OS.area > 0.) ? UNUR_SUCCESS : UNUR_ERR_DISTR_DATA;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_corder_debug( const struct unur_distr *os, const char *genid )
{
  FILE *LOG;
  CHECK_NULL(os,RETURN_VOID);
  COOKIE_CHECK(os,CK_DISTR_CONT,RETURN_VOID);
  CHECK_NULL(os->base,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = order statistics of continuous univariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,os->name);
  fprintf(LOG,"%s:\tsample size\tn = %d\n",genid,(int)(OS.params[0]+0.5));
  fprintf(LOG,"%s:\trank\t\tk = %d\n",genid,(int)(OS.params[1]+0.5));
  fprintf(LOG,"%s:\n",genid);
  if (os->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %g\n",genid,OS.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)",genid,OS.domain[0],OS.domain[1]);
  _unur_print_if_default(os,UNUR_DISTR_SET_DOMAIN);
  fprintf(LOG,"\n%s:\tarea below PDF = %g",genid,OS.area);
  _unur_print_if_default(os,UNUR_DISTR_SET_PDFAREA);
  fprintf(LOG,"\n%s:\n",genid);
  fprintf(LOG,"%s: Underlying distribution:\n",genid);
  _unur_distr_cont_debug(os->base, genid);
} 
#endif    
