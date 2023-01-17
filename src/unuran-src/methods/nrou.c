/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <utils/fmax_source.h>
#include <utils/unur_fp_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "nrou.h"
#include "nrou_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define NROU_RECT_SCALING (1.e-4)
#define BD_MAX    (DBL_MAX/1000.)
#define NROU_VARFLAG_VERIFY   0x002u   
#define NROU_DEBUG_REINIT    0x00000010u   
#define NROU_SET_U       0x001u     
#define NROU_SET_V       0x002u     
#define NROU_SET_CENTER  0x004u     
#define NROU_SET_R       0x008u     
#define GENTYPE "NROU"         
static struct unur_gen *_unur_nrou_init( struct unur_par *par );
static int _unur_nrou_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_nrou_create( struct unur_par *par );
static int _unur_nrou_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_nrou_clone( const struct unur_gen *gen );
static void _unur_nrou_free( struct unur_gen *gen);
static double _unur_nrou_sample( struct unur_gen *gen );
static double _unur_nrou_sample_check( struct unur_gen *gen );
static double _unur_aux_bound_umax(double x, void *p);
static double _unur_aux_bound_umin(double x, void *p);
static int _unur_nrou_rectangle( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_nrou_debug_init( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_nrou_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_nrou_par*)par->datap) 
#define GEN       ((struct unur_nrou_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    
#define _unur_nrou_getSAMPLE(gen) \
   ( ((gen)->variant & NROU_VARFLAG_VERIFY) \
     ? _unur_nrou_sample_check : _unur_nrou_sample )
struct unur_par *
unur_nrou_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_nrou_par) );
  COOKIE_SET(par,CK_NROU_PAR);
  par->distr    = distr;          
  PAR->umin      = 0.;         
  PAR->umax      = 0.;         
  PAR->vmax      = 0.;         
  PAR->center    = 0.;         
  PAR->r         = 1.;         
  par->method   = UNUR_METH_NROU;     
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_nrou_init;
  return par;
} 
int
unur_nrou_set_u( struct unur_par *par, double umin, double umax )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );
  if (!_unur_FP_greater(umax,umin)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
    return UNUR_ERR_PAR_SET;
  }
  PAR->umin = umin;
  PAR->umax = umax;
  par->set |= NROU_SET_U;
  return UNUR_SUCCESS;
} 
int
unur_nrou_set_v( struct unur_par *par, double vmax )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->vmax = vmax;
  par->set |= NROU_SET_V;
  return UNUR_SUCCESS;
} 
int
unur_nrou_set_center( struct unur_par *par, double center )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );
  PAR->center = center;
  par->set |= NROU_SET_CENTER;
  return UNUR_SUCCESS;
} 
int
unur_nrou_set_r( struct unur_par *par, double r )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r<=0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->r = r;
  par->set |= NROU_SET_R;
  return UNUR_SUCCESS;
} 
int
unur_nrou_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );
  par->variant = (verify) ? (par->variant | NROU_VARFLAG_VERIFY) : (par->variant & (~NROU_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_nrou_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, NROU, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= NROU_VARFLAG_VERIFY;
  else
    gen->variant &= ~NROU_VARFLAG_VERIFY;
  SAMPLE = _unur_nrou_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_nrou_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_NROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NROU_PAR,NULL);
  gen = _unur_nrou_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_nrou_check_par(gen) != UNUR_SUCCESS) {
    _unur_nrou_free(gen); return NULL;
  }
  if (_unur_nrou_rectangle(gen)!=UNUR_SUCCESS) {
    _unur_error(gen->genid , UNUR_ERR_GEN_CONDITION, "Cannot compute bounding rectangle");  
    _unur_nrou_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_nrou_debug_init(gen);
#endif
  return gen;
} 
int
_unur_nrou_reinit( struct unur_gen *gen )
{
  int rcode;
  gen->set &= ~(NROU_SET_V | NROU_SET_U);
  if ( (rcode = _unur_nrou_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  SAMPLE = _unur_nrou_getSAMPLE(gen);
  return _unur_nrou_rectangle(gen);
} 
struct unur_gen *
_unur_nrou_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NROU_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_nrou_gen) );
  COOKIE_SET(gen,CK_NROU_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_nrou_getSAMPLE(gen);
  gen->destroy = _unur_nrou_free;
  gen->clone = _unur_nrou_clone;
  gen->reinit = _unur_nrou_reinit;
  GEN->umin  = PAR->umin;             
  GEN->umax  = PAR->umax;             
  GEN->vmax  = PAR->vmax;             
  GEN->center = PAR->center;          
  GEN->r = PAR->r;                    
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_nrou_info;
#endif
  return gen;
} 
int
_unur_nrou_check_par( struct unur_gen *gen )
{
  if (!(gen->set & NROU_SET_CENTER))
    GEN->center = unur_distr_cont_get_center(gen->distr) ;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_nrou_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_nrou_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NROU_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_nrou_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_NROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NROU_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_nrou_sample( struct unur_gen *gen )
{ 
  double U,V,X;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_NROU_GEN,UNUR_INFINITY);
  while (1) {
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    U = GEN->umin + _unur_call_urng(gen->urng) * (GEN->umax - GEN->umin);
    if (_unur_isone(GEN->r))
      X = U/V + GEN->center;
    else
      X = U/pow(V,GEN->r) + GEN->center;
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    if (_unur_isone(GEN->r)) {
      if (V*V <= PDF(X)) 
        return X;
    }
    else {
      if (V <= pow(PDF(X), 1./(1.+GEN->r)) )
        return X;
    }
  }
} 
double
_unur_nrou_sample_check( struct unur_gen *gen )
{ 
  double U,V,X,fx,sfx,xfx;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_NROU_GEN,UNUR_INFINITY);
  while (1) {
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    U = GEN->umin + _unur_call_urng(gen->urng) * (GEN->umax - GEN->umin);
    if (_unur_isone(GEN->r))
      X = U/V + GEN->center;
    else
      X = U/pow(V,GEN->r) + GEN->center;
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    fx = PDF(X);
    if (_unur_isone(GEN->r)) {
      sfx = sqrt(fx);
      xfx = (X-GEN->center) * sfx;
    }
    else {
      sfx = pow(fx, 1./(1.+GEN->r));
      xfx = (X-GEN->center) * pow(fx, GEN->r/(1.+GEN->r));
    }
    if ( ( sfx > (1.+DBL_EPSILON) * GEN->vmax )   
	 || (xfx < (1.+UNUR_EPSILON) * GEN->umin) 
	 || (xfx > (1.+UNUR_EPSILON) * GEN->umax) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    if (_unur_isone(GEN->r)) {
      if (V*V <= PDF(X))
        return X;
    }
    else {
      if (V <= pow(PDF(X), 1./(1.+GEN->r)) )
        return X;
    }
  }
} 
double
_unur_aux_bound_umax(double x, void *p) 
{
  struct unur_gen *gen;
  gen = p; 
  if (_unur_isone(GEN->r)) 
    return (x-GEN->center) * sqrt( _unur_cont_PDF((x),(gen->distr)) );
  else
    return (x-GEN->center) * pow( _unur_cont_PDF((x),(gen->distr)),
				 GEN->r / (1.+ GEN->r) );
}
double
_unur_aux_bound_umin(double x, void *p)
{
  return (- _unur_aux_bound_umax(x,p)) ;
}
int
_unur_nrou_rectangle( struct unur_gen *gen )
{ 
  struct unur_funct_generic faux; 
  double mode; 
  double x, cx, sx, bx;  
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_NROU_GEN, UNUR_ERR_COOKIE );
  if ((gen->set & NROU_SET_U) && (gen->set & NROU_SET_V)) {
    return UNUR_SUCCESS;
  }
  cx=GEN->center;  
  if (!(gen->set & NROU_SET_V)) {
    mode = unur_distr_cont_get_mode(gen->distr);
    if (!_unur_isfinite(mode))
      return UNUR_ERR_GENERIC;
    GEN->vmax = pow(PDF(mode), 1./(1.+GEN->r));
    GEN->vmax = GEN->vmax * ( 1. + NROU_RECT_SCALING);
    if (! _unur_isfinite(GEN->vmax)) {
      _unur_error(gen->genid , UNUR_ERR_GENERIC, "vmax not finite");  
      return UNUR_ERR_GENERIC;
    }
  }
  if (!(gen->set & NROU_SET_U)) {
    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umin;
    faux.params = gen;
    sx = _unur_isfinite(DISTR.BD_LEFT) ? (cx+DISTR.BD_LEFT)/2. : (cx-1.); 
    bx = _unur_isfinite(DISTR.BD_LEFT) ? DISTR.BD_LEFT : (-BD_MAX);
    x = (_unur_FP_same(DISTR.BD_LEFT,cx)) 
      ? cx : _unur_util_find_max(faux, bx, cx, sx);
    while (!_unur_isfinite(x) && (fabs(bx) >= UNUR_EPSILON) ) { 
       bx = bx/10.; sx = bx/2.;  
       x = _unur_util_find_max(faux, bx, cx, sx);
    }
    GEN->umin = -faux.f(x,faux.params);
    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umax;
    faux.params = gen;
    sx = _unur_isfinite(DISTR.BD_RIGHT) ? (cx+DISTR.BD_RIGHT)/2. : (cx+1.); 
    bx = _unur_isfinite(DISTR.BD_RIGHT) ? DISTR.BD_RIGHT : BD_MAX;
    x = (_unur_FP_same(DISTR.BD_RIGHT,cx)) 
      ? cx: _unur_util_find_max(faux, cx, bx, sx);
    while (!_unur_isfinite(x) && (fabs(bx) >= UNUR_EPSILON) ) { 
       bx = bx/10.; sx = bx/2.; 
       x = _unur_util_find_max(faux, cx, bx, sx);
    }
    GEN->umax = faux.f(x,faux.params);
    GEN->umin = GEN->umin - (GEN->umax-GEN->umin)*NROU_RECT_SCALING/2.;
    GEN->umax = GEN->umax + (GEN->umax-GEN->umin)*NROU_RECT_SCALING/2.;
    if (! (_unur_isfinite(GEN->umin) && _unur_isfinite(GEN->umax))) {
       _unur_error(gen->genid , UNUR_ERR_GENERIC, "umin or umax not finite");  
       return UNUR_ERR_GENERIC;
    }
  }
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_nrou_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = nrou (naive ratio-of-uniforms)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_nrou_sample",gen->genid);
  if (gen->variant & NROU_VARFLAG_VERIFY) fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);
  fprintf(LOG,"%s: r-parameter = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,NROU_SET_R);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: center = %g\n",gen->genid,GEN->center);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: Rectangle:\n",gen->genid);
  fprintf(LOG,"%s:    left  upper point = (%g,%g)\n",gen->genid,GEN->umin,GEN->vmax);
  fprintf(LOG,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN->umax,GEN->vmax);
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_nrou_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  double harea;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX ) 
      _unur_string_append(info,"\n[ Hint: %s\n\t%s ]\n",
			  "You may provide the \"mode\" or at least",
			  "the \"center\" (a point near the mode)."); 
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: NROU (Naive Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"   r = %g\n\n", GEN->r);
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   bounding rectangle = (%g,%g) x (%g,%g)\n",
		      GEN->umin,GEN->umax, 0.,GEN->vmax);
  harea = (GEN->umax - GEN->umin) * GEN->vmax;
  _unur_string_append(info,"   area(hat) = %g\n", harea);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", 2. * harea / DISTR.area);
  else
    _unur_string_append(info,"= %.2f [approx.]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   r = %g  %s\n", GEN->r,
			(gen->set & NROU_SET_R) ? "" : "[default]");
    _unur_string_append(info,"   center = %g  %s\n",GEN->center,
			(gen->set & NROU_SET_CENTER) ? "" : "[default]");
    _unur_string_append(info,"   v = %g  %s\n", GEN->vmax,
			(gen->set & NROU_SET_V) ? "" : "[numeric.]");
    _unur_string_append(info,"   u = (%g, %g)  %s\n", GEN->umin,GEN->umax,
			(gen->set & NROU_SET_U) ? "" : "[numeric.]");
    if (gen->variant & NROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & NROU_SET_V) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"v\" to avoid numerical estimate." );
    if ( !(gen->set & NROU_SET_U) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"u\" to avoid slow (and inexact) numerical estimates." );
    _unur_string_append(info,"\n");
  }
} 
#endif   
