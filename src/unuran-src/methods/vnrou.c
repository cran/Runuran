/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h> 
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "vnrou.h"
#include "vnrou_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define VNROU_VARFLAG_VERIFY   0x002u   
#define VNROU_DEBUG_REINIT   0x00000010u   
#define VNROU_SET_U       0x001u     
#define VNROU_SET_V       0x002u     
#define VNROU_SET_R       0x008u     
#define GENTYPE "VNROU"         
static struct unur_gen *_unur_vnrou_init( struct unur_par *par );
static int _unur_vnrou_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_vnrou_create( struct unur_par *par );
static struct unur_gen *_unur_vnrou_clone( const struct unur_gen *gen );
static void _unur_vnrou_free( struct unur_gen *gen);
static int _unur_vnrou_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_vnrou_sample_check( struct unur_gen *gen, double *vec );
static int _unur_vnrou_rectangle( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_vnrou_debug_init( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_vnrou_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cvec      
#define PAR       ((struct unur_vnrou_par*)par->datap) 
#define GEN       ((struct unur_vnrou_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec  
#define SAMPLE    gen->sample.cvec            
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    
#define _unur_vnrou_getSAMPLE(gen) \
   ( ((gen)->variant & VNROU_VARFLAG_VERIFY) \
     ? _unur_vnrou_sample_check : _unur_vnrou_sample_cvec )
struct unur_par *
unur_vnrou_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_vnrou_par) );
  COOKIE_SET(par,CK_VNROU_PAR);
  par->distr    = distr;      
  PAR->r		= 1.; 	      
  PAR->vmax      = 0.;         
  PAR->umin 	= NULL;       
  PAR->umax 	= NULL;       
  par->method   = UNUR_METH_VNROU;    
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_vnrou_init;
  return par;
} 
int
unur_vnrou_set_u( struct unur_par *par, double *umin, double *umax )
{
  int d; 
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );
  for (d=0; d<par->distr->dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
  }
  PAR->umin = umin;
  PAR->umax = umax;
  par->set |= VNROU_SET_U;
  return UNUR_SUCCESS;
} 
int
unur_vnrou_chg_u( struct unur_gen *gen, double *umin, double *umax )
{
  int d; 
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );
  for (d=0; d<GEN->dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
  }
  memcpy(GEN->umin, umin, GEN->dim * sizeof(double));
  memcpy(GEN->umax, umax, GEN->dim * sizeof(double));
  gen->set |= VNROU_SET_U;
  return UNUR_SUCCESS;
} 
int
unur_vnrou_set_v( struct unur_par *par, double vmax )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->vmax = vmax;
  par->set |= VNROU_SET_V;
  return UNUR_SUCCESS;
} 
int
unur_vnrou_chg_v( struct unur_gen *gen, double vmax )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  GEN->vmax = vmax;
  gen->set |= VNROU_SET_V;
  return UNUR_SUCCESS;
} 
int
unur_vnrou_set_r( struct unur_par *par, double r )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r<=0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->r = r;
  par->set |= VNROU_SET_R;
  return UNUR_SUCCESS;
} 
int
unur_vnrou_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );
  par->variant = (verify) ? (par->variant | VNROU_VARFLAG_VERIFY) : (par->variant & (~VNROU_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_vnrou_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cvec_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= VNROU_VARFLAG_VERIFY;
  else
    gen->variant &= ~VNROU_VARFLAG_VERIFY;
  SAMPLE = _unur_vnrou_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
double 
unur_vnrou_get_volumehat( const struct unur_gen *gen )
{
  double vol;
  int d;
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, VNROU, INFINITY );
  vol = GEN->vmax;  
  for (d=0; d<GEN->dim; d++) {
    vol *= (GEN->umax[d]-GEN->umin[d]);
  }
  vol *= (GEN->r*GEN->dim+1);
  return vol;
} 
struct unur_gen *
_unur_vnrou_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_VNROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VNROU_PAR,NULL);
  gen = _unur_vnrou_create(par);
  _unur_par_free(par); 
  if (!gen) return NULL; 
  if (_unur_vnrou_rectangle(gen)!=UNUR_SUCCESS) {
    _unur_vnrou_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_vnrou_debug_init(gen);
#endif
  return gen;
} 
int
_unur_vnrou_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_vnrou_rectangle(gen))!=UNUR_SUCCESS) {
    return rcode;
  }
  SAMPLE = _unur_vnrou_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & VNROU_DEBUG_REINIT) _unur_vnrou_debug_init(gen);
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_vnrou_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VNROU_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_vnrou_gen) );
  COOKIE_SET(gen,CK_VNROU_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_vnrou_getSAMPLE(gen);
  gen->destroy = _unur_vnrou_free;
  gen->clone = _unur_vnrou_clone;
  gen->reinit = _unur_vnrou_reinit;
  GEN->dim   = gen->distr->dim;       
  GEN->r     = PAR->r;                  
  GEN->vmax  = PAR->vmax;             
  GEN->umin = _unur_xmalloc( GEN->dim * sizeof(double)); 
  GEN->umax = _unur_xmalloc( GEN->dim * sizeof(double)); 
  if (PAR->umin != NULL) memcpy(GEN->umin, PAR->umin, GEN->dim * sizeof(double));
  if (PAR->umax != NULL) memcpy(GEN->umax, PAR->umax, GEN->dim * sizeof(double));
  GEN->center = unur_distr_cvec_get_center(gen->distr);
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_vnrou_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_vnrou_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_vnrou_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VNROU_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->umin = _unur_xmalloc( GEN->dim * sizeof(double));
  CLONE->umax = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy(CLONE->umin, GEN->umin, GEN->dim * sizeof(double));
  memcpy(CLONE->umax, GEN->umax, GEN->dim * sizeof(double));
  CLONE->center = unur_distr_cvec_get_center(clone->distr);
  return clone;
#undef CLONE
} 
int
_unur_vnrou_sample_cvec( struct unur_gen *gen, double *vec )
{ 
  double U, V;
  int d, dim; 
  CHECK_NULL(gen,UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_VNROU_GEN,UNUR_ERR_COOKIE); 
  dim = GEN->dim;
  while (1) {
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    for (d=0; d<dim; d++) {
      U = GEN->umin[d] + _unur_call_urng(gen->urng) * (GEN->umax[d] - GEN->umin[d]);
      vec[d] = U/pow(V,GEN->r) + GEN->center[d];
    }
    if (V <= pow(PDF(vec),1./(GEN->r * dim + 1.)))
      return UNUR_SUCCESS;
  }
} 
int
_unur_vnrou_sample_check( struct unur_gen *gen, double *vec )
{ 
  double U, V;
  int d, dim; 
  int hat_error;
  double fx,sfx,xfx;
  CHECK_NULL(gen,UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_VNROU_GEN,UNUR_ERR_COOKIE); 
  dim = GEN->dim;
  while (1) {
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    for (d=0; d<dim; d++) {
      U = GEN->umin[d] + _unur_call_urng(gen->urng) * (GEN->umax[d] - GEN->umin[d]);
      vec[d] = U/pow(V,GEN->r) + GEN->center[d];
    }
    fx = PDF(vec);
    sfx = pow( fx, 1./(GEN->r * dim+1.) );
    hat_error=0;
    if ( sfx > (1.+DBL_EPSILON) * GEN->vmax ) hat_error++;  
    sfx = pow( fx, GEN->r/(GEN->r * dim + 1.) );
    for (d=0; d<dim; d++) {
     xfx = (vec[d]-GEN->center[d]) * sfx;
     if ( (xfx < (1.+UNUR_EPSILON) * GEN->umin[d]) 
       || (xfx > (1.+UNUR_EPSILON) * GEN->umax[d]))
       hat_error++;
    }
    if (hat_error>0) _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    if (V <= pow(PDF(vec),1./( GEN->r * dim + 1.)))
      return UNUR_SUCCESS;
  }
} 
void
_unur_vnrou_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_VNROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->umin) free(GEN->umin); 
  if (GEN->umax) free(GEN->umax);
  _unur_generic_free(gen);
} 
int
_unur_vnrou_rectangle( struct unur_gen *gen )
{ 
  int d; 
  struct MROU_RECTANGLE *rr;
  int rectangle_compute;
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_VNROU_GEN, UNUR_ERR_COOKIE );
  if ((gen->set & VNROU_SET_U) && (gen->set & VNROU_SET_V)) {
    return UNUR_SUCCESS;
  }
  rr = _unur_mrou_rectangle_new();
  rr->distr  = gen->distr;
  rr->dim    = GEN->dim;
  rr->umin   = GEN->umin;
  rr->umax   = GEN->umax;
  rr->r      = GEN->r;
  rr->center = GEN->center; 
  rr->genid  = gen->genid;
  rectangle_compute = _unur_mrou_rectangle_compute(rr);
  if (!(gen->set & VNROU_SET_V)) {
     GEN->vmax = rr->vmax;
  }
  if (!(gen->set & VNROU_SET_U)) {
    for (d=0; d<GEN->dim; d++) {
      GEN->umin[d] = rr->umin[d];
      GEN->umax[d] = rr->umax[d];
    }
  }
  free(rr);
  if (rectangle_compute != UNUR_SUCCESS)
    return UNUR_ERR_INF;
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_vnrou_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  int d, dim; 
  double vol;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  dim = GEN->dim;
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = vnrou (naive ratio-of-uniforms)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cvec_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_vnrou_sample",gen->genid);
  if (gen->variant & VNROU_VARFLAG_VERIFY) fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);
  fprintf(LOG,"%s: r-parameter = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,VNROU_SET_R);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_matrix_print_vector( GEN->dim, GEN->center, "center =", LOG, gen->genid, "\t   ");
  fprintf(LOG,"%s: Rectangle:",gen->genid);
  if (!((gen->set & VNROU_SET_U) && (gen->set & VNROU_SET_V)))
    fprintf(LOG,"\t[computed]");
  else 
    fprintf(LOG,"\t[input]");
  fprintf(LOG,"\n");
  vol = GEN->vmax;
  fprintf(LOG,"%s:\tvmax = %g\n",gen->genid, GEN->vmax);
  for (d=0; d<dim; d++) {
    vol *= (GEN->umax[d]-GEN->umin[d]);
    fprintf(LOG,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid, 
	    d, d, GEN->umin[d], GEN->umax[d]);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:\tvolume = %g\t(hat = %g)\n",gen->genid, vol, vol*(GEN->r*GEN->dim+1));
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_vnrou_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  int i;
  double hvol;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_distr_cvec_info_domain(gen);
  if ( distr->set & UNUR_DISTR_SET_MODE ) {
    _unur_string_append(info,"   mode      = ");
    _unur_distr_info_vector( gen, DISTR.mode, GEN->dim);
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   center    = ");
  _unur_distr_info_vector( gen, GEN->center, GEN->dim);
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]");
    else
      _unur_string_append(info,"  [default]");
  }
  _unur_string_append(info,"\n\n");
  _unur_string_append(info,"method: VNROU (Naive Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"   r = %g\n", GEN->r);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   bounding rectangle = ");
  for (i=0; i<GEN->dim; i++)
    _unur_string_append(info,"%s(%g,%g)", i?"x":"", GEN->umin[i], GEN->umax[i]);
  _unur_string_append(info," x (0,%g)\n", GEN->vmax);
  hvol = GEN->vmax;
  for (i=0; i<GEN->dim; i++)
    hvol *= GEN->umax[i] - GEN->umin[i];
  _unur_string_append(info,"   volume(hat) = %g\n", hvol);
  _unur_string_append(info,"   rejection constant ");
  if ((distr->set & UNUR_DISTR_SET_PDFVOLUME) && _unur_isone(GEN->r))
    _unur_string_append(info,"= %g\n", (GEN->dim + 1.) * hvol / DISTR.volume);
  else
    _unur_string_append(info,"= %.2f  [approx.]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/((1.+GEN->dim)*samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   r = %g  %s\n", GEN->r,
			(gen->set & VNROU_SET_R) ? "" : "[default]");
    _unur_string_append(info,"   v = %g  %s\n", GEN->vmax,
 			(gen->set & VNROU_SET_V) ? "" : "[numeric.]");
    _unur_string_append(info,"   u = ");
    _unur_distr_info_vector( gen, GEN->umin, GEN->dim);
    _unur_string_append(info," -- ");
    _unur_distr_info_vector( gen, GEN->umax, GEN->dim);
    _unur_string_append(info,"%s\n",(gen->set & VNROU_SET_U) ? "" : "  [numeric.]"); 
    if (gen->variant & VNROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & VNROU_SET_V) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"v\" to avoid numerical estimate." );
    if ( !(gen->set & VNROU_SET_U) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"u\" to avoid slow (and inexact) numerical estimates." );
    _unur_string_append(info,"\n");
  }
} 
#endif   
