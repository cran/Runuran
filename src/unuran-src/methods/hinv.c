/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "hinv.h"
#include "hinv_struct.h"
#define HINV_MAX_ITER      (300)
#define HINV_MAX_U_LENGTH  (0.05)
#define HINV_TAILCUTOFF    (1.e-10) 
#define HINV_XDEVIATION    (0.05)
#define HINV_DEBUG_REINIT    0x00000002u   
#define HINV_DEBUG_TABLE     0x00000010u   
#define HINV_DEBUG_CHG       0x00001000u   
#define HINV_DEBUG_SAMPLE    0x01000000u   
#define HINV_SET_ORDER          0x001u  
#define HINV_SET_U_RESOLUTION   0x002u  
#define HINV_SET_STP            0x004u  
#define HINV_SET_BOUNDARY       0x008u  
#define HINV_SET_GUIDEFACTOR    0x010u  
#define HINV_SET_MAX_IVS        0x020u  
#define GENTYPE "HINV"         
static struct unur_gen *_unur_hinv_init( struct unur_par *par );
static int _unur_hinv_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_hinv_create( struct unur_par *par );
static int _unur_hinv_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_hinv_clone( const struct unur_gen *gen );
static void _unur_hinv_free( struct unur_gen *gen );
static double _unur_hinv_sample( struct unur_gen *gen );
static double _unur_hinv_eval_approxinvcdf( const struct unur_gen *gen, double u );
static int _unur_hinv_find_boundary( struct unur_gen *gen );
static int _unur_hinv_create_table( struct unur_gen *gen );
static struct unur_hinv_interval *_unur_hinv_interval_new( struct unur_gen *gen, double p, double u );
static struct unur_hinv_interval *_unur_hinv_interval_adapt( struct unur_gen *gen, 
							     struct unur_hinv_interval *iv, 
							     int *error_count_shortinterval );
static int _unur_hinv_interval_is_monotone( struct unur_gen *gen, struct unur_hinv_interval *iv );
static int _unur_hinv_interval_parameter( struct unur_gen *gen, struct unur_hinv_interval *iv );
static double _unur_hinv_eval_polynomial( double x, double *coeff, int order );
static int _unur_hinv_list_to_array( struct unur_gen *gen );
static int _unur_hinv_make_guide_table( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_hinv_debug_init( const struct unur_gen *gen, int ok);
static void _unur_hinv_debug_intervals( const struct unur_gen *gen );
static void _unur_hinv_debug_chg_truncated( const struct unur_gen *gen);
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_hinv_par*)par->datap) 
#define GEN       ((struct unur_hinv_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define SAMPLE    gen->sample.cont      
#define CDF(x)  ((_unur_cont_CDF((x),(gen->distr))-GEN->CDFmin)/(GEN->CDFmax-GEN->CDFmin))
#define PDF(x)  (_unur_cont_PDF((x),(gen->distr))/(GEN->CDFmax-GEN->CDFmin)) 
#define dPDF(x) (_unur_cont_dPDF((x),(gen->distr))/(GEN->CDFmax-GEN->CDFmin))
#define _unur_hinv_getSAMPLE(gen)  (_unur_hinv_sample)
struct unur_par *
unur_hinv_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_hinv_par) );
  COOKIE_SET(par,CK_HINV_PAR);
  par->distr   = distr;           
  PAR->order = (DISTR_IN.pdf) ? 3 : 1;  
  PAR->u_resolution = 1.0e-10;    
  PAR->guide_factor = 1.;         
  PAR->bleft = -1.e20;            
  PAR->bright = 1.e20;            
  PAR->max_ivs = 1000000;         
  PAR->stp = NULL;                
  PAR->n_stp = 0;                 
  par->method   = UNUR_METH_HINV; 
  par->variant  = 0u;             
  par->set      = 0u;                      
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_hinv_init;
  return par;
} 
int
unur_hinv_set_order( struct unur_par *par, int order)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HINV );
  if (order!=1 && order!=3 && order!=5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"order");
    return UNUR_ERR_PAR_SET;
  }
  if (order > 1 && par->distr->data.cont.pdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    return UNUR_ERR_DISTR_REQUIRED;
  }
  if (order > 3 && par->distr->data.cont.dpdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dPDF");
    return UNUR_ERR_DISTR_REQUIRED;
  }
  PAR->order = order;
  par->set |= HINV_SET_ORDER;
  return UNUR_SUCCESS;
} 
int
unur_hinv_set_u_resolution( struct unur_par *par, double u_resolution )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HINV );
  if (u_resolution < 5.*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution");
    return UNUR_ERR_PAR_SET;
  }
  if (u_resolution < UNUR_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution so small that problems may occur");
  }
  PAR->u_resolution = u_resolution;
  par->set |= HINV_SET_U_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_hinv_set_cpoints( struct unur_par *par, const double *stp, int n_stp )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HINV );
  if (n_stp < 1 || stp==NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 1");
    return UNUR_ERR_PAR_SET;
  }
  for( i=1; i<n_stp; i++ )
    if (stp[i] <= stp[i-1]) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
      return UNUR_ERR_PAR_SET;
    }
  PAR->stp = stp;
  PAR->n_stp = n_stp;
  par->set |= HINV_SET_STP;
  return UNUR_SUCCESS;
} 
int
unur_hinv_set_boundary( struct unur_par *par, double left, double right )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HINV );
  if (left >= right) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (left <= -INFINITY || right >= INFINITY) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }
  PAR->bleft = left;
  PAR->bright = right;
  par->set |= HINV_SET_BOUNDARY;
  return UNUR_SUCCESS;
} 
int
unur_hinv_set_guidefactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HINV );
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->guide_factor = factor;
  par->set |= HINV_SET_GUIDEFACTOR;
  return UNUR_SUCCESS;
} 
int
unur_hinv_set_max_intervals( struct unur_par *par, int max_ivs )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HINV );
  if (max_ivs < 1000 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1000");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ivs = max_ivs;
  par->set |= HINV_SET_MAX_IVS;
  return UNUR_SUCCESS;
} 
int
unur_hinv_get_n_intervals( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, HINV, 0 );
  return GEN->N;
} 
int 
unur_hinv_chg_truncated( struct unur_gen *gen, double left, double right )
{
  double Umin, Umax, Uminbound, Umaxbound;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object(gen, HINV, UNUR_ERR_GEN_INVALID);
  if (left < GEN->bleft) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, increase left boundary");
    left = GEN->bleft;
  }
  if (right > GEN->bright) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, decrease right boundary");
    right = GEN->bright;
  }
  if (!_unur_FP_less(left,right)) {
    _unur_error(gen->genid,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  Uminbound = _unur_max(0.,GEN->intervals[0]);
  Umaxbound = _unur_min(1.,GEN->intervals[(GEN->N-1)*(GEN->order+2)]);
  Umin = (left > -INFINITY) ? CDF(left)  : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;
  if (Umin > Umax) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  if (_unur_FP_equal(Umin,Umax)) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return UNUR_ERR_DISTR_SET;
    }
  }
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = _unur_max(Umin, Uminbound);
  GEN->Umax = _unur_min(Umax, Umaxbound);
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & HINV_DEBUG_CHG) 
    _unur_hinv_debug_chg_truncated( gen );
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hinv_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_HINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HINV_PAR,NULL);
  gen = _unur_hinv_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_hinv_check_par(gen) != UNUR_SUCCESS) {
    _unur_hinv_free(gen); return NULL;
  }
  if (_unur_hinv_create_table(gen)!=UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) {
      _unur_hinv_list_to_array( gen );
      _unur_hinv_debug_init(gen,FALSE);
    }
#endif
    _unur_hinv_free(gen); return NULL;
  }
  _unur_hinv_list_to_array( gen );
  GEN->Umin = _unur_max(0.,GEN->intervals[0]);
  GEN->Umax = _unur_min(1.,GEN->intervals[(GEN->N-1)*(GEN->order+2)]);
  _unur_hinv_make_guide_table(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_hinv_debug_init(gen,TRUE);
#endif
  GEN->stp = NULL;
  GEN->n_stp = 0;
  return gen;
} 
int
_unur_hinv_reinit( struct unur_gen *gen )
{
  int rcode;
  SAMPLE = _unur_hinv_getSAMPLE(gen);
  if ( (rcode = _unur_hinv_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  if ( (rcode = _unur_hinv_create_table(gen)) != UNUR_SUCCESS)
    return rcode;
  _unur_hinv_list_to_array( gen );
  GEN->Umin = _unur_max(0.,GEN->intervals[0]);
  GEN->Umax = _unur_min(1.,GEN->intervals[(GEN->N-1)*(GEN->order+2)]);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & HINV_DEBUG_REINIT) _unur_hinv_debug_init(gen,TRUE);
#endif
  _unur_hinv_make_guide_table(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hinv_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HINV_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_hinv_gen) );
  COOKIE_SET(gen,CK_HINV_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_hinv_getSAMPLE(gen);
  gen->destroy = _unur_hinv_free;
  gen->clone = _unur_hinv_clone;
  gen->reinit = _unur_hinv_reinit;
  GEN->order = PAR->order;            
  GEN->u_resolution = PAR->u_resolution; 
  GEN->guide_factor = PAR->guide_factor; 
  GEN->bleft_par  = PAR->bleft;          
  GEN->bright_par = PAR->bright;
  GEN->max_ivs = PAR->max_ivs;           
  GEN->stp = PAR->stp;               
  GEN->n_stp = PAR->n_stp;           
  GEN->tailcutoff_left  = -1.;       
  GEN->tailcutoff_right = 10.;
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  GEN->Umin = 0.;
  GEN->Umax = 1.;
  GEN->N = 0;
  GEN->iv = NULL;
  GEN->intervals = NULL;
  GEN->guide_size = 0; 
  GEN->guide = NULL;
  return gen;
} 
int
_unur_hinv_check_par( struct unur_gen *gen )
{
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];
  GEN->CDFmin = (DISTR.domain[0] > -INFINITY) ? _unur_cont_CDF((DISTR.domain[0]),(gen->distr)) : 0.;
  GEN->CDFmax = (DISTR.domain[1] < INFINITY)  ? _unur_cont_CDF((DISTR.domain[1]),(gen->distr)) : 1.;
  if (!_unur_FP_less(GEN->CDFmin,GEN->CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    return UNUR_ERR_GEN_DATA;
  }
  if (DISTR.domain[0] <= -INFINITY || 
      (DISTR.pdf!=NULL && _unur_cont_PDF((DISTR.domain[0]),(gen->distr))<=0.) ) {
    GEN->tailcutoff_left = _unur_min(HINV_TAILCUTOFF, 0.1*GEN->u_resolution);
    GEN->tailcutoff_left = _unur_max(GEN->tailcutoff_left,2*DBL_EPSILON);
  }
  if (DISTR.domain[1] >= INFINITY || 
      (DISTR.pdf!=NULL && _unur_cont_PDF((DISTR.domain[1]),(gen->distr))<=0.) ) {
    GEN->tailcutoff_right = _unur_min(HINV_TAILCUTOFF, 0.1*GEN->u_resolution);
    GEN->tailcutoff_right = _unur_max(GEN->tailcutoff_right,2*DBL_EPSILON);
    GEN->tailcutoff_right = 1. - GEN->tailcutoff_right;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hinv_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_hinv_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->intervals = _unur_xmalloc( GEN->N*(GEN->order+2) * sizeof(double) );
  memcpy( CLONE->intervals, GEN->intervals, GEN->N*(GEN->order+2) * sizeof(double) );
  CLONE->guide = _unur_xmalloc( GEN->guide_size * sizeof(int) );
  memcpy( CLONE->guide, GEN->guide, GEN->guide_size * sizeof(int) );
  return clone;
#undef CLONE
} 
void
_unur_hinv_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_HINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HINV_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->iv) {
    struct unur_hinv_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }
  if (GEN->intervals) free (GEN->intervals);
  if (GEN->guide)     free (GEN->guide);
  _unur_generic_free(gen);
} 
double
_unur_hinv_sample( struct unur_gen *gen )
{ 
  double U,X;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HINV_GEN,INFINITY);
  U = GEN->Umin + _unur_call_urng(gen->urng) * (GEN->Umax - GEN->Umin);
  X = _unur_hinv_eval_approxinvcdf(gen,U);
  if (X<DISTR.trunc[0]) return DISTR.trunc[0];
  if (X>DISTR.trunc[1]) return DISTR.trunc[1];
  return X;
} 
double
_unur_hinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
{ 
  int i;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HINV_GEN,INFINITY);
  i =  GEN->guide[(int) (GEN->guide_size*u)];
  while (u > GEN->intervals[i+GEN->order+2])
    i += GEN->order+2;
  u = (u-GEN->intervals[i])/(GEN->intervals[i+GEN->order+2] - GEN->intervals[i]);
  return _unur_hinv_eval_polynomial( u, GEN->intervals+i+1, GEN->order );
} 
double
unur_hinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
{ 
  double x;
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_HINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_HINV_GEN,INFINITY);
  if ( u<0. || u>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }
  if (u<=0.) return DISTR.domain[0];
  if (u>=1.) return DISTR.domain[1];
  x = _unur_hinv_eval_approxinvcdf(gen,u);
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];
  return x;
} 
int
unur_hinv_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE )
{ 
  double U, ualt, X;
  double max=0., average=0., uerror, errorat=0.;
  int j, outside_interval=0;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,UNUR_ERR_COOKIE);
  for(j=0;j<samplesize;j++) {  
    U = GEN->Umin + _unur_call_urng(gen->urng) * (GEN->Umax - GEN->Umin);
    ualt=U;
    X = _unur_hinv_eval_approxinvcdf(gen,U);
    if (X<DISTR.trunc[0]) { X = DISTR.trunc[0]; outside_interval++; }
    if (X>DISTR.trunc[1]) { X = DISTR.trunc[1]; outside_interval++; }
    uerror = fabs(ualt-CDF(X));
    average += uerror;
    if(uerror>max) {
      max = uerror;
      errorat = X;
    }
  }
  *max_error = max;
  *MAE = average/samplesize;
  return UNUR_SUCCESS;
} 
int
_unur_hinv_find_boundary( struct unur_gen *gen )
{
  double x,u;
  int i;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,UNUR_ERR_COOKIE);
  GEN->N = 0;
  if (GEN->bleft  < DISTR.domain[0]) GEN->bleft  = DISTR.domain[0];
  if (GEN->bright > DISTR.domain[1]) GEN->bright = DISTR.domain[1];
  for (x = GEN->bleft, i=0; i<HINV_MAX_ITER; i++) {
    GEN->bleft = x;
    u = CDF(GEN->bleft);
    if (u <= GEN->tailcutoff_left || GEN->tailcutoff_left < 0.) 
      break;
    if (DISTR.domain[0] <= -INFINITY) {
      x = (GEN->bleft > -1.) ? -1. : 10.*GEN->bleft;
      if (! _unur_isfinite(x) )  
	i = HINV_MAX_ITER;
    }
    else {
      x = _unur_arcmean(GEN->bleft, DISTR.domain[0]);
      if (_unur_FP_equal(x,DISTR.domain[0])) 
	i = HINV_MAX_ITER;
    }
  }
  if (i >= HINV_MAX_ITER)
    _unur_warning(gen->genid,UNUR_ERR_DISTR_PROP,"cannot find l.h.s. of domain");
  GEN->iv = _unur_hinv_interval_new(gen,GEN->bleft,u);
  if (GEN->iv == NULL) return UNUR_ERR_GEN_DATA;
  for (x = GEN->bright, i=0; i<HINV_MAX_ITER; i++) {
    GEN->bright = x;
    u = CDF(GEN->bright);
    if (u >= GEN->tailcutoff_right || GEN->tailcutoff_right > 1.1) 
      break;
    if (DISTR.domain[1] >= INFINITY) {
      x = (GEN->bright < 1.) ? 1. : 10.*GEN->bright;
      if (! _unur_isfinite(x) )  
	i = HINV_MAX_ITER;
    }
    else {
      x = _unur_arcmean(GEN->bright, DISTR.domain[1]);
      if (_unur_FP_equal(x,DISTR.domain[1])) 
	i = HINV_MAX_ITER;
    }
  }
  if (i >= HINV_MAX_ITER)
    _unur_warning(gen->genid,UNUR_ERR_DISTR_PROP,"cannot find r.h.s. of domain");
  GEN->iv->next = _unur_hinv_interval_new(gen,GEN->bright,u);
  if (GEN->iv->next == NULL) return UNUR_ERR_GEN_DATA;
  return UNUR_SUCCESS;
} 
int
_unur_hinv_create_table( struct unur_gen *gen )
{
  struct unur_hinv_interval *iv, *iv_new;
  int i, error_count_shortinterval=0;
  double Fx;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,UNUR_ERR_COOKIE);
  if (_unur_hinv_find_boundary(gen) != UNUR_SUCCESS)
    return UNUR_ERR_GEN_DATA;
  if (GEN->stp) {
    iv = GEN->iv;
    for (i=0; i<GEN->n_stp; i++) {
      if (!_unur_FP_greater(GEN->stp[i],GEN->bleft)) continue; 
      if (!_unur_FP_less(GEN->stp[i],GEN->bright))   break;    
      Fx = CDF(GEN->stp[i]);
      iv_new = _unur_hinv_interval_new(gen,GEN->stp[i],Fx);
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;
      iv_new->next = iv->next;
      iv->next = iv_new;
      iv = iv_new;
      if (Fx > GEN->tailcutoff_right)
	break;
    }
  }
  else 
    if( (gen->distr->set & UNUR_DISTR_SET_MODE) &&
        _unur_FP_greater(DISTR.mode, GEN->bleft) &&
        _unur_FP_less(DISTR.mode, GEN->bright) ) {
      iv = GEN->iv;
      iv_new = _unur_hinv_interval_new(gen,DISTR.mode,CDF(DISTR.mode));
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;
      iv_new->next = iv->next;
      iv->next = iv_new;
    }
  for (iv=GEN->iv; iv->next!=NULL; ) {
    COOKIE_CHECK(iv,CK_HINV_IV,UNUR_ERR_COOKIE);
    if (GEN->N >= GEN->max_ivs) {
      _unur_error(GENTYPE,UNUR_ERR_GEN_CONDITION,"too many intervals");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv = _unur_hinv_interval_adapt(gen,iv, &error_count_shortinterval);
    if (iv == NULL) return UNUR_ERR_GEN_DATA;
  }
  iv->spline[0] = iv->p;
  return UNUR_SUCCESS;
}  
struct unur_hinv_interval *
_unur_hinv_interval_new( struct unur_gen *gen, double p, double u )
{
  struct unur_hinv_interval *iv;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,NULL);
  if (u<0.) {
    if (u < -UNUR_SQRT_DBL_EPSILON) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF(x) < 0.");
      return NULL;
    }
    else { 
      u = 0.;
    }
  }
  if (u>1.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF(x) > 1.");
    return NULL;
  }
  iv = _unur_xmalloc( sizeof(struct unur_hinv_interval) );
  COOKIE_SET(iv,CK_HINV_IV);
  switch (GEN->order) {
  case 5:
    iv->df = dPDF(p);
  case 3:
    iv->f = PDF(p);
  case 1:
    iv->p = p;
    iv->u = u;
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    free(iv);
    return NULL;
  }
  iv->next = NULL;  
  ++(GEN->N);   
  return iv;
} 
struct unur_hinv_interval *
_unur_hinv_interval_adapt( struct unur_gen *gen, struct unur_hinv_interval *iv,
                           int *error_count_shortinterval )
{
  double p_new;   
  struct unur_hinv_interval *iv_new, *iv_tmp;
  double x, Fx;
  CHECK_NULL(gen,NULL);       COOKIE_CHECK(gen,CK_HINV_GEN,NULL);
  CHECK_NULL(iv,NULL);        COOKIE_CHECK(iv,CK_HINV_IV,NULL);
  CHECK_NULL(iv->next,NULL);  COOKIE_CHECK(iv->next,CK_HINV_IV,NULL);
  iv_tmp = iv->next->next;
  if(iv_tmp && iv->next->u > GEN->tailcutoff_right) {
    free (iv_tmp);
    iv->next->next = NULL;
    GEN->N--;
    GEN->bright = iv->next->p;
    return iv;
  }
  if (iv==GEN->iv && iv->next->next && iv->next->u < GEN->tailcutoff_left) {
    iv_tmp = GEN->iv;
    GEN->iv = iv->next;
    free (iv_tmp);
    GEN->N--;
    GEN->bleft = GEN->iv->p;
    return GEN->iv;
  }
  p_new = 0.5 * (iv->next->p + iv->p);
  if (_unur_FP_equal(p_new,iv->p) || _unur_FP_equal(p_new,iv->next->p)) {
    if(!(*error_count_shortinterval)){ 
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
		     "one or more intervals very short; possibly due to numerical problems with a pole");
      (*error_count_shortinterval)++;
    } 
    _unur_hinv_interval_parameter(gen,iv);
    return iv->next;
  }
  if ( (iv->next->u - iv->u > HINV_MAX_U_LENGTH) ||
       (! _unur_hinv_interval_is_monotone(gen,iv)) ) {
    iv_new = _unur_hinv_interval_new(gen,p_new,CDF(p_new));
    if (iv_new == NULL) return NULL;
    iv_new->next = iv->next;
    iv->next = iv_new;
    return iv;
  }
  _unur_hinv_interval_parameter(gen,iv);
  x = _unur_hinv_eval_polynomial( 0.5, iv->spline, GEN->order );
  Fx = CDF(x);
  if (!(fabs(Fx - 0.5*(iv->next->u + iv->u)) < GEN->u_resolution)) {
    if(fabs(p_new-x)< HINV_XDEVIATION * (iv->next->p - iv->p))
      iv_new = _unur_hinv_interval_new(gen,x,Fx);
    else
      iv_new = _unur_hinv_interval_new(gen,p_new,CDF(p_new));
    if (iv_new == NULL) return NULL;
    iv_new->next = iv->next;
    iv->next = iv_new;
    return iv;
  }
  return iv->next;
} 
int 
_unur_hinv_interval_is_monotone( struct unur_gen *gen, struct unur_hinv_interval *iv )
{
  double bound;
  switch (GEN->order) {
  case 5:
  case 3:
    if (_unur_iszero(iv->u) || _unur_FP_approx(iv->u,iv->next->u))
      return TRUE;
    bound = 3.*(iv->next->p - iv->p)/(iv->next->u - iv->u);
    return (1./iv->next->f > bound || 1./iv->f > bound) ? FALSE : TRUE;
  case 1:
  default:  
    return TRUE;
  }
} 
int
_unur_hinv_interval_parameter( struct unur_gen *gen, struct unur_hinv_interval *iv )
{
  double delta_u, delta_p;
  double f1, fs0, fs1, fss0, fss1;
  delta_u = iv->next->u - iv->u;
  delta_p = iv->next->p - iv->p;
  switch (GEN->order) {
  case 5:    
    if (iv->f > 0. && iv->next->f > 0. &&
	iv->df < INFINITY && iv->df > -INFINITY && 
	iv->next->df < INFINITY && iv->next->df > -INFINITY ) {
      f1   = delta_p;
      fs0  = delta_u / iv->f;      
      fs1  = delta_u / iv->next->f;
      fss0 = -delta_u * delta_u * iv->df / (iv->f * iv->f * iv->f);
      fss1 = -delta_u * delta_u * iv->next->df / (iv->next->f * iv->next->f * iv->next->f);
      iv->spline[0] = iv->p;
      iv->spline[1] = fs0;
      iv->spline[2] = 0.5*fss0;
      iv->spline[3] = 10.*f1 - 6.*fs0 - 4.*fs1 - 1.5*fss0 + 0.5*fss1;
      iv->spline[4] = -15.*f1 + 8.*fs0 + 7.*fs1 + 1.5*fss0 - fss1;
      iv->spline[5] = 6.*f1 - 3.*fs0 - 3.*fs1 - 0.5*fss0 + 0.5*fss1;
      return UNUR_SUCCESS;
    }
    else {
      iv->spline[4] = 0.;
      iv->spline[5] = 0.;
    }
  case 3:    
    if (iv->f > 0. && iv->next->f > 0.) {
      iv->spline[0] = iv->p;
      iv->spline[1] = delta_u / iv->f;
      iv->spline[2] = 3.* delta_p - delta_u * (2./iv->f + 1./iv->next->f);
      iv->spline[3] = -2.* delta_p + delta_u * (1./iv->f + 1./iv->next->f);
      return UNUR_SUCCESS;
    }
    else {
      iv->spline[2] = 0.;
      iv->spline[3] = 0.;
    }
  case 1:    
    iv->spline[0] = iv->p;
    iv->spline[1] = delta_p;
    return UNUR_SUCCESS;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
} 
double
_unur_hinv_eval_polynomial( double x, double *coeff, int order )
{
  int i;
  double poly;
  poly = coeff[order];
  for (i=order-1; i>=0; i--)
    poly = x*poly + coeff[i];
  return poly;
} 
int
_unur_hinv_list_to_array( struct unur_gen *gen )
{
  int i; 
  struct unur_hinv_interval *iv, *next;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,UNUR_ERR_COOKIE);
  GEN->intervals = 
    _unur_xrealloc( GEN->intervals, GEN->N*(GEN->order+2)*sizeof(double) );
  i = 0;
  for (iv=GEN->iv; iv!=NULL; iv=next) {
    GEN->intervals[i] = iv->u;
    memcpy( GEN->intervals+(i+1), &(iv->spline), (GEN->order+1)*sizeof(double) );
    i += GEN->order+2;
    next = iv->next;
    free(iv);
  }
  GEN->iv = NULL;
  return UNUR_SUCCESS;
} 
int
_unur_hinv_make_guide_table( struct unur_gen *gen )
{
  int i,j, imax;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_HINV_GEN,UNUR_ERR_COOKIE);
  GEN->guide_size = GEN->N * GEN->guide_factor;
  if (GEN->guide_size <= 0) GEN->guide_size = 1; 
  GEN->guide = _unur_xrealloc( GEN->guide, GEN->guide_size * sizeof(int) );
  imax = (GEN->N-2) * (GEN->order+2);
# define u(i)  (GEN->intervals[(i)+GEN->order+2])
  i = 0;
  GEN->guide[0] = 0;
  for( j=1; j<GEN->guide_size ;j++ ) {
    while( u(i) < (j/(double)GEN->guide_size) && i <= imax)
      i += GEN->order+2;
    if (i > imax) break;
    GEN->guide[j]=i;
  }
# undef u
  i = _unur_min(i,imax);
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = i;
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_hinv_debug_init( const struct unur_gen *gen, int ok )
{
  FILE *log;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HINV_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = HINV (Spline approximation of INVerse CDF)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(log,"%s: sampling routine = _unur_hinv_sample\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: order of polynomial = %d",gen->genid,GEN->order);
  _unur_print_if_default(gen,HINV_SET_ORDER);
  fprintf(log,"\n%s: u-resolution = %g",gen->genid,GEN->u_resolution);
  _unur_print_if_default(gen,HINV_SET_U_RESOLUTION);
  fprintf(log,"\n%s: tail cut-off points = ",gen->genid);
  if (GEN->tailcutoff_left < 0.)  fprintf(log,"none, ");
  else                            fprintf(log,"%g, ",GEN->tailcutoff_left);
  if (GEN->tailcutoff_right > 1.) fprintf(log,"none\n");
  else                            fprintf(log,"1.-%g\n",1.-GEN->tailcutoff_right);
  fprintf(log,"%s: domain of computation = [%g,%g]\n",gen->genid,GEN->bleft,GEN->bright);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);
  fprintf(log,"%s:\n",gen->genid);
  if (GEN->stp && gen->set & HINV_SET_STP) {
    fprintf(log,"%s: starting points: (%d)",gen->genid,GEN->n_stp);
    for (i=0; i<GEN->n_stp; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",GEN->stp[i]);
    }
  fprintf(log,"\n%s:\n",gen->genid);
  }
  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN->guide_factor);
  _unur_print_if_default(gen,HINV_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);
  _unur_hinv_debug_intervals(gen);
  fprintf(log,"%s: initialization %s\n",gen->genid,((ok)?"successful":"failed")); 
  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} 
void
_unur_hinv_debug_intervals( const struct unur_gen *gen )
{
  int i,n;
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HINV_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: Intervals: %d\n",gen->genid,GEN->N-1);
  if (gen->debug & HINV_DEBUG_TABLE) {
    fprintf(log,"%s:   Nr.      u=CDF(p)     p=spline[0]   spline[1]    ...\n",gen->genid);
    for (n=0; n<GEN->N-1; n++) {
      i = n*(GEN->order+2);
      fprintf(log,"%s:[%4d]: %#12.6g  %#12.6g  %#12.6g", gen->genid, n,
	      GEN->intervals[i], GEN->intervals[i+1], GEN->intervals[i+2]);
      if (GEN->order>1)
	fprintf(log,"  %#12.6g  %#12.6g", GEN->intervals[i+3], GEN->intervals[i+4]);
      if (GEN->order>3)
	fprintf(log,"  %#12.6g  %#12.6g", GEN->intervals[i+5], GEN->intervals[i+6]);
      fprintf(log,"\n");
    }
  }
  fprintf(log,"%s:\n",gen->genid);
} 
void 
_unur_hinv_debug_chg_truncated( const struct unur_gen *gen )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HINV_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(log,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);
} 
#endif   
