/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <distr/cvec.h>
#include <urng/urng.h>
#include <utils/matrix_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#include "unur_methods_source.h"
#include "arou.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "hitro.h"
#include "hitro_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define HITRO_MIN_MULTIPLIER  (1.0001)
#define HITRO_START_UVMIN  (1.e-3)
#define HITRO_DEFAULT_ADAPTIVE_MULTIPLIER  (1.1)
#define HITRO_VARMASK_VARIANT     0x000fu   
#define HITRO_VARIANT_COORD       0x0001u   
#define HITRO_VARIANT_RANDOMDIR   0x0002u   
#define HITRO_VARFLAG_ADAPTLINE   0x0010u   
#define HITRO_VARFLAG_ADAPTRECT   0x0020u   
#define HITRO_VARFLAG_BOUNDRECT   0x0040u   
#define HITRO_VARFLAG_BOUNDDOMAIN 0x0080u   
#define HITRO_SET_R          0x0001u   
#define HITRO_SET_X0         0x0002u   
#define HITRO_SET_THINNING   0x0004u   
#define HITRO_SET_BURNIN     0x0008u   
#define HITRO_SET_U          0x0010u   
#define HITRO_SET_V          0x0020u   
#define HITRO_SET_ADAPTLINE  0x0100u   
#define HITRO_SET_ADAPTRECT  0x0200u   
#define HITRO_SET_BOUNDRECT  0x0400u   
#define HITRO_SET_ADAPTMULT  0x0800u   
#define GENTYPE "HITRO"        
static struct unur_gen *_unur_hitro_init( struct unur_par *par );
static struct unur_gen *_unur_hitro_create( struct unur_par *par );
static struct unur_gen *_unur_hitro_clone( const struct unur_gen *gen );
static void _unur_hitro_free( struct unur_gen *gen);
static int _unur_hitro_coord_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_hitro_randomdir_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_hitro_rectangle( struct unur_gen *gen );
static void _unur_hitro_xy_to_vu( const struct unur_gen *gen, const double *x, double y, double *vu );
static void _unur_hitro_vu_to_x( const struct unur_gen *gen, const double *vu, double *x );
static double _unur_hitro_xv_to_u( const struct unur_gen *gen, double x, double v, int k );
static int _unur_hitro_vu_is_inside_region( const struct unur_gen *gen, const double *vu );
static struct unur_gen *_unur_hitro_normalgen( struct unur_gen *gen );
static void _unur_hitro_random_unitvector( struct unur_gen *gen, double *direction );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_hitro_debug_init_start( const struct unur_gen *gen );
static void _unur_hitro_debug_init_finished( const struct unur_gen *gen );
static void _unur_hitro_debug_free( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_hitro_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cvec      
#define PAR       ((struct unur_hitro_par*)par->datap) 
#define GEN       ((struct unur_hitro_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec 
#define SAMPLE    gen->sample.cvec           
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    
#define GEN_NORMAL    gen->gen_aux
static UNUR_SAMPLING_ROUTINE_CVEC *
_unur_hitro_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant & HITRO_VARMASK_VARIANT) {
  case HITRO_VARIANT_COORD:
    return _unur_hitro_coord_sample_cvec;
  case HITRO_VARIANT_RANDOMDIR:
  default:
    return _unur_hitro_randomdir_sample_cvec;
  }
} 
struct unur_par *
unur_hitro_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_hitro_par) );
  COOKIE_SET(par,CK_HITRO_PAR);
  par->distr    = distr;      
  par->method   = UNUR_METH_HITRO ;             
  par->variant  = ( HITRO_VARIANT_COORD |       
		    HITRO_VARFLAG_ADAPTLINE );  
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  PAR->r        = 1.;               
  PAR->thinning = 1;                
  PAR->burnin   = 0;                
  PAR->x0       = NULL;             
  PAR->adaptive_mult = HITRO_DEFAULT_ADAPTIVE_MULTIPLIER; 
  PAR->vmax     = -1.;        
  PAR->umin     = NULL;       
  PAR->umax     = NULL;       
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_hitro_init;
  return par;
} 
int
unur_hitro_set_variant_coordinate( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  par->variant = (par->variant & ~HITRO_VARMASK_VARIANT) | HITRO_VARIANT_COORD;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_variant_random_direction( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  par->variant = (par->variant & ~HITRO_VARMASK_VARIANT) | HITRO_VARIANT_RANDOMDIR;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_use_adaptiveline( struct unur_par *par, int adaptive )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  par->variant = (adaptive) 
    ? (par->variant | HITRO_VARFLAG_ADAPTLINE) 
    : (par->variant & (~HITRO_VARFLAG_ADAPTLINE));
  par->set |= HITRO_SET_ADAPTLINE;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_use_adaptiverectangle( struct unur_par *par, int adaptive )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  par->variant = (adaptive) 
    ? (par->variant | HITRO_VARFLAG_ADAPTRECT) 
    : (par->variant & (~HITRO_VARFLAG_ADAPTRECT));
  par->set |= HITRO_SET_ADAPTRECT;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_use_boundingrectangle( struct unur_par *par, int rectangle )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  par->variant = (rectangle) 
    ? (par->variant | HITRO_VARFLAG_BOUNDRECT) 
    : (par->variant & (~HITRO_VARFLAG_BOUNDRECT));
  par->set |= HITRO_SET_BOUNDRECT;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_r( struct unur_par *par, double r )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r <= 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->r = r;
  par->set |= HITRO_SET_R;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_u( struct unur_par *par, const double *umin, const double *umax )
{
  int d;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );
  for (d = 0; d < par->distr->dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
    if (! (_unur_isfinite(umax[d]) && _unur_isfinite(umin[d])) ) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"rectangle not bounded");
      return UNUR_ERR_PAR_SET;
    }
  }
  PAR->umin = umin;
  PAR->umax = umax;
  par->set |= HITRO_SET_U;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_v( struct unur_par *par, double vmax )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (! _unur_isfinite(vmax) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"rectangle not bounded");
    return UNUR_ERR_PAR_SET;
  }
  PAR->vmax = vmax;
  par->set |= HITRO_SET_V;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_adaptive_multiplier( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  if (factor < HITRO_MIN_MULTIPLIER) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"multiplier too small (<= 1.0001)");
    return UNUR_ERR_PAR_SET;
  }
  PAR->adaptive_mult = factor;
  par->set |= HITRO_SET_ADAPTMULT;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_startingpoint( struct unur_par *par, const double *x0)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  PAR->x0 = x0;
  par->set |= HITRO_SET_X0;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_thinning( struct unur_par *par, int thinning )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  if (thinning < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"thinning < 1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->thinning = thinning;
  par->set |= HITRO_SET_THINNING;
  return UNUR_SUCCESS;
} 
int
unur_hitro_set_burnin( struct unur_par *par, int burnin )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  if (burnin < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"burnin < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->burnin = burnin;
  par->set |= HITRO_SET_BURNIN;
  return UNUR_SUCCESS;
} 
const double *
unur_hitro_get_state( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, NULL );
  if (gen->method != UNUR_METH_HITRO) {
    _unur_error(gen->genid, UNUR_ERR_GEN_INVALID,"");
    return NULL;
  }
  return GEN->state;
} 
int
unur_hitro_chg_state( struct unur_gen *gen, const double *state )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, HITRO, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( gen->genid, state, UNUR_ERR_NULL );
  if ( ! _unur_hitro_vu_is_inside_region(gen,state) ) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"invalid state");
    return UNUR_ERR_PAR_SET;
  }
  memcpy( GEN->state, state, GEN->dim * sizeof(double));
  return UNUR_SUCCESS;
} 
int
unur_hitro_reset_state( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, HITRO, UNUR_ERR_GEN_INVALID );
  memcpy( GEN->state, GEN->x0, GEN->dim * sizeof(double));
  _unur_hitro_xy_to_vu(gen, GEN->x0, GEN->fx0/2., GEN->state );
  memcpy( GEN->vu, GEN->state, (GEN->dim + 1) * sizeof(double) );
  GEN->vumax[0] = pow(GEN->fx0, 1./(GEN->r * GEN->dim + 1.)) * (1. + DBL_EPSILON); 
  if (gen->variant & HITRO_VARIANT_COORD) GEN->coord = 0;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hitro_init( struct unur_par *par )
{
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_HITRO ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HITRO_PAR,NULL);
  if ( par->variant & HITRO_VARIANT_COORD ) {
    if (_unur_distr_cvec_has_boundeddomain(par->distr))
      par->variant |= HITRO_VARFLAG_BOUNDDOMAIN;
    else
      par->variant |= HITRO_VARFLAG_BOUNDRECT;
    if (!(par->set & HITRO_SET_ADAPTRECT) )
      par->variant |= HITRO_VARFLAG_ADAPTRECT;
  }
  gen = _unur_hitro_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_hitro_debug_init_start(gen);
#endif
  GEN->fx0 = PDF(GEN->x0);
  if ( (GEN->fx0 / 2.) <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"x0 not in support of PDF");
    _unur_hitro_free(gen); return NULL;
  }
  _unur_hitro_xy_to_vu(gen, GEN->x0, GEN->fx0/2., GEN->state );
  memcpy( GEN->vu, GEN->state, (GEN->dim + 1) * sizeof(double) );
  GEN->vumax[0] = pow(GEN->fx0, 1./(GEN->r * GEN->dim + 1.)) * (1. + DBL_EPSILON); 
  if (gen->variant & HITRO_VARIANT_RANDOMDIR ) {
    GEN_NORMAL = _unur_hitro_normalgen( gen );
    if ( GEN_NORMAL == NULL ) {
      _unur_hitro_free(gen); return NULL;
    }
  }
  if ( !(gen->variant & HITRO_VARFLAG_ADAPTRECT) )
    if (_unur_hitro_rectangle(gen) != UNUR_SUCCESS ) {
      _unur_hitro_free(gen); return NULL;
    }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_hitro_debug_init_finished(gen);
#endif
  if (GEN->burnin > 0 ) {
    int thinning, burnin;
    double *X;
    X = _unur_xmalloc( GEN->dim * sizeof(double) );
    thinning = GEN->thinning;
    GEN->thinning = 1;
    for (burnin = GEN->burnin; burnin>0; --burnin)
      _unur_sample_vec(gen,X);
    GEN->thinning = thinning;
    free (X);
  }
  gen->status = UNUR_SUCCESS;
  return gen;
} 
static struct unur_gen *
_unur_hitro_create( struct unur_par *par )
{
  struct unur_gen *gen;
  int i;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HITRO_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_hitro_gen) );
  COOKIE_SET(gen,CK_HITRO_GEN);
  GEN->dim = gen->distr->dim;
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_hitro_getSAMPLE(gen);
  gen->destroy = _unur_hitro_free;
  gen->clone = _unur_hitro_clone;
  gen->variant = par->variant;
  GEN->thinning = PAR->thinning; 
  GEN->burnin = PAR->burnin;     
  GEN->r = PAR->r;               
  GEN->adaptive_mult = PAR->adaptive_mult; 
  GEN->center = unur_distr_cvec_get_center(gen->distr);
  GEN->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
  if (PAR->x0 == NULL)
    PAR->x0 = unur_distr_cvec_get_center(gen->distr);
  memcpy( GEN->x0, PAR->x0, GEN->dim * sizeof(double));
  GEN->vumin = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
  GEN->vumax = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
  GEN->vumin[0] = 0.;
  GEN->vumax[0] = (PAR->vmax > 0.) ? PAR->vmax : HITRO_START_UVMIN;  
  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    if (PAR->umin && PAR->umax) {
      memcpy (GEN->vumin+1, PAR->umin, GEN->dim * sizeof(double) );
      memcpy (GEN->vumax+1, PAR->umax, GEN->dim * sizeof(double) );
    }
    else {
      for (i=1; i<GEN->dim+1; i++) GEN->vumin[i] = -HITRO_START_UVMIN;
      for (i=1; i<GEN->dim+1; i++) GEN->vumax[i] =  HITRO_START_UVMIN;
    }
  }
  GEN->state = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
  GEN->x     = _unur_xmalloc( GEN->dim * sizeof(double) );
  GEN->vu    = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
  GEN->direction = _unur_xmalloc( (1 + GEN->dim) * sizeof(double));
  GEN->coord = 0;           
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_hitro_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_hitro_clone( const struct unur_gen *gen )
{
#define CLONE         ((struct unur_hitro_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HITRO_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->center = unur_distr_cvec_get_center(clone->distr);
  if (GEN->state) {
    CLONE->state = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
    memcpy( CLONE->state, GEN->state, (1 + GEN->dim) * sizeof(double) );
  }
  if (GEN->vumin) {
    CLONE->vumin = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
    memcpy( CLONE->vumin, GEN->vumin, (GEN->dim+1) * sizeof(double) );
  }
  if (GEN->vumax) {
    CLONE->vumax = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
    memcpy( CLONE->vumax, GEN->vumax, (GEN->dim+1) * sizeof(double) );
  }
  if (GEN->x0) {
    CLONE->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
    memcpy( CLONE->x0, GEN->x0, GEN->dim * sizeof(double));
  }
  if (GEN->x) {
    CLONE->x = _unur_xmalloc( GEN->dim * sizeof(double) );
    memcpy( CLONE->x, GEN->x, GEN->dim * sizeof(double) );
  }
  if (GEN->vu) {
    CLONE->vu = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
    memcpy( CLONE->vu, GEN->vu, (1 + GEN->dim) * sizeof(double) );
  }
  if (GEN->direction) {
    CLONE->direction = _unur_xmalloc( (1 + GEN->dim) * sizeof(double));
    memcpy( CLONE->direction, GEN->direction, (1 + GEN->dim) * sizeof(double));
  }
  return clone;
#undef CLONE
} 
void
_unur_hitro_free( struct unur_gen *gen )
{
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_HITRO ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);
  SAMPLE = NULL;   
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_hitro_debug_free(gen);
#endif
  if (GEN->state) free (GEN->state);
  if (GEN->x0) free (GEN->x0);
  if (GEN->x) free (GEN->x);
  if (GEN->vu) free (GEN->vu);
  if (GEN->direction) free (GEN->direction);
  if (GEN->vumin) free (GEN->vumin);
  if (GEN->vumax) free (GEN->vumax);
  _unur_generic_free(gen);
} 
int
_unur_hitro_coord_sample_cvec( struct unur_gen *gen, double *vec )
{
  int thinning;
  double lmin, lmax, lmid;   
  double *vuaux;  
  int coord;      
  double U;       
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_HITRO_GEN,UNUR_ERR_COOKIE);
  vuaux = GEN->vu;
  for (thinning = GEN->thinning; thinning > 0; --thinning) {
    coord = GEN->coord = (GEN->coord + 1) % (GEN->dim + 1);
    if (! (gen->variant & HITRO_VARFLAG_BOUNDDOMAIN) || coord == 0) {
      lmin = GEN->vumin[coord];
      lmax = GEN->vumax[coord];
    }
    else {
      int k = coord-1;
      double *domain = DISTR.domainrect;
      lmin = _unur_hitro_xv_to_u(gen, domain[2*k], vuaux[0], k );
      lmax = _unur_hitro_xv_to_u(gen, domain[2*k+1], vuaux[0], k );
      if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
	lmin = _unur_max(lmin,GEN->vumin[coord]);
	lmax = _unur_min(lmax,GEN->vumax[coord]);
      }
    }
    if ( gen->variant & HITRO_VARFLAG_ADAPTRECT ) {
      lmid = 0.5 * (lmin + lmax);
      vuaux[coord] = lmax;
      while ( _unur_hitro_vu_is_inside_region(gen,vuaux) ) {
	lmax = lmid + (lmax-lmid) * GEN->adaptive_mult;
	GEN->vumax[coord] = vuaux[coord] = lmax;
      }
      vuaux[coord] = lmin;
      while ( coord!=0 && _unur_hitro_vu_is_inside_region(gen,vuaux) ) {
	lmin = lmid + (lmin-lmid) * GEN->adaptive_mult;
	GEN->vumin[coord] = vuaux[coord] = lmin;
      }
    }
    while (1) {
      U = _unur_call_urng(gen->urng);
      vuaux[coord] = U * lmin + (1.-U) * lmax;
      if (_unur_hitro_vu_is_inside_region(gen,vuaux) )
	break;
      if ( gen->variant & HITRO_VARFLAG_ADAPTLINE ) {
	if (GEN->state[coord] < vuaux[coord])
	  lmax = vuaux[coord];
	else
	  lmin = vuaux[coord];
      }
    }
    GEN->state[coord] = vuaux[coord];
  }
  _unur_hitro_vu_to_x( gen, GEN->state, vec );
  return UNUR_SUCCESS;
} 
int
_unur_hitro_randomdir_sample_cvec( struct unur_gen *gen, double *vec )
{
#define new_point(ll)  { int j; for (j=0;j<dim+1;j++) vuaux[j] = GEN->state[j]+(ll)*GEN->direction[j]; }
  int thinning;
  int i, d, k;
  double lambda, lb[2];  
  double *vuaux;       
  double U;
  int update;
  int dim = GEN->dim;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_HITRO_GEN,UNUR_ERR_COOKIE);
  d = (gen->variant & HITRO_VARFLAG_BOUNDRECT) ? dim+1 : 1;
  vuaux = GEN->vu;
  for (thinning = GEN->thinning; thinning > 0; --thinning) {
    _unur_hitro_random_unitvector( gen, GEN->direction );
    lb[1] = UNUR_INFINITY;  
    lb[0] = -UNUR_INFINITY;
    for (i=0; i<d; i++) {
      lambda = (GEN->vumin[i] - GEN->state[i]) / GEN->direction[i];
      if (lambda>0 && lambda<lb[1]) lb[1] = lambda;
      if (lambda<0 && lambda>lb[0]) lb[0] = lambda;
      lambda = (GEN->vumax[i] - GEN->state[i]) / GEN->direction[i];
      if (lambda>0 && lambda<lb[1]) lb[1] = lambda;
      if (lambda<0 && lambda>lb[0]) lb[0] = lambda;
    }
    if (! (_unur_isfinite(lb[0]) && _unur_isfinite(lb[1])) ) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"line segment not bounded, try again");
      continue;
    }
    if ( gen->variant & HITRO_VARFLAG_ADAPTRECT ) {
      for (k=0; k<2; k++) {
	update = FALSE;
	while (1) {
	  new_point(lb[k]);
	  if (! _unur_hitro_vu_is_inside_region(gen,vuaux) )
	    break;
	  update = TRUE;
	  lb[k] *= GEN->adaptive_mult;
	}
	if (update) {
	  new_point(lb[k]);
	  for (i=0; i<d; i++) {
	    if (vuaux[i] < GEN->vumin[i] && i!=0) GEN->vumin[i] =  vuaux[i];
	    if (vuaux[i] > GEN->vumax[i])         GEN->vumax[i] =  vuaux[i];
	  }			
	}
      }
    }		
    while (1) {
      U = _unur_call_urng(gen->urng); 
      lambda = U * lb[0] + (1.-U) * lb[1];
      new_point(lambda);
      if (_unur_hitro_vu_is_inside_region(gen,vuaux) )
	break;
      if ( gen->variant & HITRO_VARFLAG_ADAPTLINE ) {
	if (lambda < 0) lb[0] = lambda;
	else            lb[1] = lambda;
      }
    }
    memcpy( GEN->state, vuaux, (dim+1)*sizeof(double) );
  }
  _unur_hitro_vu_to_x( gen, GEN->state, vec );
  return UNUR_SUCCESS;
#undef new_point
} 
int
_unur_hitro_rectangle( struct unur_gen *gen )
{
  int d; 
  struct MROU_RECTANGLE *rr;
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_HITRO_GEN, UNUR_ERR_COOKIE );
  if ((gen->set & HITRO_SET_U) && (gen->set & HITRO_SET_V))
    return UNUR_SUCCESS;
  rr = _unur_mrou_rectangle_new();
  rr->distr  = gen->distr;
  rr->dim    = GEN->dim;
  rr->umin   = GEN->vumin+1;
  rr->umax   = GEN->vumax+1;
  rr->r      = GEN->r;
  rr->center = GEN->center;
  rr->genid  = gen->genid;
  rr->bounding_rectangle = 
    ( (gen->variant & HITRO_VARFLAG_BOUNDRECT) && !(gen->set & HITRO_SET_U) )
    ? 1 : 0;
  if ( _unur_mrou_rectangle_compute(rr) != UNUR_SUCCESS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot compute bounding rectangle, try adaptive");
    gen->variant &= HITRO_VARFLAG_ADAPTRECT;
    free(rr); return UNUR_ERR_GEN_CONDITION;
  }
  if (!(gen->set & HITRO_SET_V)) {
    GEN->vumax[0] = rr->vmax;
  }
  if (rr->bounding_rectangle) {
    for (d=0; d<GEN->dim; d++) GEN->vumin[d+1] = rr->umin[d];
    for (d=0; d<GEN->dim; d++) GEN->vumax[d+1] = rr->umax[d];
  }
  free(rr);
  return UNUR_SUCCESS;
} 
double
_unur_hitro_xv_to_u( const struct unur_gen *gen, double x, double v, int k )
{
  if (_unur_isone(GEN->r)) 
    return (x - GEN->center[k]) * v;
  else
    return (x - GEN->center[k]) * pow(v,GEN->r) ;
} 
void 
_unur_hitro_xy_to_vu( const struct unur_gen *gen, const double *x, double y, double *vu )
{
  int d;
  double v;
  double *u = vu+1;
  vu[0] = v = pow(y, 1./(GEN->r * GEN->dim + 1.));
  if (_unur_isone(GEN->r)) 
    for (d=0; d<GEN->dim; d++)  u[d] = (x[d] - GEN->center[d]) * v;
  else
    for (d=0; d<GEN->dim; d++)  u[d] = (x[d] - GEN->center[d]) * pow(v,GEN->r) ;
} 
void 
_unur_hitro_vu_to_x( const struct unur_gen *gen, const double *vu, double *x )
{
  int d;
  double v = vu[0];
  const double *u = vu+1;
  if (v<=0.) {
    for (d=0; d<GEN->dim; d++)  x[d] = 0.;
    return;
  }
  if (_unur_isone(GEN->r))
    for (d=0; d<GEN->dim; d++)  x[d] = u[d]/v + GEN->center[d];
  else
    for (d=0; d<GEN->dim; d++)  x[d] = u[d]/pow(v,GEN->r) + GEN->center[d];
} 
int
_unur_hitro_vu_is_inside_region( const struct unur_gen *gen, const double *vu )
{
  double y;
  double v = vu[0];
  _unur_hitro_vu_to_x( gen, vu, GEN->x );
  y = PDF(GEN->x);
  if (y <= 0. || v <= 0.) return FALSE;
  return ( (v < pow(y,1./(GEN->r * GEN->dim + 1.))) ? TRUE : FALSE );
} 
struct unur_gen *
_unur_hitro_normalgen( struct unur_gen *gen )
{
  struct unur_gen   *normalgen;
  struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
  struct unur_par   *normalpar = unur_arou_new( normaldistr );
  unur_arou_set_usedars( normalpar, TRUE );
  normalgen = unur_init( normalpar );
  _unur_distr_free( normaldistr );
  if (normalgen == NULL) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,
		"Cannot create aux Gaussian generator");
    return NULL;
  }
  normalgen->urng = gen->urng;
  normalgen->debug = gen->debug;
  return normalgen;
} 
void
_unur_hitro_random_unitvector( struct unur_gen *gen, double *direction )
{
  int i;
  do {
    for (i=0; i<GEN->dim+1; i++)
      direction[i] = unur_sample_cont(GEN_NORMAL);
    _unur_vector_normalize(GEN->dim+1, direction);
  } while (!_unur_isfinite(direction[0]));
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_hitro_debug_init_start( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = HITRO (Markov Chain - HITRO sampler)\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  switch (gen->variant & HITRO_VARMASK_VARIANT) {
  case HITRO_VARIANT_COORD:
    fprintf(LOG,"coordinate sampling [default]\n"); break;
  case HITRO_VARIANT_RANDOMDIR:
    fprintf(LOG,"random direction sampling\n"); break;
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: r = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,HITRO_SET_R);
  fprintf(LOG,"\n%s: adaptive line sampling: %s",gen->genid,
	  (gen->variant&HITRO_VARFLAG_ADAPTLINE)?"on":"off");
  _unur_print_if_default(gen,HITRO_SET_ADAPTLINE);
  fprintf(LOG,"\n%s: use entire bounding rectangle: %s",gen->genid,
	  (gen->variant&HITRO_VARFLAG_BOUNDRECT)?"on":"off");
  _unur_print_if_default(gen,HITRO_SET_BOUNDRECT);
  fprintf(LOG,"\n%s: adaptive bounding rectangle: %s",gen->genid,
	  (gen->variant&HITRO_VARFLAG_ADAPTRECT)?"on":"off");
  _unur_print_if_default(gen,HITRO_SET_ADAPTRECT);
  if (gen->variant&HITRO_VARFLAG_ADAPTRECT) {
    fprintf(LOG,"\n%s:\tmultiplier = %g",gen->genid,GEN->adaptive_mult);
    _unur_print_if_default(gen,HITRO_SET_ADAPTMULT);
  }
  fprintf(LOG,"\n%s: use domain of distribution: %s\n",gen->genid,
	  (gen->variant&HITRO_VARFLAG_BOUNDDOMAIN)?"on":"off");
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cvec_debug( gen->distr, gen->genid );
  switch (gen->variant & HITRO_VARMASK_VARIANT) {
  case HITRO_VARIANT_COORD:
    fprintf(LOG,"%s: sampling routine = _unur_hitro_coord_sample_cvec()\n",gen->genid);
    break;
  case HITRO_VARIANT_RANDOMDIR:
    fprintf(LOG,"%s: sampling routine = _unur_hitro_randomdir_sample_cvec()\n",gen->genid);
    break;
  }
  fprintf(LOG,"%s: thinning = %d",gen->genid,GEN->thinning);
  _unur_print_if_default(gen,HITRO_SET_THINNING);
  fprintf(LOG,"\n%s: burn-in = %d",gen->genid,GEN->burnin);
  _unur_print_if_default(gen,HITRO_SET_BURNIN);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_matrix_print_vector( GEN->dim, GEN->x0, "starting point = ", LOG, gen->genid, "\t   ");
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_hitro_debug_init_finished( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    fprintf(LOG,"%s: bounding rectangle%s:\n",gen->genid,
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? " [start for adaptive rectangle]" : "" );
    fprintf(LOG,"%s: vmax = %g\n",gen->genid, GEN->vumax[0]);
    _unur_matrix_print_vector( GEN->dim, GEN->vumin+1, "umin =", LOG, gen->genid, "\t   ");
    _unur_matrix_print_vector( GEN->dim, GEN->vumax+1, "umax =", LOG, gen->genid, "\t   ");
  }
  else {
    fprintf(LOG,"%s: upper bound vmax = %g %s\n",gen->genid, GEN->vumax[0],
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[start for adaptive bound]" : "" );
  }
  _unur_matrix_print_vector( GEN->dim+1, GEN->state, "starting state = ", LOG, gen->genid, "\t   ");
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fflush(LOG);
} 
void
_unur_hitro_debug_free( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->status == UNUR_SUCCESS) {
    fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
  }
  else {
    fprintf(LOG,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    fprintf(LOG,"%s: bounding rectangle%s:\n",gen->genid,
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? " [adaptive]" : "" );
    fprintf(LOG,"%s: vmax = %g\n",gen->genid, GEN->vumax[0]);
    _unur_matrix_print_vector( GEN->dim, GEN->vumin+1, "umin =", LOG, gen->genid, "\t   ");
    _unur_matrix_print_vector( GEN->dim, GEN->vumax+1, "umax =", LOG, gen->genid, "\t   ");
  }
  else {
    fprintf(LOG,"%s: upper bound vmax = %g %s\n",gen->genid, GEN->vumax[0],
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[adaptive]" : "" );
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_hitro_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  int i;
  double rc;
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
  _unur_string_append(info,"method: HITRO (HIT-and-run sampler with Ratio-Of-uniforms [MCMC])\n");
  _unur_string_append(info,"   variant = %s\n",
		      ((gen->variant & HITRO_VARMASK_VARIANT)==HITRO_VARIANT_COORD)
		      ? "coordinate sampling [default]" : "random direction sampling");
  _unur_string_append(info,"   r = %g\n", GEN->r);
  _unur_string_append(info,"   thinning = %d\n", GEN->thinning);
  _unur_string_append(info,"   adaptive line sampling = %s\n", 
		      (gen->variant&HITRO_VARFLAG_ADAPTLINE)?"on":"off");
  _unur_string_append(info,"   use entire bounding rectangle = %s\n",
		      (gen->variant&HITRO_VARFLAG_BOUNDRECT)?"on":"off");
  if (gen->variant&HITRO_VARFLAG_ADAPTRECT)
    _unur_string_append(info,"   adaptive bounding rectangle = on  [multiplier = %g]\n",
			GEN->adaptive_mult);
  else
    _unur_string_append(info,"   adaptive bounding rectangle = off\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  rc = unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize);
  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    _unur_string_append(info,"   bounding rectangle %s= ",
			(gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[adaptive] " : "" );
    for (i=0; i<GEN->dim; i++)
      _unur_string_append(info,"%s(%g,%g)", i?"x":"", GEN->vumin[i+1], GEN->vumax[i+1]);
    _unur_string_append(info," x (0,%g)\n", GEN->vumax[0]);
  }
  else {
    _unur_string_append(info,"   upper bound vmax = %g %s\n", GEN->vumax[0],
			(gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[adaptive]" : "" );
  }
  _unur_string_append(info,"   rejection constant =  %.2f  [approx.]\n", rc);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    switch (gen->variant & HITRO_VARMASK_VARIANT) {
    case HITRO_VARIANT_COORD:
      _unur_string_append(info,"   variant_coordinate  [default]\n"); break;
    case HITRO_VARIANT_RANDOMDIR:
      _unur_string_append(info,"   variant_random_direction\n"); break;
    }
    _unur_string_append(info,"   r = %g  %s\n", GEN->r,
 			(gen->set & HITRO_SET_R) ? "" : "[default]");
    _unur_string_append(info,"   adaptiveline = %s  %s\n", 
			(gen->variant&HITRO_VARFLAG_ADAPTLINE)?"on":"off",
 			(gen->set & HITRO_SET_ADAPTLINE) ? "" : "[default]");
    _unur_string_append(info,"   boundingrectangle = %s  %s\n",
			(gen->variant&HITRO_VARFLAG_BOUNDRECT)?"on":"off",
 			(gen->set & HITRO_SET_BOUNDRECT) ? "" : "[default]");
    _unur_string_append(info,"   adaptiverectangle = %s  %s\n", 
			(gen->variant&HITRO_VARFLAG_ADAPTRECT)?"on":"off",
 			(gen->set & HITRO_SET_ADAPTRECT) ? "" : "[default]");
    if (gen->variant&HITRO_VARFLAG_ADAPTRECT)
      _unur_string_append(info,"   adaptive_multiplier = %g  %s\n", 
			  GEN->adaptive_mult,
			  (gen->set & HITRO_SET_ADAPTMULT) ? "" : "[default]");
   _unur_string_append(info,"   thinning = %d  %s\n", GEN->thinning,
 			(gen->set & HITRO_SET_THINNING) ? "" : "[default]");
   _unur_string_append(info,"   burnin = %d  %s\n", GEN->burnin,
 			(gen->set & HITRO_SET_THINNING) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }
} 
#endif   
