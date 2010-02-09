/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <distr/condi.h>
#include <distr/cvec.h>
#include <urng/urng.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "ars.h"
#include "tdr.h"
#include "gibbs.h"
#include "gibbs_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define GIBBS_VARMASK_VARIANT     0x000fu    
#define GIBBS_VARIANT_COORD       0x0001u    
#define GIBBS_VARIANT_RANDOMDIR   0x0002u    
#define GIBBS_VARMASK_T           0x00f0u    
#define GIBBS_VAR_T_SQRT          0x0010u    
#define GIBBS_VAR_T_LOG           0x0020u    
#define GIBBS_VAR_T_POW           0x0030u    
#define GIBBS_DEBUG_CONDI   0x01000000u
#define GIBBS_SET_C          0x001u    
#define GIBBS_SET_X0         0x002u    
#define GIBBS_SET_THINNING   0x004u    
#define GIBBS_SET_BURNIN     0x008u    
#define GENTYPE "GIBBS"        
static struct unur_gen *_unur_gibbs_init( struct unur_par *par );
static int _unur_gibbs_coord_init( struct unur_gen *gen );
static int _unur_gibbs_randomdir_init( struct unur_gen *gen );
static struct unur_gen *_unur_gibbs_create( struct unur_par *par );
static struct unur_gen *_unur_gibbs_clone( const struct unur_gen *gen );
static void _unur_gibbs_free( struct unur_gen *gen);
static int _unur_gibbs_coord_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_gibbs_randomdir_sample_cvec( struct unur_gen *gen, double *vec );
static struct unur_gen *_unur_gibbs_normalgen( struct unur_gen *gen );
static void _unur_gibbs_random_unitvector( struct unur_gen *gen, double *direction );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_gibbs_debug_init_start( const struct unur_gen *gen );
static void _unur_gibbs_debug_init_condi( const struct unur_gen *gen );
static void _unur_gibbs_debug_burnin_failed( const struct unur_gen *gen );
static void _unur_gibbs_debug_init_finished( const struct unur_gen *gen, int success );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_gibbs_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cvec      
#define PAR       ((struct unur_gibbs_par*)par->datap) 
#define GEN       ((struct unur_gibbs_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec 
#define SAMPLE    gen->sample.cvec           
#define GEN_CONDI     gen->gen_aux_list     
#define GEN_NORMAL    gen->gen_aux
static UNUR_SAMPLING_ROUTINE_CVEC *
_unur_gibbs_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_RANDOMDIR:
    return _unur_gibbs_randomdir_sample_cvec;
  case GIBBS_VARIANT_COORD:
  default:
    return _unur_gibbs_coord_sample_cvec;
  }
} 
struct unur_par *
unur_gibbs_new( const struct unur_distr *distr )
{
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);
  if (DISTR_IN.logpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"logPDF");
    return NULL;
  }
  if (DISTR_IN.dlogpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dlogPDF");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_gibbs_par) );
  COOKIE_SET(par,CK_GIBBS_PAR);
  par->distr    = distr;      
  PAR->c_T      = 0.;        
  par->method   = UNUR_METH_GIBBS ;       
  par->variant  = GIBBS_VARIANT_COORD;    
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  PAR->thinning = 1;                
  PAR->burnin   = 0;                
  PAR->x0       = NULL;             
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_gibbs_init;
  return par;
} 
int
unur_gibbs_set_variant_coordinate( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  par->variant = (par->variant & ~GIBBS_VARMASK_VARIANT) | GIBBS_VARIANT_COORD;
  return UNUR_SUCCESS;
} 
int
unur_gibbs_set_variant_random_direction( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  par->variant = (par->variant & ~GIBBS_VARMASK_VARIANT) | GIBBS_VARIANT_RANDOMDIR;
  return UNUR_SUCCESS;
} 
int
unur_gibbs_set_c( struct unur_par *par, double c )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  if (c > 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c > 0");
    return UNUR_ERR_PAR_SET;
  }
  if (c < -0.5) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"c < -0.5 not implemented yet");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_iszero(c) && c > -0.5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"-0.5 < c < 0 not recommended. using c = -0.5 instead.");
    c = -0.5;
  }
  PAR->c_T = c;
  par->set |= GIBBS_SET_C;
  return UNUR_SUCCESS;
} 
int 
unur_gibbs_set_startingpoint( struct unur_par *par, const double *x0)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  PAR->x0 = x0;
  par->set |= GIBBS_SET_X0;
  return UNUR_SUCCESS;
} 
int
unur_gibbs_set_thinning( struct unur_par *par, int thinning )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  if (thinning < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"thinning < 1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->thinning = thinning;
  par->set |= GIBBS_SET_THINNING;
  return UNUR_SUCCESS;
} 
int
unur_gibbs_set_burnin( struct unur_par *par, int burnin )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  if (burnin < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"burnin < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->burnin = burnin;
  par->set |= GIBBS_SET_BURNIN;
  return UNUR_SUCCESS;
} 
const double *
unur_gibbs_get_state( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, NULL );
  if (gen->method != UNUR_METH_GIBBS) {
    _unur_error(gen->genid, UNUR_ERR_GEN_INVALID,"");
    return NULL;
  }
  return GEN->state;
} 
int 
unur_gibbs_chg_state( struct unur_gen *gen, const double *state )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, GIBBS, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( gen->genid, state, UNUR_ERR_NULL );
  memcpy( GEN->state, state, GEN->dim * sizeof(double));
  return UNUR_SUCCESS;
} 
int 
unur_gibbs_reset_state( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, GIBBS, UNUR_ERR_GEN_INVALID );
  memcpy( GEN->state, GEN->x0, GEN->dim * sizeof(double));
  if (gen->variant & GIBBS_VARIANT_COORD)
    GEN->coord = (GEN->dim)-1;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_gibbs_init( struct unur_par *par )
{
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_GIBBS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_GIBBS_PAR,NULL);
  gen = _unur_gibbs_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_gibbs_debug_init_start(gen);
#endif
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:  
    if (_unur_gibbs_coord_init(gen)!=UNUR_SUCCESS) {
      _unur_gibbs_free(gen); return NULL;
    }
    break;
  case GIBBS_VARIANT_RANDOMDIR:  
    if (_unur_gibbs_randomdir_init(gen)!=UNUR_SUCCESS) {
      _unur_gibbs_free(gen); return NULL;
    }
    break;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_gibbs_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_gibbs_debug_init_condi(gen);
#endif
  if (GEN->burnin > 0 ) {
    int thinning, burnin;
    double *X;
    X = _unur_xmalloc( GEN->dim * sizeof(double) );
    thinning = GEN->thinning;
    GEN->thinning = 1;
    for (burnin = GEN->burnin; burnin>0; --burnin) {
      if ( _unur_sample_vec(gen,X) != UNUR_SUCCESS ) {
#ifdef UNUR_ENABLE_LOGGING
	_unur_gibbs_debug_burnin_failed(gen);
	if (gen->debug) _unur_gibbs_debug_init_finished(gen,FALSE);
#endif
	_unur_gibbs_free(gen); free (X); return NULL;
      }
    }
    GEN->thinning = thinning;
    free (X);
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_gibbs_debug_init_finished(gen,TRUE);
#endif
  return gen;
} 
int
_unur_gibbs_coord_init( struct unur_gen *gen )
{
  struct unur_par *par_condi;
  struct unur_gen *gen_condi;
  int i;
  int errorcode = UNUR_SUCCESS;
  GEN->distr_condi = unur_distr_condi_new( gen->distr, GEN->state, NULL, 0);
  for (i=0; i<GEN->dim; i++) {
    if ( (errorcode=unur_distr_condi_set_condition(GEN->distr_condi,GEN->state,NULL,i))
	 != UNUR_SUCCESS )
      break;
    switch( gen->variant & GIBBS_VARMASK_T ) {
    case GIBBS_VAR_T_LOG:
      par_condi = unur_ars_new(GEN->distr_condi);
      unur_ars_set_reinit_percentiles(par_condi,2,NULL);
      break;
    case GIBBS_VAR_T_SQRT:
      par_condi = unur_tdr_new(GEN->distr_condi);
      unur_tdr_set_reinit_percentiles(par_condi,2,NULL);
      unur_tdr_set_c(par_condi,-0.5);
      unur_tdr_set_usedars(par_condi,FALSE);
      unur_tdr_set_variant_gw(par_condi);
      break;
    case GIBBS_VAR_T_POW:
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    unur_set_use_distr_privatecopy( par_condi, FALSE );
    unur_set_debug( par_condi, (gen->debug&GIBBS_DEBUG_CONDI)?gen->debug:(gen->debug?1u:0u));
    unur_set_urng( par_condi, gen->urng );
    gen_condi = unur_init(par_condi);
    if (gen_condi == NULL) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "Cannot create generator for conditional distributions");
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug) _unur_gibbs_debug_init_finished(gen,FALSE);
#endif
      errorcode = UNUR_ERR_GEN_CONDITION;
      break;
    }
    GEN_CONDI[i] = gen_condi;
    if (i==0 && DISTR.domainrect==NULL) {
      for (i=1; i<GEN->dim; i++)
	GEN_CONDI[i] = unur_gen_clone(gen_condi);
      break;
    }
  }
  return errorcode;
} 
int
_unur_gibbs_randomdir_init( struct unur_gen *gen )
{
  struct unur_par *par_condi;
  struct unur_gen *gen_condi;
  GEN_NORMAL = _unur_gibbs_normalgen( gen );
  if ( GEN_NORMAL == NULL ) return UNUR_FAILURE;
  _unur_gibbs_random_unitvector( gen, GEN->direction );
  GEN->distr_condi = unur_distr_condi_new( gen->distr, GEN->state, GEN->direction, 0);
  switch( gen->variant & GIBBS_VARMASK_T ) {
  case GIBBS_VAR_T_LOG:
    par_condi = unur_ars_new(GEN->distr_condi);
    unur_ars_set_reinit_percentiles(par_condi,2,NULL);
    break;
  case GIBBS_VAR_T_SQRT:
    par_condi = unur_tdr_new(GEN->distr_condi);
    unur_tdr_set_reinit_percentiles(par_condi,2,NULL);
    unur_tdr_set_c(par_condi,-0.5);
    unur_tdr_set_usedars(par_condi,FALSE);
    break;
  case GIBBS_VAR_T_POW:
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  unur_set_use_distr_privatecopy( par_condi, FALSE );
  unur_set_debug( par_condi, (gen->debug&GIBBS_DEBUG_CONDI)?gen->debug:(gen->debug?1u:0u));
  unur_set_urng( par_condi, gen->urng );
  gen_condi = unur_init(par_condi);
  if (gen_condi == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"Cannot create generator for conditional distributions");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_gibbs_debug_init_finished(gen,FALSE);
#endif
    return UNUR_ERR_GEN_CONDITION;
  }
  *GEN_CONDI = gen_condi;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_gibbs_create( struct unur_par *par )
{
  struct unur_gen *gen;
  int i;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_GIBBS_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_gibbs_gen) );
  COOKIE_SET(gen,CK_GIBBS_GEN);
  GEN->dim = gen->distr->dim;
  gen->genid = _unur_set_genid(GENTYPE);
  if ( _unur_iszero(PAR->c_T) )
    par->variant = (par->variant & (~GIBBS_VARMASK_T)) | GIBBS_VAR_T_LOG;
  else if (_unur_FP_same(PAR->c_T, -0.5))
    par->variant = (par->variant & (~GIBBS_VARMASK_T)) | GIBBS_VAR_T_SQRT;
  else
    par->variant = (par->variant & (~GIBBS_VARMASK_T)) | GIBBS_VAR_T_POW;
  SAMPLE = _unur_gibbs_getSAMPLE(gen);
  gen->destroy = _unur_gibbs_free;
  gen->clone = _unur_gibbs_clone;
  gen->variant = par->variant;        
  GEN->thinning = PAR->thinning;           
  GEN->burnin = PAR->burnin;               
  GEN->c_T = PAR->c_T;                     
  GEN->state = _unur_xmalloc( GEN->dim * sizeof(double));
  GEN->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
  if (PAR->x0 == NULL) 
    PAR->x0 = unur_distr_cvec_get_center(gen->distr);
  memcpy( GEN->state, PAR->x0, GEN->dim * sizeof(double));
  memcpy( GEN->x0, PAR->x0, GEN->dim * sizeof(double));
  GEN->distr_condi = NULL;
  GEN_CONDI = _unur_xmalloc( GEN->dim * sizeof(struct unur_gen *) );
  gen->n_gen_aux_list = GEN->dim;   
  for (i=0; i<GEN->dim; i++) GEN_CONDI[i] = NULL;
  GEN->direction = _unur_xmalloc( GEN->dim * sizeof(double));
  GEN->coord = (GEN->dim)-1;      
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_gibbs_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_gibbs_clone( const struct unur_gen *gen )
{
#define CLONE         ((struct unur_gibbs_gen*)clone->datap)
#define CLONE_CONDI   clone->gen_aux_list     
  int i;
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_GIBBS_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->state = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy( CLONE->state, GEN->state, GEN->dim * sizeof(double));
  CLONE->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy( CLONE->x0, GEN->x0, GEN->dim * sizeof(double));
  if (GEN->distr_condi) CLONE->distr_condi = _unur_distr_clone( GEN->distr_condi );
  if (CLONE_CONDI) {
    for (i=0; i<GEN->dim; i++)
      if (CLONE_CONDI[i])
	CLONE_CONDI[i]->distr = CLONE->distr_condi;
  }
  CLONE->direction = _unur_xmalloc( GEN->dim * sizeof(double));
  return clone;
#undef CLONE
#undef CLONE_CONDI
} 
void
_unur_gibbs_free( struct unur_gen *gen )
{
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_GIBBS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->state) free (GEN->state);
  if (GEN->x0) free (GEN->x0);
  if (GEN->direction) free (GEN->direction);
  if (GEN->distr_condi) _unur_distr_free (GEN->distr_condi);
  _unur_generic_free(gen);
} 
int
_unur_gibbs_coord_sample_cvec( struct unur_gen *gen, double *vec )
{
  double X;
  int thinning;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_GIBBS_GEN,UNUR_ERR_COOKIE);
  for (thinning = GEN->thinning; thinning > 0; --thinning) {
    GEN->coord = (GEN->coord + 1) % GEN->dim;
    if (!_unur_isfinite(GEN->state[GEN->coord]))
      continue;
    unur_distr_condi_set_condition( GEN->distr_condi, GEN->state, NULL, GEN->coord);
    if (unur_reinit(GEN_CONDI[GEN->coord]) == UNUR_SUCCESS) {
      X = unur_sample_cont(GEN_CONDI[GEN->coord]);
      if (_unur_isfinite(X)) {
	GEN->state[GEN->coord] = X;
	continue;
      }
    }
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"reset chain");
    unur_gibbs_reset_state(gen);
    return UNUR_FAILURE;
  }
  memcpy(vec, GEN->state, GEN->dim * sizeof(double)); 
  return UNUR_SUCCESS;
} 
int
_unur_gibbs_randomdir_sample_cvec( struct unur_gen *gen, double *vec )
{
  int i;
  double X;
  int thinning;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_GIBBS_GEN,UNUR_ERR_COOKIE);
  for (thinning = GEN->thinning; thinning > 0; --thinning) {
    if (!_unur_isfinite(GEN->state[0]))
      break;
    _unur_gibbs_random_unitvector( gen, GEN->direction );
    unur_distr_condi_set_condition( GEN->distr_condi, GEN->state, GEN->direction, 0);
    if (unur_reinit(*GEN_CONDI) == UNUR_SUCCESS) {
      X = unur_sample_cont(*GEN_CONDI);
      if (_unur_isfinite(X)) {
	for (i=0; i<GEN->dim; i++)
	  GEN->state[i] += X * GEN->direction[i];	  
	continue;
      }
    }
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"reset chain");
    unur_gibbs_reset_state(gen);
    return UNUR_FAILURE;
  }
  memcpy(vec, GEN->state, GEN->dim * sizeof(double)); 
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_gibbs_normalgen( struct unur_gen *gen )
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
_unur_gibbs_random_unitvector( struct unur_gen *gen, double *direction )
{
  int i;
  do {
    for (i=0; i<GEN->dim; i++) 
      direction[i] = unur_sample_cont(GEN_NORMAL);
    _unur_vector_normalize(GEN->dim, direction);
  } while (!_unur_isfinite(direction[0]));
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_gibbs_debug_init_start( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = GIBBS (Markov Chain - GIBBS sampler)\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:
    fprintf(LOG,"coordinate sampling (original Gibbs sampler)  [default]\n"); break;
  case GIBBS_VARIANT_RANDOMDIR:
    fprintf(LOG,"random directions\n"); break;
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: transformation T_c(x) for TDR method = ",gen->genid);
  switch( gen->variant & GIBBS_VARMASK_T ) {
  case GIBBS_VAR_T_LOG:
    fprintf(LOG,"log(x)  ... c = 0");                   break;
  case GIBBS_VAR_T_SQRT:
    fprintf(LOG,"-1/sqrt(x)  ... c = -1/2");            break;
  case GIBBS_VAR_T_POW:
    fprintf(LOG,"-x^(%g)  ... c = %g",GEN->c_T,GEN->c_T); break;
  }
  _unur_print_if_default(gen,GIBBS_SET_C);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_distr_cvec_debug( gen->distr, gen->genid );
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:
    fprintf(LOG,"%s: sampling routine = _unur_gibbs_coord_sample()\n",gen->genid);
    break;
  case GIBBS_VARIANT_RANDOMDIR:
    fprintf(LOG,"%s: sampling routine = _unur_gibbs_randomdir_sample()\n",gen->genid);
    break;
  }
  fprintf(LOG,"%s: thinning = %d",gen->genid,GEN->thinning);
  _unur_print_if_default(gen,GIBBS_SET_THINNING);
  fprintf(LOG,"\n%s: burn-in = %d",gen->genid,GEN->burnin);
  _unur_print_if_default(gen,GIBBS_SET_BURNIN);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_matrix_print_vector( GEN->dim, GEN->x0, "starting point = ", LOG, gen->genid, "\t   ");
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_gibbs_debug_init_finished( const struct unur_gen *gen, int success )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  if (success) 
    fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  else
    fprintf(LOG,"%s: INIT failed **********************\n",gen->genid);
} 
void
_unur_gibbs_debug_burnin_failed( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: Burn-in failed --> INIT failed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_gibbs_debug_init_condi( const struct unur_gen *gen )
{
  int i;
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:
    fprintf(LOG,"%s: generators for full conditional distributions = \n",gen->genid);
    fprintf(LOG,"%s:\t",gen->genid);
    for (i=0; i<GEN->dim; i++)
      fprintf(LOG,"[%s] ", GEN_CONDI[i]->genid);
    fprintf(LOG,"\n%s:\n",gen->genid);
    break;
  case GIBBS_VARIANT_RANDOMDIR:
    fprintf(LOG,"%s: generators for full conditional distributions = [%s]\n",gen->genid,
	    GEN_CONDI[0]->genid);
    fprintf(LOG,"%s: generator for random directions = [%s]\n",gen->genid,
	    GEN_NORMAL->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    break;
  }
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_gibbs_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_distr_cvec_info_domain(gen);
  _unur_string_append(info,"   center    = ");
  _unur_distr_info_vector( gen, unur_distr_cvec_get_center(gen->distr), GEN->dim);
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]");
    else
      _unur_string_append(info,"  [default]");
  }
  _unur_string_append(info,"\n\n");
  _unur_string_append(info,"method: GIBBS (GIBBS sampler [MCMC])\n");
  _unur_string_append(info,"   variant = %s\n",
		      ((gen->variant & GIBBS_VARMASK_VARIANT)==GIBBS_VARIANT_COORD)
		      ? "coordinate sampling [default]" : "random direction sampling");
  _unur_string_append(info,"   T_c(x) = ");
  switch( gen->variant & GIBBS_VARMASK_T ) {
  case GIBBS_VAR_T_LOG:
    _unur_string_append(info,"log(x)  ... c = 0\n"); break;
  case GIBBS_VAR_T_SQRT:
    _unur_string_append(info,"-1/sqrt(x)  ... c = -1/2\n"); break;
  case GIBBS_VAR_T_POW:
    _unur_string_append(info,"-x^(%g)  ... c = %g\n",GEN->c_T,GEN->c_T); break;
  }
  _unur_string_append(info,"   thinning = %d\n", GEN->thinning);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   rejection constant = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    switch (gen->variant & GIBBS_VARMASK_VARIANT) {
    case GIBBS_VARIANT_COORD:
      _unur_string_append(info,"   variant_coordinate  [default]\n"); break;
    case GIBBS_VARIANT_RANDOMDIR:
      _unur_string_append(info,"   variant_random_direction\n"); break;
    }
    _unur_string_append(info,"   c = %g  %s\n", GEN->c_T,
 			(gen->set & GIBBS_SET_C) ? "" : "[default]");
    _unur_string_append(info,"   thinning = %d  %s\n", GEN->thinning,
 			(gen->set & GIBBS_SET_THINNING) ? "" : "[default]");
    _unur_string_append(info,"   burnin = %d  %s\n", GEN->burnin,
 			(gen->set & GIBBS_SET_THINNING) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }
} 
#endif   
