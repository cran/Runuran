/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifdef DEBUG_STORE_IP 
#  undef DEBUG_STORE_IP
#endif
#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "ars.h"
#include "ars_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define ARS_VARFLAG_VERIFY    0x0100u   
#define ARS_VARFLAG_PEDANTIC  0x0800u   
#define ARS_DEBUG_REINIT    0x00000002u  
#define ARS_DEBUG_IV        0x00000010u
#define ARS_DEBUG_SPLIT     0x00010000u
#define ARS_SET_CPOINTS        0x001u
#define ARS_SET_N_CPOINTS      0x002u
#define ARS_SET_PERCENTILES    0x004u
#define ARS_SET_N_PERCENTILES  0x008u
#define ARS_SET_RETRY_NCPOINTS 0x010u
#define ARS_SET_MAX_IVS        0x020u
#define ARS_SET_MAX_ITER       0x040u   
#define GENTYPE "ARS"          
static struct unur_gen *_unur_ars_init( struct unur_par *par );
static int _unur_ars_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_ars_create( struct unur_par *par );
static struct unur_gen *_unur_ars_clone( const struct unur_gen *gen );
static void _unur_ars_free( struct unur_gen *gen);
static double _unur_ars_sample( struct unur_gen *generator );
static double _unur_ars_sample_check( struct unur_gen *generator );
static int _unur_ars_starting_cpoints( struct unur_gen *gen );
static int _unur_ars_starting_intervals( struct unur_gen *gen );
static int _unur_ars_interval_parameter( struct unur_gen *gen, struct unur_ars_interval *iv );
static struct unur_ars_interval *_unur_ars_interval_new( struct unur_gen *gen,
							 double x, double logfx );
static int _unur_ars_tangent_intersection_point( struct unur_gen *gen,
						 struct unur_ars_interval *iv, double *ipt );
static double _unur_ars_interval_logarea( struct unur_gen *gen, struct unur_ars_interval *iv,
					  double slope, double x );
static int _unur_ars_interval_split( struct unur_gen *gen,
				     struct unur_ars_interval *iv_old, double x, double logfx );
static int _unur_ars_improve_hat( struct unur_gen *gen, struct unur_ars_interval *iv,
				  double x, double logfx);
static int _unur_ars_make_area_table( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_ars_debug_init_start( const struct unur_gen *gen );
static void _unur_ars_debug_init_finished( const struct unur_gen *gen );
static void _unur_ars_debug_reinit_start( const struct unur_gen *gen );
static void _unur_ars_debug_reinit_retry( const struct unur_gen *gen );
static void _unur_ars_debug_reinit_finished( const struct unur_gen *gen );
static void _unur_ars_debug_free( const struct unur_gen *gen );
static void _unur_ars_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
static void _unur_ars_debug_split_start( const struct unur_gen *gen,
					 const struct unur_ars_interval *iv,
					 double x, double logfx );
static void _unur_ars_debug_split_stop( const struct unur_gen *gen,
					const struct unur_ars_interval *iv_left,
					const struct unur_ars_interval *iv_right );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_ars_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_ars_par*)par->datap) 
#define GEN       ((struct unur_ars_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define logPDF(x)  _unur_cont_logPDF((x),(gen->distr))   
#define dlogPDF(x) _unur_cont_dlogPDF((x),(gen->distr))  
#define scaled_logarea(iv)  ((iv)->logAhat - GEN->logAmax)
#define scaled_area(iv)     (exp(scaled_logarea(iv)))
#define rescaled_logf(logf) ((logf) - GEN->logAmax)
#define _unur_ars_getSAMPLE(gen) \
   ( ((gen)->variant & ARS_VARFLAG_VERIFY) \
     ? _unur_ars_sample_check : _unur_ars_sample )
struct unur_par *
unur_ars_new( const struct unur_distr* distr )
{
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.logpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"logPDF"); return NULL; }
  if (DISTR_IN.dlogpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of logPDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_ars_par) );
  COOKIE_SET(par,CK_ARS_PAR);
  par->distr              = distr;  
  PAR->starting_cpoints    = NULL;   
  PAR->n_starting_cpoints  = 2;      
  PAR->percentiles         = NULL;   
  PAR->n_percentiles       = 2;      
  PAR->retry_ncpoints      = 30;     
  PAR->max_ivs             = 200;    
  PAR->max_iter            = 10000;  
  par->method   = UNUR_METH_ARS;     
  par->variant  = 0u;                
  par->set      = 0u;               
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = par->urng;               
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_ars_init;
  return par;
} 
int
unur_ars_set_max_intervals( struct unur_par *par, int max_ivs )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ivs = max_ivs;
  par->set |= ARS_SET_MAX_IVS;
  return UNUR_SUCCESS;
} 
int 
unur_ars_set_cpoints( struct unur_par *par, int n_cpoints, const double *cpoints )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  if (n_cpoints < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 2. using defaults");
    n_cpoints = 2;
    cpoints = NULL;
  }
  if (cpoints)
    for( i=1; i<n_cpoints; i++ )
      if (cpoints[i] <= cpoints[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
  PAR->starting_cpoints = cpoints;
  PAR->n_starting_cpoints = n_cpoints;
  par->set |= ARS_SET_N_CPOINTS | ((cpoints) ? ARS_SET_CPOINTS : 0);
  return UNUR_SUCCESS;
} 
int
unur_ars_set_reinit_percentiles( struct unur_par *par, int n_percentiles, const double *percentiles )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }
  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
  if (percentiles) {
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }
  PAR->percentiles = percentiles;
  PAR->n_percentiles = n_percentiles;
  par->set |= ARS_SET_N_PERCENTILES | ((percentiles) ? ARS_SET_PERCENTILES : 0);
  return UNUR_SUCCESS;
} 
int
unur_ars_chg_reinit_percentiles( struct unur_gen *gen, int n_percentiles, const double *percentiles )
{
  int i;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }
  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
  if (percentiles) {
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }
  GEN->n_percentiles = n_percentiles;
  GEN->percentiles = _unur_xrealloc( GEN->percentiles, n_percentiles * sizeof(double) );
  if (percentiles) {
    memcpy( GEN->percentiles, percentiles, n_percentiles * sizeof(double) );
  }
  else {
    if (n_percentiles == 2) {
      GEN->percentiles[0] = 0.25;
      GEN->percentiles[1] = 0.75;
    }
    else {
      for (i=0; i<n_percentiles; i++ )
	GEN->percentiles[i] = (i + 1.) / (n_percentiles + 1.);
    }
  }
  gen->set |= ARS_SET_N_PERCENTILES | ((percentiles) ? ARS_SET_PERCENTILES : 0);
  return UNUR_SUCCESS;
} 
int
unur_ars_set_reinit_ncpoints( struct unur_par *par, int ncpoints )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }
  PAR->retry_ncpoints = ncpoints;
  par->set |= ARS_SET_RETRY_NCPOINTS; 
  return UNUR_SUCCESS;
} 
int
unur_ars_chg_reinit_ncpoints( struct unur_gen *gen, int ncpoints )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }
  GEN->retry_ncpoints = ncpoints;
  gen->set |= ARS_SET_RETRY_NCPOINTS; 
  return UNUR_SUCCESS;
} 
int
unur_ars_set_max_iter( struct unur_par *par, int max_iter )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of iterations");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_iter = max_iter;
  par->set |= ARS_SET_MAX_ITER;
  return UNUR_SUCCESS;
} 
int
unur_ars_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  par->variant = (verify) ? (par->variant | ARS_VARFLAG_VERIFY) : (par->variant & (~ARS_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_ars_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  gen->variant = (verify) 
    ? (gen->variant | ARS_VARFLAG_VERIFY) 
    : (gen->variant & (~ARS_VARFLAG_VERIFY));
  SAMPLE = _unur_ars_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
int
unur_ars_set_pedantic( struct unur_par *par, int pedantic )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );
  par->variant = (pedantic) ? (par->variant | ARS_VARFLAG_PEDANTIC) : (par->variant & (~ARS_VARFLAG_PEDANTIC));
  return UNUR_SUCCESS;
} 
double
unur_ars_get_loghatarea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ARS, INFINITY );
  return log(GEN->Atotal) + GEN->logAmax;
} 
struct unur_gen *
_unur_ars_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_ARS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_ARS_PAR,NULL);
  gen = _unur_ars_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_ars_debug_init_start(gen);
#endif
  if (_unur_ars_starting_cpoints(gen)!=UNUR_SUCCESS) {
    _unur_ars_free(gen); return NULL;
  }
  if (_unur_ars_starting_intervals(gen)!=UNUR_SUCCESS) {
    _unur_ars_free(gen); return NULL;
  }
  if (GEN->n_ivs > GEN->max_ivs) {
    GEN->max_ivs = GEN->n_ivs;
  }
  _unur_ars_make_area_table(gen);
  if (GEN->Atotal <= 0. || !_unur_isfinite(GEN->Atotal)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_ars_free(gen);
    return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_ars_debug_init_finished(gen);
#endif
  gen->status = UNUR_SUCCESS;
  return gen;
} 
int
_unur_ars_reinit( struct unur_gen *gen )
{
  struct unur_ars_interval *iv,*next;
  double *bak_cpoints;
  int bak_n_cpoints;
  int i;
  int n_trials;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );
  n_trials = 1;
  if (gen->set & ARS_SET_N_PERCENTILES) {
    if (GEN->starting_cpoints==NULL || (GEN->n_starting_cpoints != GEN->n_percentiles)) {
      GEN->n_starting_cpoints = GEN->n_percentiles;
      GEN->starting_cpoints = _unur_xrealloc( GEN->starting_cpoints, GEN->n_percentiles * sizeof(double));
    }
    for (i=0; i<GEN->n_percentiles; i++) {
      GEN->starting_cpoints[i] = unur_ars_eval_invcdfhat( gen, GEN->percentiles[i] );
      if (!_unur_isfinite(GEN->starting_cpoints[i])) 
	n_trials = 2;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & ARS_DEBUG_REINIT)
    _unur_ars_debug_reinit_start(gen);
#endif
  bak_n_cpoints = GEN->n_starting_cpoints;
  bak_cpoints = GEN->starting_cpoints;
  for (;; ++n_trials) {
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
    GEN->iv = NULL;
    GEN->n_ivs = 0;
    GEN->Atotal = 0.;
    GEN->logAmax = 0.;
    if (n_trials > 2) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points for reinit");
      GEN->n_starting_cpoints = bak_n_cpoints;
      GEN->starting_cpoints = bak_cpoints;
      return UNUR_FAILURE;
    }
    if (n_trials > 1) {
      GEN->n_starting_cpoints = GEN->retry_ncpoints;
      GEN->starting_cpoints = NULL;
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug & ARS_DEBUG_REINIT)
	_unur_ars_debug_reinit_retry(gen);
#endif
    }
    if (_unur_ars_starting_cpoints(gen)!=UNUR_SUCCESS)
      continue;
    if (_unur_ars_starting_intervals(gen)!=UNUR_SUCCESS)
      continue;
    if (GEN->n_ivs > GEN->max_ivs)
      GEN->max_ivs = GEN->n_ivs;
    _unur_ars_make_area_table(gen);
    if (GEN->Atotal <= 0.)
      continue;
    break;
  }
  if (n_trials > 1) {
    GEN->n_starting_cpoints = bak_n_cpoints;
    GEN->starting_cpoints = bak_cpoints;
  }
  SAMPLE = _unur_ars_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & ARS_DEBUG_REINIT)
    _unur_ars_debug_reinit_finished(gen);
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_ars_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_ARS_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_ars_gen) );
  COOKIE_SET(gen,CK_ARS_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_ars_getSAMPLE(gen);
  gen->destroy = _unur_ars_free;
  gen->clone = _unur_ars_clone;
  gen->reinit = _unur_ars_reinit;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  GEN->percentiles = NULL;
  GEN->Atotal      = 0.;
  GEN->logAmax     = 0.;
  GEN->n_starting_cpoints = PAR->n_starting_cpoints;
  if (PAR->starting_cpoints) {
    GEN->starting_cpoints = _unur_xmalloc( PAR->n_starting_cpoints * sizeof(double) );
    memcpy( GEN->starting_cpoints, PAR->starting_cpoints, PAR->n_starting_cpoints * sizeof(double) );
  }
  else {
    GEN->starting_cpoints = NULL;
  }
  if (gen->set & ARS_SET_N_PERCENTILES)
    unur_ars_chg_reinit_percentiles( gen, PAR->n_percentiles, PAR->percentiles );
  GEN->retry_ncpoints = PAR->retry_ncpoints;   
  GEN->max_ivs = _unur_max(2*PAR->n_starting_cpoints,PAR->max_ivs);  
  GEN->max_iter = PAR->max_iter;
  gen->variant = par->variant;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_ars_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_ars_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_ars_gen*)clone->datap)
  struct unur_gen *clone;
  struct unur_ars_interval *iv, *clone_iv, *clone_prev;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = iv->next) {
    clone_iv = _unur_xmalloc( sizeof(struct unur_ars_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_ars_interval) );
    if (clone_prev == NULL) {
      CLONE->iv = clone_iv;
    }
    else {
      clone_prev->next = clone_iv;
    }
    clone_prev = clone_iv;
  }
  if (clone_iv) clone_iv->next = NULL;
  if (GEN->starting_cpoints) {
    CLONE->starting_cpoints = _unur_xmalloc( GEN->n_starting_cpoints * sizeof(double) );
    memcpy( CLONE->starting_cpoints, GEN->starting_cpoints, GEN->n_starting_cpoints * sizeof(double) );
  }
  if (GEN->percentiles) {
    CLONE->percentiles = _unur_xmalloc( GEN->n_percentiles * sizeof(double) );
    memcpy( CLONE->percentiles, GEN->percentiles, GEN->n_percentiles * sizeof(double) );
  }
  return clone;
#undef CLONE
} 
void
_unur_ars_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_ARS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  SAMPLE = NULL;   
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_ars_debug_free(gen);
#endif
  {
    struct unur_ars_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }
  if (GEN->starting_cpoints) 
    free (GEN->starting_cpoints);
  if (GEN->percentiles) 
    free (GEN->percentiles);
  _unur_generic_free(gen);
} 
double
_unur_ars_sample( struct unur_gen *gen )
{ 
  struct unur_ars_interval *iv, *cp;
  double U, logV;                   
  double X;                         
  double logfx, logsqx, loghx;      
  double x0, logfx0, dlogfx0, fx0;  
  int n_trials;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ARS_GEN,INFINITY);
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 
  for (n_trials=0; n_trials<GEN->max_iter; ++n_trials) {
    U = _unur_call_urng(gen->urng);
    iv =  GEN->iv;
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }
    U -= iv->Acum;    
    if (-U < (scaled_area(iv) * iv->Ahatr_fract)) { 
      cp = iv->next;
    }
    else {                
      cp = iv;
      U += scaled_area(iv);
    }
    x0 = cp->x;
    logfx0 = cp->logfx;
    dlogfx0 = cp->dlogfx;
    fx0 = exp(rescaled_logf(logfx0));
    if (_unur_iszero(dlogfx0))
      X = x0 + U / fx0;
    else {
      double t = dlogfx0 * U / fx0;
      if (fabs(t) > 1.e-6)
	X = x0 + log(t + 1.) * U / (fx0 * t);
      else if (fabs(t) > 1.e-8)
	X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
      else
	X = x0 + U / fx0 * (1 - t/2.);
    }
    loghx = rescaled_logf(logfx0) + dlogfx0*(X - x0);
    logV = log(_unur_call_urng(gen->urng)) + loghx;
    logsqx = rescaled_logf(iv->logfx) + iv->sq*(X - iv->x);
    if (logV <= logsqx)
      return X;
    logfx = logPDF(X);
    if (logV <= rescaled_logf(logfx))
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if (! (_unur_isfinite(X) && _unur_isfinite(logfx)) ) {
	X = _unur_arcmean(iv->x,iv->next->x);  
	logfx = logPDF(X);
      }
      if ( (_unur_ars_improve_hat( gen, iv, X, logfx) != UNUR_SUCCESS)
	   && (gen->variant & ARS_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
  }
  _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"max number of iterations exceeded");
  return UNUR_INFINITY;
} 
double
_unur_ars_sample_check( struct unur_gen *gen )
{ 
  struct unur_ars_interval *iv, *cp;
  double U, logV;                   
  double X;                         
  double logfx, logsqx, loghx;      
  double x0, logfx0, dlogfx0, fx0;  
  int n_trials;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ARS_GEN,INFINITY);
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"empty generator object");
    return INFINITY;
  } 
  for (n_trials=0; n_trials<GEN->max_iter; ++n_trials) {
    U = _unur_call_urng(gen->urng);
    iv =  GEN->iv;
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }
    U -= iv->Acum;    
    if (-U < (scaled_area(iv) * iv->Ahatr_fract)) { 
      cp = iv->next;
    }
    else {                
      cp = iv;
      U += scaled_area(iv);
    }
    x0 = cp->x;
    logfx0 = cp->logfx;
    dlogfx0 = cp->dlogfx;
    fx0 = exp(rescaled_logf(logfx0));
    if (_unur_iszero(dlogfx0))
      X = x0 + U / fx0;
    else {
      double t = dlogfx0 * U / fx0;
      if (fabs(t) > 1.e-6)
	X = x0 + log(t + 1.) * U / (fx0 * t);
      else if (fabs(t) > 1.e-8)
	X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
      else
	X = x0 + U / fx0 * (1 - t/2.);
    }
    loghx = rescaled_logf(logfx0) + dlogfx0*(X - x0);
    logsqx = rescaled_logf(iv->logfx) + iv->sq*(X - iv->x);
    logfx = logPDF(X);
    if (X < DISTR.BD_LEFT || X > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
    }
    if (_unur_FP_greater(rescaled_logf(logfx), loghx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not log-concave!");
    }
    if (_unur_FP_less(rescaled_logf(logfx), logsqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not log-concave!");
    }
    logV = log(_unur_call_urng(gen->urng)) + loghx;
    if (logV <= logsqx)
      return X;
    if (logV <= rescaled_logf(logfx))
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if (! (_unur_isfinite(X) && _unur_isfinite(logfx)) ) {
	X = _unur_arcmean(iv->x,iv->next->x);  
	logfx = logPDF(X);
      }
      if ( (_unur_ars_improve_hat( gen, iv, X, logfx) != UNUR_SUCCESS)
	   && (gen->variant & ARS_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
  }
  _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"max number of iterations exceeded");
  return UNUR_INFINITY;
} 
double
unur_ars_eval_invcdfhat( const struct unur_gen *gen, double U )
{ 
  struct unur_ars_interval *iv, *cp;
  double X;                         
  double x0, logfx0, dlogfx0, fx0;  
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_ARS ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_ARS_GEN,INFINITY);
  if ( U<0. || U>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 
  if (U<=0.) return DISTR.domain[0];
  if (U>=1.) return DISTR.domain[1];
  iv =  GEN->iv;
  U *= GEN->Atotal;
  while (iv->Acum < U) {
    iv = iv->next;
  }
  U -= iv->Acum;
  if (-U < (scaled_area(iv) * iv->Ahatr_fract)) { 
    cp = iv->next;
  }
  else {                
    cp = iv;
    U += scaled_area(iv);
  }
  x0 = cp->x;
  logfx0 = cp->logfx;
  dlogfx0 = cp->dlogfx;
  fx0 = exp(rescaled_logf(logfx0));
  if (_unur_iszero(dlogfx0))
    X = x0 + U / fx0;
  else {
    double t = dlogfx0 * U / fx0;
    if (fabs(t) > 1.e-6)
      X = x0 + log(t + 1.) * U / (fx0 * t);
    else if (fabs(t) > 1.e-8)
      X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
    else
      X = x0 + U / fx0 * (1 - t/2.);
  }
  return X;
} 
int
_unur_ars_improve_hat( struct unur_gen *gen, struct unur_ars_interval *iv,
			  double x, double logfx )
{
  int result;
  result = _unur_ars_interval_split(gen, iv, x, logfx);
  if (result!=UNUR_SUCCESS && result!=UNUR_ERR_SILENT) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
    if (gen->variant & ARS_VARFLAG_PEDANTIC) {
      SAMPLE = _unur_sample_cont_error;
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  _unur_ars_make_area_table(gen);
  return UNUR_SUCCESS;
} 
int
_unur_ars_starting_cpoints( struct unur_gen *gen )
{
  struct unur_ars_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, logfx, logfx_last;
  int is_increasing;
  int i;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  GEN->n_ivs = 0;
  if (!GEN->starting_cpoints) {
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT);
    diff_angle = (right_angle-left_angle) / (GEN->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   
  x = DISTR.BD_LEFT;
  is_increasing = TRUE;
  logfx = logfx_last = _unur_isfinite(x) ? logPDF(x) : -INFINITY;
  iv = GEN->iv = _unur_ars_interval_new( gen, x, logfx );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  
  for( i=0; i<=GEN->n_starting_cpoints; i++ ) {
    if (i < GEN->n_starting_cpoints) {
      if (GEN->starting_cpoints) {
	x = GEN->starting_cpoints[i];
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	angle += diff_angle;
	x = tan( angle );
      }
    }
    else {
      x = DISTR.BD_RIGHT;
    }
    logfx = _unur_isfinite(x) ? logPDF(x) : -INFINITY;
    if (!is_increasing && logfx > logfx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal!");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (!_unur_isfinite(logfx) && !_unur_isfinite(logfx_last) ) {
      if (is_increasing) {
	if (i<GEN->n_starting_cpoints) {
	  iv->x = x;  
	  continue;   
	}
      }
      else
	break;
    }
    iv->next = _unur_ars_interval_new( gen, x, logfx );
    if (iv->next == NULL) return UNUR_ERR_GEN_DATA;    
    iv = iv->next;
    if (is_increasing && logfx < logfx_last)
      is_increasing = FALSE;
    logfx_last = logfx;
  }
  iv->logAhat = -INFINITY;
  iv->Ahatr_fract = iv->sq = 0.;
  iv->Acum = INFINITY;
#ifdef DEBUG_STORE_IP 
  iv->ip = iv->x;
#endif
  iv->next = NULL;         
  --(GEN->n_ivs);          
  return UNUR_SUCCESS;
} 
int
_unur_ars_starting_intervals( struct unur_gen *gen )
{
  struct unur_ars_interval *iv, *iv_new, *iv_tmp;
  double x, logfx;              
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_ARS_IV,UNUR_ERR_COOKIE);
  for( iv=GEN->iv; iv->next != NULL; ) {
    switch (_unur_ars_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:      
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:      
      break;
    case UNUR_ERR_SILENT:   
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN->n_ivs);
      if (iv->next==NULL) {
	iv->logAhat = -INFINITY;
	iv->Ahatr_fract = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      continue;
    default:     
      return UNUR_ERR_GEN_CONDITION;
    }
    x = _unur_arcmean(iv->x,iv->next->x);  
    logfx = logPDF(x);
    if (GEN->n_ivs >= GEN->max_ivs) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = _unur_ars_interval_new( gen, x, logfx );
    if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  
    if (!_unur_isfinite(logfx) ) {
      if (!_unur_isfinite(iv->logfx) ) {
	iv_new->next = iv->next;
	free(iv);
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	iv = iv_new;
      }
      else if (!_unur_isfinite(iv->next->logfx) ) {
	free(iv->next);
	--(GEN->n_ivs);
	iv->next = iv_new;
      }
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	free(iv_new);
	return UNUR_ERR_GEN_CONDITION;
      }
    }
    else {
      iv_new->next = iv->next;
      iv->next = iv_new;
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_ars_interval *
_unur_ars_interval_new( struct unur_gen *gen, double x, double logfx )
{
  struct unur_ars_interval *iv;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,NULL);
  if (!(logfx < INFINITY)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"logPDF(x) overflow");
    return NULL;
  }
  iv = _unur_xmalloc( sizeof(struct unur_ars_interval) );
  iv->next = NULL; 
  ++(GEN->n_ivs);   
  COOKIE_SET(iv,CK_ARS_IV);
  iv->logAhat = -INFINITY;
  iv->Acum = iv->Ahatr_fract = 0.;
  iv->sq = 0.;
#ifdef DEBUG_STORE_IP 
  iv->ip = 0.;
#endif
  iv->x = x;              
  iv->logfx = logfx;      
  iv->dlogfx = _unur_isfinite(logfx) ? dlogPDF(x) : INFINITY;
  if ( !(iv->dlogfx > -INFINITY))
    iv->dlogfx = INFINITY;
  return iv;
} 
int
_unur_ars_interval_parameter( struct unur_gen *gen, struct unur_ars_interval *iv )
{
  double logAhatl, logAhatr;  
  double ip = 0.;             
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE); 
  CHECK_NULL(iv->next,UNUR_ERR_NULL);  COOKIE_CHECK(iv->next,CK_ARS_IV,UNUR_ERR_COOKIE); 
  if ( _unur_ars_tangent_intersection_point(gen,iv,&ip)!=UNUR_SUCCESS )
    return UNUR_ERR_GEN_CONDITION;
#ifdef DEBUG_STORE_IP 
  iv->ip = ip;
#endif
  if (_unur_isfinite(iv->logfx) && _unur_isfinite(iv->next->dlogfx) ) {
    if (_unur_FP_approx(iv->x, iv->next->x) )
      return UNUR_ERR_SILENT;   
    iv->sq = (iv->next->logfx - iv->logfx) / (iv->next->x - iv->x);
    if ( ( (iv->sq > iv->dlogfx      && (!_unur_FP_approx(iv->sq,iv->dlogfx)) ) ||
	   (iv->sq < iv->next->dlogfx && (!_unur_FP_approx(iv->sq,iv->next->dlogfx)) ) )
	 && iv->next->dlogfx < INFINITY ) {   
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. PDF not T-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  else {  
    iv->sq = -INFINITY;
  }
  logAhatl = _unur_ars_interval_logarea( gen, iv, iv->dlogfx, ip);
  logAhatr = _unur_ars_interval_logarea( gen, iv->next, iv->next->dlogfx, ip);
  if (! (logAhatl < INFINITY && logAhatr < INFINITY) )
    return UNUR_ERR_INF;
  iv->logAhat = (logAhatl > logAhatr) 
    ? logAhatl+log(1+exp(logAhatr-logAhatl)) 
    : logAhatr+log(1+exp(logAhatl-logAhatr)) ;       
  iv->Ahatr_fract = 1./(1.+exp(logAhatl-logAhatr));  
  return UNUR_SUCCESS;
} 
int
_unur_ars_interval_split( struct unur_gen *gen, struct unur_ars_interval *iv_oldl, double x, double logfx )
{
  struct unur_ars_interval *iv_newr;  
  struct unur_ars_interval iv_bak;    
  int success, success_r;
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_oldl,UNUR_ERR_NULL);  COOKIE_CHECK(iv_oldl,CK_ARS_IV,UNUR_ERR_COOKIE);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & ARS_DEBUG_SPLIT)
    _unur_ars_debug_split_start( gen,iv_oldl,x,logfx );
#endif
  if (x < iv_oldl->x || x > iv_oldl->next->x) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }
  memcpy(&iv_bak, iv_oldl, sizeof(struct unur_ars_interval));
  if (!_unur_isfinite(logfx)) {
    if (!_unur_isfinite(iv_oldl->logfx)) {
      iv_oldl->x = x;
    }
    else if (!_unur_isfinite(iv_oldl->next->logfx)) {
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not log-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    success = _unur_ars_interval_parameter(gen, iv_oldl);
    iv_newr = NULL;
  }
  else {
    iv_newr = _unur_ars_interval_new( gen, x, logfx );
    if (iv_newr == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    iv_newr->next = iv_oldl->next;
    iv_oldl->next = iv_newr;
    success   = _unur_ars_interval_parameter(gen, iv_oldl);
    success_r = _unur_ars_interval_parameter(gen, iv_newr);
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  if (success!=UNUR_SUCCESS) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not log-concave");
    memcpy(iv_oldl, &iv_bak, sizeof(struct unur_ars_interval));
    if (iv_newr) {
      --(GEN->n_ivs);
      free( iv_newr );
    }
  return ( (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
	   ? UNUR_ERR_GEN_CONDITION : UNUR_SUCCESS );
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & ARS_DEBUG_SPLIT) {
    GEN->Atotal = ( GEN->Atotal - scaled_area(&iv_bak) + scaled_area(iv_oldl) + ((iv_newr) ? scaled_area(iv_newr) : 0.) );
    _unur_ars_debug_split_stop( gen,iv_oldl,iv_newr );
  }
#endif
  return UNUR_SUCCESS;
} 
int
_unur_ars_tangent_intersection_point( struct unur_gen *gen, struct unur_ars_interval *iv, double *ipt )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE); 
  if ( iv->dlogfx > 1.e+140 ) {
    *ipt = iv->x;        
    return UNUR_SUCCESS;
  }
  if ( iv->next->dlogfx < -1.e+140 || _unur_FP_is_infinity(iv->next->dlogfx)) {
    *ipt = iv->next->x;   
    return UNUR_SUCCESS;
  }
  if ( _unur_FP_less( iv->dlogfx, iv->next->dlogfx ) ) {
    if ( fabs(iv->dlogfx) < DBL_EPSILON * fabs(iv->next->dlogfx) ) {
      *ipt = iv->x;        
      iv->dlogfx = INFINITY;
      return UNUR_SUCCESS;
    }
    else if ( fabs(iv->next->dlogfx) < DBL_EPSILON * fabs(iv->dlogfx) ) {
      *ipt = iv->next->x;   
      iv->next->dlogfx = INFINITY;
      return UNUR_SUCCESS;
    }
    else {
      if (_unur_FP_approx(iv->dlogfx, iv->next->dlogfx)) {
        *ipt = 0.5 * (iv->x + iv->next->x);
        return UNUR_SUCCESS;
      }
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). PDF not log-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  if (_unur_FP_approx(iv->dlogfx, iv->next->dlogfx)) {
    *ipt = 0.5 * (iv->x + iv->next->x);
    return UNUR_SUCCESS;
  }
  *ipt = ( (iv->next->logfx - iv->logfx - iv->next->dlogfx * iv->next->x + iv->dlogfx * iv->x) /
	   (iv->dlogfx - iv->next->dlogfx) );
  if (_unur_FP_less(*ipt, iv->x) || _unur_FP_greater(*ipt, iv->next->x))
    *ipt = 0.5 * (iv->x + iv->next->x);
  return UNUR_SUCCESS;
} 
double
_unur_ars_interval_logarea( struct unur_gen *gen ATTRIBUTE__UNUSED, 
			      struct unur_ars_interval *iv, double slope, double x )
{
  double x0, logfx0;
  double logxdiff;
  double t, logt;
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_ARS_IV,INFINITY); 
  if (_unur_FP_same(x, iv->x))
    return -INFINITY;
  if (!_unur_isfinite(iv->x)) 
    return INFINITY;
  if ( !_unur_isfinite(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   
    return INFINITY;
  x0 = iv->x;
  logfx0 = iv->logfx;
  logxdiff = log(fabs(x - x0));
  if (_unur_iszero(slope))
    return (_unur_isfinite(x) ? logfx0 + logxdiff : INFINITY);
  if (!_unur_isfinite(x))
    return (logfx0 - log(fabs(slope)));
  t = slope * (x - x0);
  logt = log(fabs(slope)) + logxdiff;
  if (fabs(t) > 1.e-6) {
    if (t > MAXLOG / 10.)
      return ( logfx0 + logxdiff + t - logt );
    else
      return ( logfx0 + logxdiff + log( fabs(exp(t) - 1.) ) - log(fabs(t)) );
  }
  else  
    return (logfx0 + logxdiff + log1p(t/2. + t*t/6.));
} 
int
_unur_ars_make_area_table( struct unur_gen *gen )
{
  struct unur_ars_interval *iv;
  double Acum;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  GEN->logAmax = -INFINITY;
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE);
    if (GEN->logAmax < iv->logAhat)
      GEN->logAmax = iv->logAhat;
  }
  Acum = 0.;            
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE);
    Acum += scaled_area(iv);
    iv->Acum = Acum;
  }
  GEN->Atotal = Acum;
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_ars_debug_init_start( const struct unur_gen *gen )
{
  FILE *LOG;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ARS (Adaptive Rejection Sampling)\n",gen->genid);
  fprintf(LOG,"%s: transformation T_c(x) = log(x)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->distr_is_privatecopy)
    fprintf(LOG,"%s: use private copy of distribution object\n",gen->genid);
  else
    fprintf(LOG,"%s: use pointer to external distribution object (dangerous!)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_ars_sample",gen->genid);
  if (gen->variant & ARS_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: maximum number of intervals        = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(gen,ARS_SET_MAX_IVS);
  fprintf(LOG,"\n");
  fprintf(LOG,"%s: maximum number of iterations       = %d",gen->genid,GEN->max_iter);
  _unur_print_if_default(gen,ARS_SET_MAX_ITER);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: number of starting points = %d",gen->genid,GEN->n_starting_cpoints);
  _unur_print_if_default(gen,ARS_SET_N_CPOINTS);
  fprintf(LOG,"\n%s: starting points:",gen->genid);
  if (gen->set & ARS_SET_CPOINTS)
    for (i=0; i<GEN->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",GEN->starting_cpoints[i]);
    }
  else
    fprintf(LOG," use \"equidistribution\" rule [default]");
  fprintf(LOG,"\n%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_ars_debug_init_finished( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  _unur_ars_debug_intervals(gen,"INIT completed",TRUE);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_ars_debug_reinit_start( const struct unur_gen *gen )
{
  int i;
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: *** Re-Initialize generator object ***\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->set & ARS_SET_N_PERCENTILES) {
    fprintf(LOG,"%s: use percentiles of old hat as starting points for new hat:",gen->genid);
    for (i=0; i<GEN->n_percentiles; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",GEN->percentiles[i]);
    }
    fprintf(LOG,"\n%s: starting points:",gen->genid);
    for (i=0; i<GEN->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",GEN->starting_cpoints[i]);
    }
    fprintf(LOG,"\n");
  }
  else {
    fprintf(LOG,"%s: use starting points given at init\n",gen->genid);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_ars_debug_reinit_retry( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: *** Re-Initialize failed  -->  second trial ***\n",gen->genid);
  fprintf(LOG,"%s: use equal-area-rule with %d points\n",gen->genid,GEN->retry_ncpoints);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_ars_debug_reinit_finished( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  _unur_ars_debug_intervals(gen," *** Generator reinitialized ***",TRUE);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void 
_unur_ars_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas )
{
  FILE *LOG;
  struct unur_ars_interval *iv;
  double Ahat, Ahatl, Ahatr;
  double sAhatl, sAhatr, Atotal, logAmax;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  if (header) fprintf(LOG,"%s:%s\n",gen->genid,header);
  fprintf(LOG,"%s:Intervals: %d\n",gen->genid,GEN->n_ivs);
  if (GEN->iv) {
    if (gen->debug & ARS_DEBUG_IV) {
#ifdef DEBUG_STORE_IP 
      fprintf(LOG,"%s: Nr.            tp            ip       logf(tp)     dlogf(tp)       squeeze\n",gen->genid);
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
	fprintf(LOG,"%s:[%3d]: %#12.6g  %#12.6g  %#12.6g  %#12.6g  %#12.6g\n", gen->genid, i,
		iv->x, iv->ip, iv->logfx, iv->dlogfx, iv->sq);
      }
#else
      fprintf(LOG,"%s: Nr.            tp       logf(tp)     dlogf(tp)       squeeze\n",gen->genid);
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
	fprintf(LOG,"%s:[%3d]: %#12.6g  %#12.6g  %#12.6g  %#12.6g\n", gen->genid, i,
		iv->x, iv->logfx, iv->dlogfx, iv->sq);
      }
#endif
      COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
      fprintf(LOG,"%s:[...]: %#12.6g                %#12.6g  %#12.6g\n", gen->genid,
	      iv->x, iv->logfx, iv->dlogfx);
    }
    fprintf(LOG,"%s:\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: No intervals !\n",gen->genid);
  if (!print_areas || GEN->Atotal <= 0.) return;
  Atotal = GEN->Atotal;
  logAmax = GEN->logAmax;
  if (gen->debug & ARS_DEBUG_IV) {
    fprintf(LOG,"%s:Areas in intervals relative to maximum:\t[ log(A_max) = %g ]\n",gen->genid, logAmax);
    fprintf(LOG,"%s: Nr.\tbelow hat (left and right)\t\t   cumulated\n",gen->genid);
    sAhatl = sAhatr = 0.;
    if (GEN->iv) {
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
	Ahat = scaled_area(iv);
	sAhatr += Ahatr = Ahat * iv->Ahatr_fract;
	sAhatl += Ahatl = Ahat - Ahatr;
	fprintf(LOG,"%s:[%3d]: %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		Ahatl, Ahatr, Ahat * 100. / Atotal,
		iv->Acum, iv->Acum * 100. / Atotal);
      }
      fprintf(LOG,"%s:       ------------------------  ---------  +\n",gen->genid);
      fprintf(LOG,"%s: Sum :        %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / Atotal);
      fprintf(LOG,"%s:\n",gen->genid);
    }
  }
  fprintf(LOG,"%s: A(total) = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_ars_debug_free( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->status == UNUR_SUCCESS) {
    fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    _unur_ars_debug_intervals(gen,NULL,TRUE);
  }
  else {
    fprintf(LOG,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
    _unur_ars_debug_intervals(gen,"Intervals after failure:",FALSE);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_ars_debug_split_start( const struct unur_gen *gen, 
			     const struct unur_ars_interval *iv,
			     double x, double logfx )
{
  FILE *LOG;
  double Ahat;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  CHECK_NULL(iv,RETURN_VOID);   COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID);
  LOG = unur_get_stream();
  Ahat = scaled_area(iv);
  fprintf(LOG,"%s: split interval at x = %g \t\tlogf(x) = %g\n",gen->genid,x,logfx);
  fprintf(LOG,"%s: old interval:\n",gen->genid);
  fprintf(LOG,"%s:   left  construction point = %-12.6g\tlogf(x) = %-12.6g\n",gen->genid,iv->x,iv->logfx);
  fprintf(LOG,"%s:   right construction point = %-12.6g\tlogf(x) = %-12.6g\n",gen->genid,iv->next->x,iv->next->logfx);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\t[ relative to A_max ]\n",gen->genid,
	  Ahat * (1.-iv->Ahatr_fract), Ahat * iv->Ahatr_fract, Ahat*100./GEN->Atotal);
  fflush(LOG);
} 
void
_unur_ars_debug_split_stop( const struct unur_gen *gen, 
			    const struct unur_ars_interval *iv_left, 
			    const struct unur_ars_interval *iv_right )
{
  FILE *LOG;
  double Ahat;
  CHECK_NULL(gen,RETURN_VOID);       COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  CHECK_NULL(iv_left,RETURN_VOID);   COOKIE_CHECK(iv_left,CK_ARS_IV,RETURN_VOID);
  if (iv_right == NULL) iv_right = iv_left;
  LOG = unur_get_stream();
  fprintf(LOG,"%s: inserted point:\n",gen->genid);
  fprintf(LOG,"%s: x = %g, logf(x) = %g, dlogf(x) = %g, squeeze = %g:\n",
	  gen->genid, iv_right->x, iv_right->logfx, iv_right->dlogfx, iv_right->sq);
  fprintf(LOG,"%s: new intervals:\n",gen->genid);
  fprintf(LOG,"%s:   left   construction point = %g\n",gen->genid, iv_left->x);
  if (iv_left != iv_right)
    fprintf(LOG,"%s:   middle construction point = %g\n",gen->genid, iv_right->x);
  fprintf(LOG,"%s:   right  construction point = %g\n",gen->genid, iv_right->next->x);
  fprintf(LOG,"%s: left interval:\n",gen->genid);
  Ahat = scaled_area(iv_left);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\t[ relative to A_max ]\n",gen->genid,
	  Ahat * (1.-iv_left->Ahatr_fract), Ahat * iv_left->Ahatr_fract,
	  Ahat * 100./GEN->Atotal);
  if (iv_left == iv_right)
    fprintf(LOG,"%s: interval chopped.\n",gen->genid);
  else {
    fprintf(LOG,"%s: right interval:\n",gen->genid);
    Ahat = scaled_area(iv_right);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\t[ relative to A_max ]\n",gen->genid,
	    Ahat * (1.-iv_right->Ahatr_fract), Ahat * iv_right->Ahatr_fract,
	    Ahat * 100./GEN->Atotal);
  }
  fprintf(LOG,"%s: total areas:\n",gen->genid);
  fprintf(LOG,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
#endif    
#ifdef UNUR_ENABLE_INFO
void
_unur_ars_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  int n_ivs_bak;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = logPDF dlogPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: ARS (Transformed Density Rejection -- Gilks&Wild variant)\n");
  _unur_string_append(info,"   T_c(x) = log(x)  ... c = 0\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g  [ log = %g ]\n", 
		      GEN->Atotal*exp(GEN->logAmax), log(GEN->Atotal) + GEN->logAmax);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", GEN->Atotal*exp(GEN->logAmax)/DISTR.area);
  else {
    n_ivs_bak = GEN->n_ivs;
    GEN->n_ivs = GEN->max_ivs+1;
    _unur_string_append(info,"= %.3f  [approx.]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
    GEN->n_ivs = n_ivs_bak;
  }
  _unur_string_append(info,"   # intervals = %d\n", GEN->n_ivs);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   cpoints = %d  %s\n", GEN->n_starting_cpoints,
			(gen->set & ARS_SET_N_CPOINTS) ? "" : "[default]");
    if (gen->variant & ARS_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    if (gen->variant & ARS_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");
    _unur_string_append(info,"\n");
  }
} 
#endif
