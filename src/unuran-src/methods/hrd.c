/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "hrd.h"
#include "hrd_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define HRD_VARFLAG_VERIFY     0x01u    
#define HRD_DEBUG_REINIT    0x00000010u   
#define HRD_DEBUG_SAMPLE       0x01000000u    
#define GENTYPE "HRD"          
static struct unur_gen *_unur_hrd_init( struct unur_par *par );
static int _unur_hrd_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_hrd_create( struct unur_par *par );
static int _unur_hrd_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_hrd_clone( const struct unur_gen *gen );
static void _unur_hrd_free( struct unur_gen *gen );
static double _unur_hrd_sample( struct unur_gen *gen );
static double _unur_hrd_sample_check( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_hrd_debug_init( const struct unur_gen *gen );
static void _unur_hrd_debug_sample( const struct unur_gen *gen, double x, int i );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_hrd_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_hrd_par*)par->datap) 
#define GEN       ((struct unur_hrd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define SAMPLE    gen->sample.cont      
#define HR(x)     _unur_cont_HR((x),(gen->distr))   
#define _unur_hrd_getSAMPLE(gen) \
   ( ((gen)->variant & HRD_VARFLAG_VERIFY) \
     ? _unur_hrd_sample_check : _unur_hrd_sample )
struct unur_par *
unur_hrd_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.hr == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"HR"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_hrd_par) );
  COOKIE_SET(par,CK_HRD_PAR);
  par->distr   = distr;         
  par->method   = UNUR_METH_HRD;  
  par->variant  = 0u;             
  par->set      = 0u;                      
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_hrd_init;
  return par;
} 
int
unur_hrd_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HRD );
  par->variant = (verify) ? (par->variant | HRD_VARFLAG_VERIFY) : (par->variant & (~HRD_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_hrd_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, HRD, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  gen->variant = (verify) 
    ? (gen->variant | HRD_VARFLAG_VERIFY) 
    : (gen->variant & (~HRD_VARFLAG_VERIFY));
  SAMPLE = _unur_hrd_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hrd_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_HRD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HRD_PAR,NULL);
  gen = _unur_hrd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_hrd_check_par(gen) != UNUR_SUCCESS) {
    _unur_hrd_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_hrd_debug_init(gen);
#endif
  return gen;
} 
int
_unur_hrd_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_hrd_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  SAMPLE = _unur_hrd_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hrd_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HRD_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_hrd_gen) );
  COOKIE_SET(gen,CK_HRD_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_hrd_getSAMPLE(gen);
  gen->destroy = _unur_hrd_free;
  gen->clone = _unur_hrd_clone;
  gen->reinit = _unur_hrd_reinit;
  GEN->left_border = 0.;             
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_hrd_info;
#endif
  return gen;
} 
int
_unur_hrd_check_par( struct unur_gen *gen )
{
  if (DISTR.domain[0] < 0.)            DISTR.domain[0] = 0.;
  if (DISTR.domain[1] < UNUR_INFINITY) DISTR.domain[1] = UNUR_INFINITY;
  GEN->left_border = DISTR.domain[0];
  GEN->upper_bound = HR(GEN->left_border);
  if (GEN->upper_bound <= 0. || _unur_FP_is_infinity(GEN->upper_bound)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"no valid upper bound for HR at left boundary");
    return UNUR_ERR_GEN_CONDITION;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hrd_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_hrd_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HRD_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_hrd_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_HRD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HRD_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_hrd_sample( struct unur_gen *gen )
{ 
  double U,V,E,X,hrx;
  double lambda;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_HRD_GEN,UNUR_INFINITY);
  lambda = GEN->upper_bound;
  X = GEN->left_border;
  for(;;) {
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );
    E = -log(U) / lambda;
    X += E;
    hrx = HR(X);
    V =  lambda * _unur_call_urng(gen->urng);
    if( V <= hrx ) 
      return X;      
    if (hrx > 0.)
      lambda = hrx;   
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not valid");
      return UNUR_INFINITY;
    } 
  }
} 
double
_unur_hrd_sample_check( struct unur_gen *gen )
{ 
  double U,V,E,X,hrx;
  double lambda;
  int i;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_HRD_GEN,UNUR_INFINITY);
  lambda = GEN->upper_bound;
  X = GEN->left_border;
  for(i=1;;i++) {
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );
    E = -log(U) / lambda;
    X += E;
    hrx = HR(X);
    if ( (1.+UNUR_EPSILON) * lambda < hrx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not decreasing");
    V =  lambda * _unur_call_urng(gen->urng);
    if( V <= hrx ) {
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug & HRD_DEBUG_SAMPLE)
	_unur_hrd_debug_sample( gen, X, i );
#endif
      return X;
    }
    if (hrx > 0.)
      lambda = hrx;   
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not valid");
      return UNUR_INFINITY;
    } 
  }
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_hrd_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HRD_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = HRD (Hazard Rate Decreasing)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_hrd_sample",gen->genid);
  if (gen->variant & HRD_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_hrd_debug_sample( const struct unur_gen *gen, double x, int i )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HRD_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: X = %g\t #iterations = %d\n",gen->genid,x,i);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_hrd_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = HR\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: HRD (Hazard Rate Decreasing)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E[#iterations] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");
    if (help) {
      _unur_string_append(info,"parameters:\n");
      if (gen->variant & HRD_VARFLAG_VERIFY)
	_unur_string_append(info,"   verify = on\n");
      _unur_string_append(info,"\n");
    }
} 
#endif   
