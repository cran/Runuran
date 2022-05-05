/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distributions/unur_stddistr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "mvstd.h"
#include "mvstd_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define MVSTD_DEBUG_REINIT   0x00000010u   
#define GENTYPE "MVSTD"        
static struct unur_gen *_unur_mvstd_init( struct unur_par *par );
static int _unur_mvstd_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_mvstd_create( struct unur_par *par );
static int _unur_mvstd_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_mvstd_clone( const struct unur_gen *gen );
static void _unur_mvstd_free( struct unur_gen *gen);
#ifdef UNUR_ENABLE_LOGGING
static void _unur_mvstd_debug_init( struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_mvstd_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cvec      
#define PAR       ((struct unur_mvstd_par*)par->datap) 
#define GEN       ((struct unur_mvstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec 
#define SAMPLE    gen->sample.cvec      
struct unur_par *
unur_mvstd_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL(GENTYPE,distr,NULL);
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);
  if (distr->id == UNUR_DISTR_GENERIC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"standard distribution");
    return NULL;
  }
  if (DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"init() for special generators");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_mvstd_par) );
  COOKIE_SET(par,CK_MVSTD_PAR);
  par->distr    = distr;            
  par->method   = UNUR_METH_MVSTD;  
  par->variant  = 0u;               
  par->set      = 0u;                   
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_mvstd_init;
  return par;
} 
struct unur_gen *
_unur_mvstd_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if (par->DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,"");
    return NULL;
  }
  if ( par->method != UNUR_METH_MVSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_MVSTD_PAR,NULL);
  gen = _unur_mvstd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if ( DISTR.init(gen)!=UNUR_SUCCESS ) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    _unur_mvstd_free(gen); return NULL; 
  }
  if (_unur_mvstd_check_par(gen) != UNUR_SUCCESS) {
    _unur_mvstd_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_mvstd_debug_init(gen);
#endif
  return gen;
} 
int
_unur_mvstd_reinit( struct unur_gen *gen )
{
  if ( DISTR.init(gen)!=UNUR_SUCCESS ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return UNUR_ERR_GEN_DATA;
  }
  return _unur_mvstd_check_par(gen);
} 
struct unur_gen *
_unur_mvstd_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MVSTD_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_mvstd_gen) );
  COOKIE_SET(gen,CK_MVSTD_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = NULL;      
  gen->destroy = _unur_mvstd_free;
  gen->clone = _unur_mvstd_clone;
  gen->reinit = _unur_mvstd_reinit;
  GEN->sample_routine_name = NULL ;  
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_mvstd_info;
#endif
  return gen;
} 
int
_unur_mvstd_check_par( struct unur_gen *gen )
{
  if (gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"truncated domain");
    return UNUR_ERR_GEN_CONDITION;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_mvstd_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_mvstd_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MVSTD_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_mvstd_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  COOKIE_CHECK(gen,CK_MVSTD_GEN,RETURN_VOID);
  if ( gen->method != UNUR_METH_MVSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_mvstd_debug_init( struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVSTD_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cvec_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = ",gen->genid);
  if (GEN->sample_routine_name)
    fprintf(LOG,"%s()\n",GEN->sample_routine_name);
  else
    fprintf(LOG,"(Unknown)\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_mvstd_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  int dim = gen->distr->dim;
  int samplesize = 10000;
  double E_urn;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",dim);
  _unur_distr_cvec_info_domain(gen);
  _unur_string_append(info,"\n\n");
  _unur_string_append(info,"method: MVSTD (special generator for MultiVariate continuous STandarD distribution)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  E_urn = unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize);
  _unur_string_append(info,"   E [#urn] = %.2f x %d = %.2f  [approx.]\n",
		      E_urn / dim, dim, E_urn);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
