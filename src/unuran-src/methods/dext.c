/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dext.h"
#include "dext_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define GENTYPE "DEXT"         
static struct unur_gen *_unur_dext_init( struct unur_par *par );
static int _unur_dext_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dext_create( struct unur_par *par );
static struct unur_gen *_unur_dext_clone( const struct unur_gen *gen );
static void _unur_dext_free( struct unur_gen *gen);
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dext_debug_init( struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dext_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr     
#define PAR       ((struct unur_dext_par*)par->datap) 
#define GEN       ((struct unur_dext_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define SAMPLE    gen->sample.discr     
struct unur_par *
unur_dext_new( const struct unur_distr *distr )
{
  struct unur_par *par;
  if (distr != NULL) {
    if (distr->type != UNUR_DISTR_DISCR) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
    COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  }
  par = _unur_par_new( sizeof(struct unur_dext_par) );
  COOKIE_SET(par,CK_DEXT_PAR);
  par->distr    = distr;            
  PAR->init     = NULL;   
  PAR->sample   = NULL;   
  par->method   = UNUR_METH_DEXT;   
  par->variant  = 0u;               
  par->set      = 0u;               
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_dext_init;
  return par;
} 
int
unur_dext_set_init( struct unur_par *par, int (*init)(struct unur_gen *gen) )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DEXT );
  PAR->init = init;
  return UNUR_SUCCESS;
} 
int
unur_dext_set_sample( struct unur_par *par, int (*sample)(struct unur_gen *gen) )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, sample, UNUR_ERR_NULL );
  _unur_check_par_object( par, DEXT );
  PAR->sample = sample;
  return UNUR_SUCCESS;
} 
void *
unur_dext_get_params( struct unur_gen *gen, size_t size )
{
  _unur_check_NULL( GENTYPE, gen, NULL );
  COOKIE_CHECK(gen, CK_DEXT_GEN, NULL);
  if (size && size != GEN->size_param) {
    GEN->param = _unur_xrealloc(GEN->param, size);
    GEN->size_param = size;
  }
  return GEN->param;
} 
double  *
unur_dext_get_distrparams( struct unur_gen *gen )
{
  CHECK_NULL(gen, NULL);
  COOKIE_CHECK(gen, CK_DEXT_GEN, NULL);
  return DISTR.params;
} 
int
unur_dext_get_ndistrparams( struct unur_gen *gen )
{
  CHECK_NULL(gen, 0);
  COOKIE_CHECK(gen, CK_DEXT_GEN, 0);
  return DISTR.n_params;
} 
struct unur_gen *
_unur_dext_init( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_DEXT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_DEXT_PAR,NULL);
  if (PAR->sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_CONDITION,"sampling routine missing");
    return NULL;
  }
  gen = _unur_dext_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (GEN->init != NULL) {
    if (GEN->init(gen) != UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_FAILURE,"init for external generator failed");
      _unur_dext_free(gen); return NULL;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_dext_debug_init(gen);
#endif
  return gen;
} 
int
_unur_dext_reinit( struct unur_gen *gen )
{
  if (GEN->init != NULL) {
    if (GEN->init(gen) != UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_FAILURE,"init for external generator failed");
      return UNUR_FAILURE;
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dext_create( struct unur_par *par )
{
  struct unur_gen *gen;
  struct unur_distr *distr = NULL;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DEXT_PAR,NULL);
  if (par->distr == NULL)
    par->distr = distr = unur_distr_discr_new();
  gen = _unur_generic_create( par, sizeof(struct unur_dext_gen) );
  COOKIE_SET(gen,CK_DEXT_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = PAR->sample;      
  gen->destroy = _unur_dext_free;
  gen->clone = _unur_dext_clone;
  gen->reinit = _unur_dext_reinit;
  GEN->init = PAR->init;
  GEN->sample = PAR->sample;
  GEN->param    = NULL;   
  GEN->size_param  = 0;   
  if (distr) _unur_distr_free(distr);
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dext_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_dext_clone( const struct unur_gen *gen )
{
#define CLONE  ((struct unur_dext_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DEXT_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  if (GEN->param) {
    CLONE->param = _unur_xmalloc( GEN->size_param );
    memcpy( CLONE->param, GEN->param, GEN->size_param );
  }
  return clone;
#undef CLONE
} 
void
_unur_dext_free( struct unur_gen *gen )
{
  if( !gen ) 
    return;
  COOKIE_CHECK(gen,CK_DEXT_GEN,RETURN_VOID);
  if ( gen->method != UNUR_METH_DEXT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }
  SAMPLE = NULL;   
  if (GEN->param)  free(GEN->param);
  _unur_generic_free(gen);
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_dext_debug_init( struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DEXT_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = wrapper for external generator\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dext_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: DEXT (wrapper for Discrete EXTernal generators)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#urn] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
