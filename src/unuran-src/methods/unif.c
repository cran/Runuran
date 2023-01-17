/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "unif.h"
#include "unif_struct.h"
#define GENTYPE "UNIF"         
static struct unur_gen *_unur_unif_init( struct unur_par *par );
static int _unur_unif_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_unif_create( struct unur_par *par );
static struct unur_gen *_unur_unif_clone( const struct unur_gen *gen );
static void _unur_unif_free( struct unur_gen *gen);
static double _unur_unif_sample( struct unur_gen *gen );
#ifdef UNUR_ENABLE_INFO
static void _unur_unif_info( struct unur_gen *gen, int help );
#endif
#define PAR       ((struct unur_unif_par*)par->datap) 
#define GEN       ((struct unur_unif_gen*)gen->datap) 
#define SAMPLE  gen->sample.cont
#define _unur_unif_getSAMPLE(gen)  ( _unur_unif_sample )
struct unur_par *
unur_unif_new( const struct unur_distr *distr_dummy )
{ 
  struct unur_par *par;
  if (distr_dummy != NULL) {
    if (distr_dummy->type != UNUR_DISTR_CONT) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
    COOKIE_CHECK(distr_dummy,CK_DISTR_CONT,NULL);
  }
  par = _unur_par_new( sizeof(struct unur_unif_par) );
  COOKIE_SET(par,CK_UNIF_PAR);
  par->distr    = distr_dummy;     
  par->method   = UNUR_METH_UNIF;  
  par->variant  = 0u;              
  par->set      = 0u;                  
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_unif_init;
  return par;
} 
struct unur_gen *
_unur_unif_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_UNIF ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_UNIF_PAR,NULL);
  gen = _unur_unif_create(par);
  _unur_par_free(par);
  return gen;
} 
int
_unur_unif_reinit( struct unur_gen *gen )
{
  SAMPLE = _unur_unif_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
static struct unur_gen *
_unur_unif_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_UNIF_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_unif_gen) );
  COOKIE_SET(gen,CK_UNIF_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_unif_getSAMPLE(gen);
  gen->destroy = _unur_unif_free;
  gen->clone = _unur_unif_clone;
  gen->reinit = _unur_unif_reinit;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_unif_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_unif_clone( const struct unur_gen *gen )
{ 
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_UNIF_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
} 
void
_unur_unif_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_UNIF ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_UNIF_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_unif_sample( struct unur_gen *gen )
{ 
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_UNIF_GEN,UNUR_INFINITY);
  return _unur_call_urng(gen->urng);
} 
#ifdef UNUR_ENABLE_LOGGING
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_unif_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution: uniform (0,1)\n\n");
  _unur_string_append(info,"method: UNIF (wrapper for UNIForm random number generator)\n\n");
  if (help) {
    _unur_string_append(info,"[Remark: allows using uniform random number generator in UNU.RAN framework]\n");
  }
} 
#endif   
