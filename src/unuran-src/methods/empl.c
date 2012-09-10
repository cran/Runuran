/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "empl.h"
#include "empl_struct.h"
#define EMPL_DEBUG_PRINTDATA   0x00000100u
#define GENTYPE "EMPL"         
static struct unur_gen *_unur_empl_init( struct unur_par *par );
static struct unur_gen *_unur_empl_create( struct unur_par *par );
static struct unur_gen *_unur_empl_clone( const struct unur_gen *gen );
static void _unur_empl_free( struct unur_gen *gen);
static double _unur_empl_sample( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_empl_debug_init( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_empl_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cemp      
#define PAR       ((struct unur_empl_par*)par->datap) 
#define GEN       ((struct unur_empl_gen*)gen->datap) 
#define DISTR     gen->distr->data.cemp 
#define SAMPLE    gen->sample.cont           
inline static int 
compare_doubles (const void *a, const void *b)
{ 
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}
#define _unur_empl_getSAMPLE(gen)   (_unur_empl_sample)
struct unur_par *
unur_empl_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CEMP,NULL);
  if (DISTR_IN.sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"observed sample"); return NULL; }
  if (DISTR_IN.n_sample < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"number of observed sample"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_empl_par) );
  COOKIE_SET(par,CK_EMPL_PAR);
  par->distr    = distr;          
  par->method   = UNUR_METH_EMPL; 
  par->variant  = 0u;             
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init     = _unur_empl_init;
  return par;
} 
struct unur_gen *
_unur_empl_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_EMPL ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_EMPL_PAR,NULL);
  gen = _unur_empl_create(par);
  _unur_par_free(par);
  if (!gen) return NULL; 
  qsort( GEN->observ, (size_t)GEN->n_observ, sizeof(double), compare_doubles);
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_empl_debug_init(gen);
#endif
  return gen;
} 
struct unur_gen *
_unur_empl_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_EMPL_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_empl_gen) );
  COOKIE_SET(gen,CK_EMPL_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_empl_getSAMPLE(gen);
  gen->destroy = _unur_empl_free;
  gen->clone = _unur_empl_clone;
  GEN->observ   = DISTR.sample;          
  GEN->n_observ = DISTR.n_sample;        
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_empl_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_empl_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_empl_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_EMPL_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->observ = clone->distr->data.cemp.sample;   
  return clone;
#undef CLONE
} 
void
_unur_empl_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_EMPL ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_EMPL_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_empl_sample( struct unur_gen *gen )
{ 
  double U,X;
  int J;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_EMPL_GEN,UNUR_INFINITY);
  U = _unur_call_urng(gen->urng) * (GEN->n_observ-1);
  J = (int) (U);
  X = GEN->observ[J] + (U-J)*(GEN->observ[J+1] - GEN->observ[J]);
  return X;
} 
#ifdef UNUR_ENABLE_LOGGING
static void
_unur_empl_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_EMPL_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = EMPL (EMPirical distribution with Linear interpolation)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cemp_debug( gen->distr, gen->genid, (gen->debug & EMPL_DEBUG_PRINTDATA));
  fprintf(LOG,"%s: sampling routine = _unur_empl_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_empl_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = DATA  [length=%d]\n", GEN->n_observ);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: EMPL (EMPirical distribution with Linear interpolation)\n");
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
