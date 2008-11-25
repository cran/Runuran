/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <distributions/unur_stddistr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dstd.h"
#include "dstd_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define DSTD_DEBUG_REINIT    0x00000010u   
#define DSTD_DEBUG_CHG       0x00001000u   
#define DSTD_SET_VARIANT          0x01u
#define GENTYPE "DSTD"         
static struct unur_gen *_unur_dstd_init( struct unur_par *par );
static int _unur_dstd_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dstd_create( struct unur_par *par );
static int _unur_dstd_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_dstd_clone( const struct unur_gen *gen );
static void _unur_dstd_free( struct unur_gen *gen);
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dstd_debug_init( const struct unur_gen *gen );
static void _unur_dstd_debug_chg_pmfparams( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dstd_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr     
#define PAR       ((struct unur_dstd_par*)par->datap) 
#define GEN       ((struct unur_dstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.discr     
struct unur_par *
unur_dstd_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  if (!(distr->id & UNUR_DISTR_STD) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"standard distribution");
    return NULL;
  }
  if (DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"init() for special generators");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_dstd_par) );
  COOKIE_SET(par,CK_DSTD_PAR);
  par->distr    = distr;            
  par->method   = UNUR_METH_DSTD;   
  par->variant  = 0u;               
  par->set      = 0u;                   
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_dstd_init;
  return par;
} 
int 
unur_dstd_set_variant( struct unur_par *par, unsigned variant )
{
  unsigned old_variant;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, par->distr, UNUR_ERR_NULL );
  _unur_check_par_object( par, DSTD );
  old_variant = par->variant;
  par->variant = variant;
  if (par->DISTR_IN.init != NULL && par->DISTR_IN.init(par,NULL)==UNUR_SUCCESS ) {
    par->set |= DSTD_SET_VARIANT;    
    return UNUR_SUCCESS;
  }
  _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
  par->variant = old_variant;
  return UNUR_ERR_PAR_VARIANT;
} 
struct unur_gen *
_unur_dstd_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if (par->DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,"");
    return NULL;
  }
  if ( par->method != UNUR_METH_DSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_DSTD_PAR,NULL);
  gen = _unur_dstd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_dstd_check_par(gen) != UNUR_SUCCESS) {
    _unur_dstd_free(gen); return NULL;
  }
  GEN->is_inversion = FALSE;   
  if ( DISTR.init(NULL,gen)!=UNUR_SUCCESS ) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    _unur_dstd_free(gen); return NULL; 
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_dstd_debug_init(gen);
#endif
  return gen;
} 
int
_unur_dstd_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_dstd_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  GEN->is_inversion = FALSE;   
  if ( DISTR.init(NULL,gen)!=UNUR_SUCCESS ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return UNUR_ERR_GEN_DATA;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & DSTD_DEBUG_REINIT)
      _unur_dstd_debug_chg_pmfparams( gen );
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dstd_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSTD_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_dstd_gen) );
  COOKIE_SET(gen,CK_DSTD_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = NULL;    
  gen->destroy = _unur_dstd_free;
  gen->clone = _unur_dstd_clone;
  gen->reinit = _unur_dstd_reinit;
  GEN->gen_param = NULL;  
  GEN->n_gen_param = 0;   
  GEN->gen_iparam = NULL; 
  GEN->n_gen_iparam = 0;
  GEN->is_inversion = FALSE;    
  GEN->sample_routine_name = NULL ;  
  GEN->umin        = 0;    
  GEN->umax        = 1;    
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dstd_info;
#endif
  return gen;
} 
int
_unur_dstd_check_par( struct unur_gen *gen )
{
  if (!(gen->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed");
    return UNUR_ERR_GEN_DATA;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dstd_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_dstd_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DSTD_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  if (GEN->gen_param) {
    CLONE->gen_param = _unur_xmalloc( GEN->n_gen_param * sizeof(double) );
    memcpy( CLONE->gen_param, GEN->gen_param, GEN->n_gen_param * sizeof(double) );
  }
  if (GEN->gen_iparam) {
    CLONE->gen_iparam = _unur_xmalloc( GEN->n_gen_iparam * sizeof(int) );
    memcpy( CLONE->gen_iparam, GEN->gen_iparam, GEN->n_gen_iparam * sizeof(int) );
  }
  return clone;
#undef CLONE
} 
void
_unur_dstd_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  COOKIE_CHECK(gen,CK_DSTD_GEN,RETURN_VOID);
  if ( gen->method != UNUR_METH_DSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }
  SAMPLE = NULL;   
  if (GEN->gen_param)   free(GEN->gen_param);
  if (GEN->gen_iparam)  free(GEN->gen_iparam);
  _unur_generic_free(gen);
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_dstd_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSTD_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );
  fprintf(LOG,"%s: sampling routine = ",gen->genid);
  if (GEN->sample_routine_name)
    fprintf(LOG,"%s()",GEN->sample_routine_name);
  else
    fprintf(LOG,"(Unknown)");
  if (GEN->is_inversion)
    fprintf(LOG,"   (Inversion)");
  fprintf(LOG,"\n%s:\n",gen->genid);
  if (!(gen->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    fprintf(LOG,"%s: domain has been changed. U in (%g,%g)\n",gen->genid,GEN->umin,GEN->umax);
    fprintf(LOG,"%s:\n",gen->genid);
  }
} 
void 
_unur_dstd_debug_chg_pmfparams( const struct unur_gen *gen )
{
  FILE *LOG;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSTD_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: parameters of distribution changed:\n",gen->genid);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\tparam[%d] = %g\n",gen->genid,i,DISTR.params[i]);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dstd_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: DSTD (special generator for Discrete STandarD distribution)\n");
  _unur_string_append(info,"   variant = %d  %s\n", gen->variant,
		      (GEN->is_inversion)?"[implements inversion method]" : "");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#urn] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   variant = %d  %s\n", gen->variant,
			(gen->set & DSTD_SET_VARIANT) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }
} 
#endif   
