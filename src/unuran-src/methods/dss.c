/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dss.h"
#include "dss_struct.h"
#define DSS_VARIANT_NONE       0x000u     
#define DSS_VARIANT_PV         0x001u     
#define DSS_VARIANT_PMF        0x002u     
#define DSS_VARIANT_CDF        0x004u     
#define DSS_DEBUG_REINIT        0x00000010u   
#define DSS_DEBUG_PRINTVECTOR   0x00000100u
#define GENTYPE "DSS"         
static struct unur_gen *_unur_dss_init( struct unur_par *par );
static int _unur_dss_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dss_create( struct unur_par *par );
static int _unur_dss_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_dss_clone( const struct unur_gen *gen );
static void _unur_dss_free( struct unur_gen *gen);
static int _unur_dss_sample( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dss_debug_init( struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dss_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr      
#define PAR       ((struct unur_dss_par*)par->datap) 
#define GEN       ((struct unur_dss_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define SAMPLE    gen->sample.discr     
#define PMF(x)    _unur_discr_PMF((x),(gen->distr))   
#define CDF(x)    _unur_discr_CDF((x),(gen->distr))   
#define _unur_dss_getSAMPLE(gen)   (_unur_dss_sample)
struct unur_par *
unur_dss_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  unsigned variant = DSS_VARIANT_NONE;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  if (DISTR_IN.pv && (distr->set & UNUR_DISTR_SET_PMFSUM))
    variant = DSS_VARIANT_PV;
  else if (DISTR_IN.pmf && (distr->set & UNUR_DISTR_SET_PMFSUM)) 
    variant = DSS_VARIANT_PMF;
  else if (DISTR_IN.cdf)
    variant = DSS_VARIANT_CDF;
  if (variant == DSS_VARIANT_NONE) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV+sum, PMF+sum, or CDF");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_dss_par) );
  COOKIE_SET(par,CK_DSS_PAR);
  par->distr       = distr;          
  par->method      = UNUR_METH_DSS;  
  par->variant     = variant;        
  par->set         = 0u;                 
  par->urng        = unur_get_default_urng(); 
  par->urng_aux    = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_dss_init;
  return par;
} 
struct unur_gen *
_unur_dss_init( struct unur_par *par )
{ 
  struct unur_gen *gen;         
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_DSS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DSS_PAR,NULL);
  gen = _unur_dss_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_dss_debug_init(gen);
#endif
  return gen;
} 
int
_unur_dss_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_dss_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  SAMPLE = _unur_dss_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & DSS_DEBUG_REINIT) _unur_dss_debug_init(gen);
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dss_create( struct unur_par *par )
{
  struct unur_gen *gen;       
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSS_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_dss_gen) );
  COOKIE_SET(gen,CK_DSS_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_dss_getSAMPLE(gen);
  gen->destroy = _unur_dss_free;
  gen->clone = _unur_dss_clone;
  gen->reinit = _unur_dss_reinit;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dss_info;
#endif
  return gen;
} 
int
_unur_dss_check_par( struct unur_gen *gen )
{
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    if (DISTR.pv != NULL) break;
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV");
    return UNUR_ERR_DISTR_REQUIRED;
  case DSS_VARIANT_PMF:
    if (DISTR.pmf != NULL) break;
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF");
    return UNUR_ERR_DISTR_REQUIRED;
  case DSS_VARIANT_CDF:
    if (DISTR.cdf != NULL) return UNUR_SUCCESS;
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF");
    return UNUR_ERR_DISTR_REQUIRED;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (unur_distr_discr_upd_pmfsum(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dss_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_dss_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DSS_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_dss_free( struct unur_gen *gen )
{ 
  if (!gen) 
    return;
  if ( gen->method != UNUR_METH_DSS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DSS_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
int
_unur_dss_sample( struct unur_gen *gen )
{ 
  int J;
  double U;
  double sum;
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DSS_GEN,INT_MAX);
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    U = DISTR.sum * _unur_call_urng(gen->urng);
    sum = 0.;
    for (J=0; J<DISTR.n_pv; J++) {
      sum += DISTR.pv[J];
      if (sum >= U) break;
    }
    return (J + DISTR.domain[0]);
  case DSS_VARIANT_PMF:
    U = DISTR.sum * _unur_call_urng(gen->urng);
    sum = 0.;
    for (J=DISTR.domain[0]; J<=DISTR.domain[1]; J++) {
      sum += PMF(J);
      if (sum >= U) break;
    }
    return J;
  case DSS_VARIANT_CDF:
    U = _unur_call_urng(gen->urng);
    for (J=DISTR.domain[0]; J<=DISTR.domain[1]; J++) {
      if (CDF(J) >= U) break;
    }
    return J;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INT_MAX;
  }
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_dss_debug_init( struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSS_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = sequential search\n",gen->genid);
  _unur_distr_discr_debug( gen->distr,gen->genid,(gen->debug & DSS_DEBUG_PRINTVECTOR));
  fprintf(LOG,"%s: sampling routine = _unur_dss_sample()\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    fprintf(LOG,"use PV\n");  break;
  case DSS_VARIANT_PMF:
    fprintf(LOG,"use PMF\n"); break;
  case DSS_VARIANT_CDF:
    fprintf(LOG,"use CDF\n"); break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  }
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dss_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  switch(gen->variant) {
  case DSS_VARIANT_PV:
    _unur_string_append(info,"   functions = PV  [length=%d]\n",DISTR.domain[1]-DISTR.domain[0]+1);
    break;
  case DSS_VARIANT_PMF:
    _unur_string_append(info,"   functions = PMF\n");
    break;
  case DSS_VARIANT_CDF:
    _unur_string_append(info,"   functions = CDF\n");
    break;
  }
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: DSS (Simple Sequential Search)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics: slow\n");
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
