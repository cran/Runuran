/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include <utils/unur_fp_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "dgt.h"
#include "dgt_struct.h"
#include "mixt.h"
#include "mixt_struct.h"
#define MIXT_VARFLAG_INVERSION   0x004u    
#define MIXT_SET_USEINVERSION     0x001u    
#define GENTYPE "MIXT"          
static struct unur_gen *_unur_mixt_init( struct unur_par *par );
static struct unur_gen *_unur_mixt_create( struct unur_par *par );
static int _unur_mixt_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_mixt_clone( const struct unur_gen *gen );
static void _unur_mixt_free( struct unur_gen *gen);
static double _unur_mixt_sample( struct unur_gen *gen );
static double _unur_mixt_sample_inv( struct unur_gen *gen );
static struct unur_gen *_unur_mixt_indexgen( const double *prob, int n_prob );
static int _unur_mixt_get_boundary( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_mixt_debug_init( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_mixt_info( struct unur_gen *gen, int help );
#endif
#define PAR       ((struct unur_mixt_par*)par->datap) 
#define GEN       ((struct unur_mixt_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define SAMPLE    gen->sample.cont      
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define INDEX     gen_aux
#define PROB      gen_aux->distr->data.discr.pv
#define COMP      gen_aux_list
#define N_COMP    n_gen_aux_list
#define _unur_mixt_getSAMPLE(gen) \
   ( ((gen)->variant & MIXT_VARFLAG_INVERSION) \
     ? _unur_mixt_sample_inv : _unur_mixt_sample )
struct unur_par *
unur_mixt_new( int n, const double *prob, struct unur_gen **comp )
{
  struct unur_par *par;
  _unur_check_NULL( GENTYPE, prob, NULL );
  _unur_check_NULL( GENTYPE, comp, NULL );
  if (n<1) { _unur_error(GENTYPE,UNUR_ERR_DISTR_DOMAIN,"n < 1"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_mixt_par) );
  COOKIE_SET(par,CK_MIXT_PAR);
  par->distr    = NULL;      
  PAR->n_comp   = n;         
  PAR->prob     = prob;      
  PAR->comp     = comp;      
  par->method   = UNUR_METH_MIXT;   
  par->variant  = 0u;               
  par->set      = 0u;               
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_mixt_init;
  return par;
} 
int
unur_mixt_set_useinversion( struct unur_par *par, int useinversion )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MIXT );
  par->variant = (useinversion)
    ? (par->variant | MIXT_VARFLAG_INVERSION)
    : (par->variant & (~MIXT_VARFLAG_INVERSION));
  par->set |= MIXT_SET_USEINVERSION;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_mixt_init( struct unur_par *par )
{
  struct unur_gen *gen;
  int i;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_MIXT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MIXT_PAR,NULL);
  gen = _unur_mixt_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }
  gen->INDEX = _unur_mixt_indexgen(PAR->prob,PAR->n_comp);
  gen->N_COMP = PAR->n_comp;    
  gen->COMP = _unur_xmalloc( gen->N_COMP * sizeof(struct unur_gen *));
  for (i=0; i<gen->N_COMP; i++)
    gen->COMP[i] = unur_gen_clone(PAR->comp[i]);
  _unur_par_free(par);
  if (_unur_mixt_check_par(gen) != UNUR_SUCCESS) {
    _unur_mixt_free(gen); return NULL;
  }
  if ( _unur_mixt_get_boundary(gen) != UNUR_SUCCESS ) {
    _unur_mixt_free(gen); return NULL;
  }
  unur_distr_set_name(gen->distr, "(mixture)");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_mixt_debug_init(gen);
#endif
  return gen;
} 
struct unur_gen *
_unur_mixt_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MIXT_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_mixt_gen) );
  COOKIE_SET(gen,CK_MIXT_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  gen->distr = unur_distr_cont_new();
  SAMPLE = _unur_mixt_getSAMPLE(gen);
  gen->destroy = _unur_mixt_free;
  gen->clone = _unur_mixt_clone;
  gen->reinit = NULL;    
  GEN->is_inversion = (gen->variant & MIXT_VARFLAG_INVERSION) ? TRUE : FALSE;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_mixt_info;
#endif
  return gen;
} 
int
_unur_mixt_check_par( struct unur_gen *gen )
{
  int i;
  int type;
  if (gen->INDEX == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"invalid probabilities");
    return UNUR_ERR_GEN_DATA;
  }
  for (i=0; i<gen->N_COMP; i++) {
    if (gen->COMP[i] == NULL) {
      _unur_error(gen->genid,UNUR_ERR_NULL,"component is NULL");
      return UNUR_ERR_NULL;
    }
    type = gen->COMP[i]->method & UNUR_MASK_TYPE;
    if ( type != UNUR_METH_DISCR && 
	 type != UNUR_METH_CONT  &&
	 type != UNUR_METH_CEMP  ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"component not univariate");
      return UNUR_ERR_GEN_INVALID;
    }
    if (GEN->is_inversion && (! unur_gen_is_inversion (gen->COMP[i]))) {
      _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"component does not implement inversion");
      return UNUR_ERR_GEN_INVALID;
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_mixt_clone( const struct unur_gen *gen )
{
#define CLONE  ((struct unur_mixt_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MIXT_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_mixt_free( struct unur_gen *gen )
{
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_MIXT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MIXT_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_mixt_sample( struct unur_gen *gen )
{
  struct unur_gen *comp;
  int J;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_MIXT_GEN,UNUR_INFINITY);
  J = unur_sample_discr(gen->INDEX);
  comp = gen->COMP[J];
  switch(comp->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    return ((double) comp->sample.discr(comp));
  case UNUR_METH_CONT:
  case UNUR_METH_CEMP:
  default:
    return (comp->sample.cont(comp));
  }
} 
double
_unur_mixt_sample_inv( struct unur_gen *gen )
{
  double U, recycle;
  int J;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_MIXT_GEN,UNUR_INFINITY);
  U = _unur_call_urng(gen->urng);
  J =unur_dgt_eval_invcdf_recycle( gen->INDEX, U, &recycle );
  if (_unur_iszero(recycle)) recycle = DBL_MIN;
  if (_unur_isone(recycle))  recycle = 1. - DBL_EPSILON;
  return unur_quantile(gen->COMP[J], recycle);
} 
double
unur_mixt_eval_invcdf( const struct unur_gen *gen, double u )
{
  double recycle;
  int J;
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  if ( ! (gen->method == UNUR_METH_MIXT && GEN->is_inversion) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return UNUR_INFINITY;
  }
  COOKIE_CHECK(gen,CK_MIXT_GEN,UNUR_INFINITY);
  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  
  }
  J =unur_dgt_eval_invcdf_recycle( gen->INDEX, u, &recycle );
  if (_unur_iszero(recycle)) recycle = DBL_MIN;
  if (_unur_isone(recycle))  recycle = 1. - DBL_EPSILON;
  return unur_quantile(gen->COMP[J], recycle);
} 
struct unur_gen *
_unur_mixt_indexgen( const double *prob, int n_prob )
{
  struct unur_distr *distr;
  struct unur_par *par;
  struct unur_gen *igen;
  distr = unur_distr_discr_new();
  unur_distr_discr_set_pv(distr, prob, n_prob);
  par = unur_dgt_new(distr);
  igen = unur_init(par);
  unur_distr_free(distr);
  return igen;
} 
int
_unur_mixt_get_boundary( struct unur_gen *gen )
{
  int i;
  int overlap = FALSE;
  double comp_left, comp_right;
  double bd_left, bd_right;
  struct unur_gen *comp;
  bd_left = UNUR_INFINITY;
  bd_right = -UNUR_INFINITY;
  for (i=0; i<gen->N_COMP; i++) {
    comp = gen->COMP[i];
    CHECK_NULL(comp,UNUR_ERR_NULL);
    switch (comp->method & UNUR_MASK_TYPE) {
    case UNUR_METH_CONT:
      comp_left  = comp->distr->data.cont.BD_LEFT;
      comp_right = comp->distr->data.cont.BD_RIGHT;
      break;
    case UNUR_METH_DISCR:
      comp_left  = (double) (comp->distr->data.discr.BD_LEFT);
      comp_right = (double) (comp->distr->data.discr.BD_RIGHT);
      break;
    default:
      comp_left = -UNUR_INFINITY;
      comp_right = UNUR_INFINITY;
    }
    if ( _unur_FP_less(comp_left,bd_right) )
      overlap = TRUE;
    bd_left = _unur_min(bd_left, comp_left);
    bd_right = _unur_max(bd_right, comp_right);
  }
  if (GEN->is_inversion && overlap) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"domains of components overlap or are unsorted");
    return UNUR_ERR_GEN_INVALID;
  }
  unur_distr_cont_set_domain(gen->distr, bd_left, bd_right);
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_mixt_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  struct unur_gen *comp;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MIXT_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = MIXT (MIXTure of distributions -- meta method)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_mixt_sample",gen->genid);
  if (GEN->is_inversion) fprintf(LOG,"_inv");
  fprintf(LOG,"()\n%s:\n",gen->genid);
  fprintf(LOG,"%s: use inversion = %s",gen->genid,
	  (GEN->is_inversion) ? "on" : "off");
  _unur_print_if_default(gen,MIXT_SET_USEINVERSION);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: probabilities (%d) = \n",gen->genid, gen->N_COMP);
  fprintf(LOG,"%s:   %g",gen->genid, (gen->PROB)[0]);
  for (i=1; i<gen->N_COMP; i++)
    fprintf(LOG,", %g", (gen->PROB)[i]);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: components (%d):\n",gen->genid, gen->N_COMP);
  for (i=0; i<gen->N_COMP; i++) {
    comp = gen->COMP[i];
    fprintf(LOG,"%s:   [%d]: %s\n",gen->genid, i, comp->genid);
    fprintf(LOG,"%s:\t type = ",gen->genid); 
    switch (comp->distr->type) {
    case UNUR_DISTR_CONT:
    case UNUR_DISTR_CEMP:
      fprintf(LOG,"continuous\n");
      break;
    case UNUR_DISTR_DISCR:
      fprintf(LOG,"discrete\n");
      break;
    default:
      fprintf(LOG,"[unknown]\n");
    }
    fprintf(LOG,"%s:\t name = %s\n",gen->genid, comp->distr->name);
  }
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_mixt_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_gen *comp;
  int i;
  double sum;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   # components = %d\n", gen->N_COMP);
  if (help) {
    sum = ((struct unur_dgt_gen*)gen->INDEX->datap)->sum;
    _unur_string_append(info,"   probabilities = (%g", gen->PROB[0] / sum);
    for (i=1; i<gen->N_COMP; i++)
      _unur_string_append(info,", %g", gen->PROB[i] / sum);
    _unur_string_append(info,")\n");
    _unur_string_append(info,"   components = \n");
    for (i=0; i<gen->N_COMP; i++) {
      comp = gen->COMP[i];
      _unur_string_append(info,"\t[%d] %s - ",i, comp->genid);
      switch (comp->distr->type) {
      case UNUR_DISTR_CONT:
      case UNUR_DISTR_CEMP:
	_unur_string_append(info,"continuous");
	break;
      case UNUR_DISTR_DISCR:
	_unur_string_append(info,"discrete");
	break;
      default:
	_unur_string_append(info,"[unknown]");
      }
      _unur_string_append(info,": %s\n",comp->distr->name);
    }
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: MIXT (MIXTure of distributions -- meta method)\n");
  _unur_string_append(info,"   select component = method DGT\n");
  _unur_string_append(info,"   inversion method = %s\n",
		      (GEN->is_inversion) ? "TRUE" : "FALSE");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics: depends on components\n");
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   useinversion = ");
    if (gen->variant & MIXT_VARFLAG_INVERSION)
      _unur_string_append(info,"on\n");
    else
      _unur_string_append(info,"off  [default]\n");
  }
} 
#endif   
