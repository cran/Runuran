/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/matr.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <methods/dgt.h>
#include <methods/dstd.h>
#include <methods/dstd_struct.h>
#include <methods/hinv.h>
#include <methods/mixt.h>
#include <methods/mixt_struct.h>
#include <methods/ninv.h>
#include <methods/pinv.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
struct unur_gen *unur_init( struct unur_par *par )
{                
  _unur_check_NULL(NULL,par,NULL);
  return (par->init(par));
} 
int unur_reinit( struct unur_gen *gen )
{
  int status = UNUR_SUCCESS;
  _unur_check_NULL(NULL,gen,UNUR_ERR_NULL);
  if (gen->reinit) {
    status = gen->reinit(gen);
    if (status == UNUR_SUCCESS) return status;
  }
  else {
    _unur_error(gen->genid,UNUR_ERR_NO_REINIT,"");
    status = UNUR_ERR_NO_REINIT;
  }
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    gen->sample.discr = _unur_sample_discr_error;
    break;
  case UNUR_METH_CONT:
  case UNUR_METH_CEMP:
    gen->sample.cont = _unur_sample_cont_error;
    break;
  case UNUR_METH_VEC:
  case UNUR_METH_CVEMP:
    gen->sample.cvec = _unur_sample_cvec_error;
    break;
  case UNUR_METH_MAT:
    gen->sample.matr = _unur_sample_matr_error;
    break;
  default:
    _unur_error("reinit",UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  }
  return status;
} 
int
unur_sample_discr( struct unur_gen *gen )
{
  CHECK_NULL(gen,0);
  return (gen->sample.discr(gen));
} 
double
unur_sample_cont( struct unur_gen *gen )
{
  CHECK_NULL(gen,INFINITY);
  return (gen->sample.cont(gen));
} 
int
unur_sample_vec( struct unur_gen *gen, double *vector )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  return (gen->sample.cvec(gen,vector));
} 
int
unur_sample_matr( struct unur_gen *gen, double *matrix )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  return (gen->sample.matr(gen,matrix));
} 
double
unur_quantile ( struct unur_gen *gen, double U )
{
  CHECK_NULL(gen,FALSE);
  switch (gen->method) {
  case UNUR_METH_HINV:
    return unur_hinv_eval_approxinvcdf(gen,U);
  case UNUR_METH_NINV:
    return unur_ninv_eval_approxinvcdf(gen,U);
  case UNUR_METH_PINV:
    return unur_pinv_eval_approxinvcdf(gen,U);
  case UNUR_METH_CSTD:
    if (((struct unur_cstd_gen*)gen->datap)->is_inversion)
      return unur_cstd_eval_invcdf(gen,U);
    break;
  case UNUR_METH_MIXT:
    if (((struct unur_mixt_gen*)gen->datap)->is_inversion)
      return unur_mixt_eval_invcdf(gen,U);
    break;
  case UNUR_METH_DGT:
    return ((double) unur_dgt_eval_invcdf(gen,U));
  case UNUR_METH_DSTD:
    if (((struct unur_dstd_gen*)gen->datap)->is_inversion)
      return unur_dstd_eval_invcdf(gen,U);
    break;
  }
  _unur_error(gen->genid,UNUR_ERR_NO_QUANTILE,"");
  return UNUR_INFINITY;
} 
int
_unur_gen_is_inversion ( struct unur_gen *gen )
{
  CHECK_NULL(gen,FALSE);
  switch (gen->method) {
  case UNUR_METH_HINV:
  case UNUR_METH_NINV:
  case UNUR_METH_PINV:
  case UNUR_METH_DGT:
    return TRUE;
  case UNUR_METH_CSTD:
    return (((struct unur_cstd_gen*)gen->datap)->is_inversion);
  case UNUR_METH_MIXT:
    return (((struct unur_mixt_gen*)gen->datap)->is_inversion);
  default:
    return FALSE;
  }
} 
int
_unur_sample_discr_error( struct unur_gen *gen ATTRIBUTE__UNUSED )
{
  unur_errno = UNUR_ERR_GEN_CONDITION;
  return 0;
} 
double
_unur_sample_cont_error( struct unur_gen *gen ATTRIBUTE__UNUSED )
{
  unur_errno = UNUR_ERR_GEN_CONDITION;
  return INFINITY;
} 
int
_unur_sample_cvec_error( struct unur_gen *gen, double *vec )
{ 
  int d;
  unur_errno = UNUR_ERR_GEN_CONDITION;
  for (d=0; d<(gen->distr->dim); d++) vec[d] = INFINITY;
  return UNUR_FAILURE;
} 
int
_unur_sample_matr_error( struct unur_gen *gen, double *mat )
{ 
  int n_rows, n_cols, dim, j;
  unur_errno = UNUR_ERR_GEN_CONDITION;
  unur_distr_matr_get_dim(gen->distr, &n_rows, &n_cols );
  dim = n_rows * n_cols;
  for (j=0; j<dim; j++)
    mat[j] = INFINITY;
  return UNUR_FAILURE;
} 
void
unur_free( struct unur_gen *gen )
{                
  if (gen) gen->destroy(gen);
} 
const char *
unur_gen_info( struct unur_gen *gen, int help )
{
#ifdef UNUR_ENABLE_INFO
  _unur_check_NULL("",gen,NULL);
  if (gen->info) {
    if (gen->infostr == NULL) 
      gen->infostr = _unur_string_new();
    else 
      _unur_string_clear(gen->infostr);
    gen->info((struct unur_gen*) gen, help);
    return gen->infostr->text;
  }
  else {
    return NULL;
  }
#else
  return "INFO string not enable";
#endif
} 
int
unur_get_dimension( const struct unur_gen *gen )
{
  CHECK_NULL(gen,0);
  return (gen->distr->dim);
} 
const char *
unur_get_genid( const struct unur_gen *gen )
{
  CHECK_NULL(gen,NULL);
  return gen->genid;
} 
struct unur_distr *
unur_get_distr( const struct unur_gen *gen )
{
  CHECK_NULL(gen,NULL);
  return gen->distr;
} 
int 
unur_set_use_distr_privatecopy( struct unur_par *par, int use_privatecopy )
{
  _unur_check_NULL("",par,UNUR_ERR_NULL);
  par->distr_is_privatecopy = use_privatecopy;
  return UNUR_SUCCESS;
} 
struct unur_gen *
unur_gen_clone( const struct unur_gen *gen )
{
  _unur_check_NULL( "Clone", gen, NULL );
  _unur_check_NULL( "Clone", gen->clone, NULL );
  return (gen->clone(gen));
} 
struct unur_par *
_unur_par_new( size_t s)
{
  struct unur_par *par = _unur_xmalloc( sizeof(struct unur_par) );
  par->datap = _unur_xmalloc(s);
  par->s_datap = s;
  par->distr_is_privatecopy = TRUE;   
  return par;
} 
struct unur_par *
_unur_par_clone( const struct unur_par *par )
{
  struct unur_par *clone;
  _unur_check_NULL("clone", par, NULL);
  clone = _unur_xmalloc( sizeof(struct unur_par) );
  memcpy (clone, par, sizeof(struct unur_par));
  clone->datap = _unur_xmalloc(par->s_datap);
  memcpy (clone->datap, par->datap, par->s_datap);
  return clone;
} 
void 
unur_par_free( struct unur_par *par)
{
  _unur_check_NULL("free", par, RETURN_VOID );
  _unur_par_free(par);
} 
struct unur_gen *
_unur_generic_create( struct unur_par *par, size_t s )
{
  struct unur_gen *gen;
  gen = _unur_xmalloc( sizeof(struct unur_gen) );
  gen->datap = _unur_xmalloc(s);
  gen->s_datap = s;
  gen->distr_is_privatecopy = par->distr_is_privatecopy;
  if (gen->distr_is_privatecopy) 
    gen->distr = (par->distr) ? _unur_distr_clone(par->distr) : NULL;
  else
    gen->distr = (struct unur_distr *) par->distr;
  gen->destroy = NULL;               
  gen->clone = NULL;                
  gen->reinit = NULL;                
  gen->method = par->method;        
  gen->variant = par->variant;      
  gen->set = par->set;              
  gen->debug = par->debug;          
  gen->urng = par->urng;            
  gen->urng_aux = par->urng_aux;    
  gen->gen_aux = NULL;              
  gen->gen_aux_list = NULL;         
  gen->n_gen_aux_list = 0;
  gen->status = UNUR_FAILURE;       
#ifdef UNUR_ENABLE_INFO
  gen->infostr = NULL;              
  gen->info = NULL;                 
#endif
  return gen;
} 
struct unur_gen *
_unur_generic_clone( const struct unur_gen *gen, const char *type )
{ 
  struct unur_gen *clone;
  clone = _unur_xmalloc( sizeof(struct unur_gen) );
  memcpy( clone, gen, sizeof(struct unur_gen) );
  clone->datap = _unur_xmalloc(gen->s_datap);
  memcpy (clone->datap, gen->datap, gen->s_datap);
  clone->genid = _unur_set_genid(type);
#ifdef UNUR_ENABLE_INFO
  clone->infostr = NULL;
#endif
  clone->distr_is_privatecopy = gen->distr_is_privatecopy;
  if (clone->distr_is_privatecopy) 
    clone->distr = (gen->distr) ? _unur_distr_clone(gen->distr) : NULL;
  else
    clone->distr = gen->distr;
  if (gen->gen_aux)
    clone->gen_aux = _unur_gen_clone( gen->gen_aux );
  if (gen->gen_aux_list && gen->n_gen_aux_list) {
    clone->gen_aux_list = _unur_gen_list_clone( gen->gen_aux_list, gen->n_gen_aux_list );
    clone->n_gen_aux_list = gen->n_gen_aux_list;
  }
  return clone;
} 
void
_unur_generic_free( struct unur_gen *gen )
{ 
  if (gen->gen_aux)
    _unur_free(gen->gen_aux);
  if (gen->gen_aux_list && gen->n_gen_aux_list)
    _unur_gen_list_free( gen->gen_aux_list, gen->n_gen_aux_list );
  if (gen->distr_is_privatecopy && gen->distr)
    _unur_distr_free( gen->distr );
  _unur_free_genid(gen);
  COOKIE_CLEAR(gen);
  free(gen->datap);
#ifdef UNUR_ENABLE_INFO
  if (gen->infostr) _unur_string_free(gen->infostr);  
#endif
  free(gen);
} 
struct unur_gen ** 
_unur_gen_list_set( struct unur_gen *gen, int n_gen_list )
{
  struct unur_gen **gen_list;
  int i;
  _unur_check_NULL( "gen_list_set", gen, NULL );
  if (n_gen_list < 1) {
    _unur_error("gen_list_set",UNUR_ERR_PAR_SET,"dimension < 1");
    return NULL;
  }
  gen_list = _unur_xmalloc (n_gen_list * sizeof(struct unur_gen *));
  for (i=0; i<n_gen_list; i++)
    gen_list[i] = gen;
  return gen_list;
} 
struct unur_gen ** 
_unur_gen_list_clone( struct unur_gen **gen_list, int n_gen_list )
{
  struct unur_gen **clone_list;
  int i;
  _unur_check_NULL( "gen_list_clone", gen_list, NULL );
  if (n_gen_list < 1) {
    _unur_error("gen_list_clone",UNUR_ERR_PAR_SET,"dimension < 1");
    return NULL;
  }
  for (i=0; i<n_gen_list; i++)
    _unur_check_NULL( "gen_list_clone", *(gen_list+i), NULL );
  clone_list = _unur_xmalloc (n_gen_list * sizeof(struct unur_gen *));
  if (n_gen_list > 1 && gen_list[0] == gen_list[1]) {
      clone_list[0] = _unur_gen_clone( gen_list[0] );
      for (i=0; i<n_gen_list; i++)
	clone_list[i] = clone_list[0];
  }
  else {
    for (i=0; i<n_gen_list; i++)
      clone_list[i] = _unur_gen_clone( gen_list[i] );
  }
  return clone_list;
} 
void
_unur_gen_list_free( struct unur_gen **gen_list, int n_gen_list )
{
  int i, i2, imax;
  if (gen_list==NULL) 
    return; 
  if (n_gen_list < 1) {
    _unur_error("gen_list_free",UNUR_ERR_PAR_SET,"dimension < 1");
    return;
  }
  i2 = (n_gen_list>1) ? 1 : 0; 
  imax = (gen_list[0] == gen_list[i2]) ? 1 : n_gen_list;
  for (i=0; i<imax; i++)
    if (gen_list[i]) _unur_free(gen_list[i]);
  free (gen_list);
} 
