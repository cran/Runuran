/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <methods/x_gen.h>
#include <methods/cstd.h>
#include <methods/mvstd.h>
#include <methods/mvstd_struct.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
static int _unur_stdgen_init_multinormal_cholesky( struct unur_gen *gen );
#define PAR       ((struct unur_mvstd_par*)par->datap) 
#define GEN       ((struct unur_mvstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec 
int 
_unur_stdgen_multinormal_init( struct unur_gen *gen )
{
  if ( gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"truncated domain not allowed");
    return UNUR_FAILURE;
  }
  gen->sample.cvec = _unur_stdgen_sample_multinormal_cholesky;
  GEN->sample_routine_name = "_unur_stdgen_sample_multinormal_cholesky";
  return _unur_stdgen_init_multinormal_cholesky(gen);
} 
#define NORMAL  gen->gen_aux        
int
_unur_stdgen_init_multinormal_cholesky( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MVSTD_GEN,UNUR_ERR_COOKIE);
  if (NORMAL==NULL) {
    struct unur_distr *distr = unur_distr_normal(NULL,0);
    NORMAL = unur_init( unur_cstd_new( distr ) );
    _unur_check_NULL( gen->genid, NORMAL, UNUR_ERR_NULL );
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
    _unur_distr_free( distr );
  }
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_multinormal_cholesky( struct unur_gen *gen, double *X )
{
#define idx(a,b) ((a)*dim+(b))
  int j,k;
  int dim = gen->distr->dim;     
  double *L = DISTR.cholesky;    
  double *mean = DISTR.mean;     
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MVSTD_GEN,UNUR_ERR_COOKIE);
  for (j=0; j<dim; j++)
    X[j] = unur_sample_cont(NORMAL);
  for (k=dim-1; k>=0; k--) {
    X[k] *= L[idx(k,k)];
    for (j=k-1; j>=0; j--)
      X[k] += X[j] * L[idx(k,j)];
    X[k] += mean[k];
  }
  return UNUR_SUCCESS;
#undef idx
} 
#undef NORMAL
