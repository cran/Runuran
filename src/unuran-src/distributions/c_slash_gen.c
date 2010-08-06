/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <methods/x_gen_source.h>
#include <distr/distr_source.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
#include "unur_distributions.h"
inline static int slash_slash_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
int 
_unur_stdgen_slash_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_slash_slash );
    return slash_slash_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define NORMAL  gen->gen_aux   
inline static int
slash_slash_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (NORMAL==NULL) {
    struct unur_distr *distr = unur_distr_normal(NULL,0);
    struct unur_par *par = unur_cstd_new( distr );
    NORMAL = (par) ? _unur_init(par) : NULL;
    _unur_check_NULL( NULL, NORMAL, UNUR_ERR_NULL );
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
    _unur_distr_free( distr );
  }
  return UNUR_SUCCESS;
} 
double
_unur_stdgen_sample_slash_slash( struct unur_gen *gen )
{
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);
  return (_unur_sample_cont(NORMAL) / uniform());
} 
#undef NORMAL
