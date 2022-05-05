/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>   
#include <methods/dstd_struct.h>
#include "unur_distributions_source.h"
inline static int logarithmic_lsk_init( struct unur_gen *gen );
#define PAR       ((struct unur_dstd_par*)par->datap) 
#define GEN       ((struct unur_dstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define uniform()  _unur_call_urng(gen->urng) 
#define theta  (DISTR.params[0])    
int 
_unur_stdgen_logarithmic_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_logarithmic_lsk );
    return logarithmic_lsk_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS  (2)
#define t   (GEN->gen_param[0])
#define h   (GEN->gen_param[1])
#define theta_limit  0.97
inline static int
logarithmic_lsk_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
    t = 0.; h = 0.;  
  }
  if (theta < theta_limit)
    t = -theta / log(1.0 - theta);
  else
    h=log(1.0 - theta);
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_logarithmic_lsk( struct unur_gen *gen )
{
  double U, V, p, q;
  int K;
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  U = uniform();
  if (theta < theta_limit) {
    K = 1;
    p = t;
    while (U > p) {
      U -= p;
      K++;
      p *= theta * (K - 1.)/((double) K);
    }
    return K;
  }
  else {
    if (U > theta) 
      return 1;
    V = uniform();
    q = 1. - exp(V * h);
    if ( U <= q * q) {
      K = 1 + (int)(log(U)/log(q));
      return K;
    }
    return ((U > q) ? 1 : 2);
  }
} 
#undef GEN_N_PARAMS
#undef t
#undef h
#undef theta_limit
