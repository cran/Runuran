/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <limits.h>
#include <unur_source.h>
#include <methods/cstd.h>   
#include <methods/dstd_struct.h>
#include "unur_distributions_source.h"
inline static int zipf_zet_init( struct unur_gen *gen );
#define PAR       ((struct unur_dstd_par*)par->datap) 
#define GEN       ((struct unur_dstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define uniform()  _unur_call_urng(gen->urng) 
#define rho  (DISTR.params[0])    
#define tau  (DISTR.params[1])
int 
_unur_stdgen_zipf_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_zipf_zet );
    return zipf_zet_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS  (2)
#define c   (GEN->gen_param[0])
#define d   (GEN->gen_param[1])
inline static int
zipf_zet_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (rho<tau) {
    c = tau - 0.5;
    d = 0.;
  }
  else {
    c = rho - 0.5;
    d = (1. + rho) * log((1. + tau)/(1. + rho));
  }
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_zipf_zet( struct unur_gen *gen )
{
  double U, V, E, X;
  int K;
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  do {
    do {
      U = uniform();
      V = uniform();
      X = (c+0.5) * exp( -log(U)/rho ) - c;
    } while (X <= 0.5 || X >= (double) INT_MAX);
    K = (long int) (X+0.5);
    E = -log(V);
  } while ( E < (1.+rho) * log( (K+tau)/(X+c)) - d );
  return K;
} 
#undef GEN_N_PARAMS
#undef c
#undef d
