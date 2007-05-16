/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define theta  (DISTR.params[0])    
#define lambda (DISTR.params[1])    
int 
_unur_stdgen_cauchy_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case UNUR_STDGEN_INVERSION:   
    if (gen) GEN->is_inversion = TRUE;
    _unur_cstd_set_sampling_routine(par,gen,_unur_stdgen_sample_cauchy_inv); 
    return UNUR_SUCCESS;
  default: 
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
} 
double _unur_stdgen_sample_cauchy_inv( struct unur_gen *gen )
{
  double U,X;
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);
  U = GEN->umin + uniform() * (GEN->umax-GEN->umin);
  X = tan( M_PI * (U - 0.5) );
  return ((DISTR.n_params==0) ? X : theta + lambda * X );
} 
