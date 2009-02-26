/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define burr_type (gen->distr->id)
#define k         (DISTR.params[1])
#define c         (DISTR.params[2])
int 
_unur_stdgen_burr_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case UNUR_STDGEN_INVERSION:   
    if (par->distr->id == UNUR_DISTR_BURR_XI) {
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (gen) GEN->is_inversion = TRUE;
    _unur_cstd_set_sampling_routine(par,gen,_unur_cstd_sample_inv); 
    return UNUR_SUCCESS;
  default: 
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
} 
#undef burr_type
#undef k
#undef c
