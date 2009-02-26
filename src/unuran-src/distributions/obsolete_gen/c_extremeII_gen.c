/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
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
#define k      (DISTR.params[0])    
#define zeta   (DISTR.params[1])    
#define theta  (DISTR.params[2])    
int 
_unur_stdgen_extremeII_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case UNUR_STDGEN_INVERSION:   
    if (gen) GEN->is_inversion = TRUE;
    _unur_cstd_set_sampling_routine(par,gen,_unur_cstd_sample_inv); 
    return UNUR_SUCCESS;
  default: 
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
} 
