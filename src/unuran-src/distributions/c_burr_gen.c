/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
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
    _unur_cstd_set_sampling_routine(par,gen,_unur_stdgen_sample_burr_inv); 
    return UNUR_SUCCESS;
  default: 
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
} 
double _unur_stdgen_sample_burr_inv( struct unur_gen *gen )
{
  double U, Y;
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);
  while (_unur_iszero(U = GEN->umin + uniform() * (GEN->umax - GEN->umin)));
  switch (burr_type) {
  case UNUR_DISTR_BURR_I:
    return U;
  case UNUR_DISTR_BURR_II:
    Y = exp( -log(U)/k );  
    return ( -log( Y - 1. ) );
  case UNUR_DISTR_BURR_III:
    Y = exp( -log(U)/k );  
    return ( exp( -log( Y - 1. )/c ) );
  case UNUR_DISTR_BURR_IV:
    Y = exp( -log(U)/k );   
    Y = exp( c * log( Y - 1. )) + 1.;
    return (c/Y);
  case UNUR_DISTR_BURR_V:
    Y = exp( -log(U)/k );   
    return atan( -log( (Y - 1.) / c ) );
  case UNUR_DISTR_BURR_VI:
    Y = exp( -log(U)/k );   
    Y = -log( (Y - 1.) / c)/k;
    return log( Y + sqrt(Y * Y +1.));
  case UNUR_DISTR_BURR_VII:
    Y = exp( log(U)/k );    
    return ( log(2. * Y / (2. - 2.*Y)) / 2. );
  case UNUR_DISTR_BURR_VIII:
    Y = exp( log(U)/k );    
    return ( log( tan( Y * M_PI/2. ) ) );
  case UNUR_DISTR_BURR_IX:
    Y = 1. + 2. * U / (c * (1.-U));
    return log( exp( log(Y) / k) - 1. );
  case UNUR_DISTR_BURR_X:
  Y = exp( log(U)/k );   
    return ( sqrt( -log( 1. - Y ) ) );
  case UNUR_DISTR_BURR_XII:
    Y = exp( -log(U)/k );   
    return ( exp( log( Y - 1.) / c) );
  case UNUR_DISTR_BURR_XI:
  default:
    _unur_error(NULL,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }
} 
#undef burr_type
#undef k
#undef c
