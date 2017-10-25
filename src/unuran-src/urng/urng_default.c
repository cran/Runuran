/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "urng.h"
#include <uniform/urng_builtin.h>
#include <uniform/urng_fvoid.h>
#include <uniform/urng_randomshift.h>
#if defined(UNURAN_HAS_GSL) && defined(UNUR_URNG_UNURAN)
#  include <uniform/urng_gsl.h>
#  include <uniform/urng_gslqrng.h>
#endif
#if defined(UNURAN_HAS_PRNG) && defined(UNUR_URNG_UNURAN)
#  include <uniform/urng_prng.h>
#endif
#if defined(UNURAN_HAS_RNGSTREAM) && defined(UNUR_URNG_UNURAN)
#  include <uniform/urng_rngstreams.h>
#endif
static UNUR_URNG *urng_default = NULL;
static UNUR_URNG *urng_aux_default = NULL;
UNUR_URNG *
unur_get_default_urng( void )
{
  if( urng_default == NULL ) {
    urng_default = UNUR_URNG_DEFAULT;
    if( urng_default == NULL ) {
      _unur_error("URNG",UNUR_ERR_NULL,"Cannot set default URNG. EXIT !!!");
#ifndef R_UNURAN
      exit(EXIT_FAILURE);
#endif
    }
  }
  return (urng_default);
} 
UNUR_URNG *
unur_set_default_urng( UNUR_URNG *urng_new )
{
  UNUR_URNG *urng_old = urng_default;
  _unur_check_NULL("URNG", urng_new, urng_default);
  urng_default = urng_new;     
  return (urng_old);
} 
UNUR_URNG *
unur_get_default_urng_aux( void )
{
  if( urng_aux_default == NULL ) {
    urng_aux_default = UNUR_URNG_AUX_DEFAULT;
    if( urng_aux_default == NULL ) {
      _unur_error("URNG",UNUR_ERR_NULL,"Cannot set default auxilliary URNG. EXIT !!!");
#ifndef R_UNURAN
      exit(EXIT_FAILURE);
#endif
    }
  }
  return (urng_aux_default);
} 
UNUR_URNG *
unur_set_default_urng_aux( UNUR_URNG *urng_aux_new )
{
  UNUR_URNG *urng_aux_old = urng_aux_default;
  _unur_check_NULL("URNG", urng_aux_new, urng_aux_default);
  urng_aux_default = urng_aux_new;     
  return (urng_aux_old);
} 
