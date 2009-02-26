/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_GSL_H_SEEN
#define URNG_GSL_H_SEEN
#include <gsl/gsl_rng.h>
UNUR_URNG *unur_urng_gsl_new( const gsl_rng_type *urngtype );
UNUR_URNG *unur_urng_gslptr_new( gsl_rng *urng );
#endif  
