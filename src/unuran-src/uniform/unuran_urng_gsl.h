/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNURAN_URNG_GSL_H_SEEN
#define UNURAN_URNG_GSL_H_SEEN
#include <gsl/gsl_rng.h>
UNUR_URNG *unur_urng_gsl_new( const gsl_rng_type *urngtype );
UNUR_URNG *unur_urng_gslptr_new( gsl_rng *urng );
#include <gsl/gsl_qrng.h>
UNUR_URNG *unur_urng_gslqrng_new( const gsl_qrng_type *qrngtype, unsigned int dim );
#endif  
