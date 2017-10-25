/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_GSLQRNG_H_SEEN
#define URNG_GSLQRNG_H_SEEN
#include <gsl/gsl_qrng.h>
UNUR_URNG *unur_urng_gslqrng_new( const gsl_qrng_type *qrngtype, unsigned int dim );
#endif  
