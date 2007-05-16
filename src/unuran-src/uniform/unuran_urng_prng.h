/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNURAN_URNG_PRNG_H_SEEN
#define UNURAN_URNG_PRNG_H_SEEN
#include <prng.h>
UNUR_URNG *unur_urng_prng_new( const char *prngstr );
UNUR_URNG *unur_urng_prngptr_new( struct prng *urng );
#endif  
