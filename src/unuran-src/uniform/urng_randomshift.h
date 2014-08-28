/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_RANDOMSHIFT_H_SEEN
#define URNG_RANDOMSHIFT_H_SEEN
UNUR_URNG *unur_urng_randomshift_new( UNUR_URNG *qrng, UNUR_URNG *srng, int dim );
int unur_urng_randomshift_nextshift( UNUR_URNG *urng );
#endif  
