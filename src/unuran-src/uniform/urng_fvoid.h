/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_FVOID_H_SEEN
#define URNG_FVOID_H_SEEN
UNUR_URNG *unur_urng_fvoid_new( double (*urand)(void *state), void (*reset)(void *state) );
#endif  
