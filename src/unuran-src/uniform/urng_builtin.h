/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_BUILTIN_H_SEEN
#define URNG_BUILTIN_H_SEEN
double unur_urng_MRG31k3p (void *dummy);
void unur_urng_MRG31k3p_seed (void *dummy, unsigned long seed);
void unur_urng_MRG31k3p_reset (void *dummy);
double unur_urng_fish (void *dummy);
void unur_urng_fish_seed (void *dummy, unsigned long seed);
void unur_urng_fish_reset (void *dummy);
double unur_urng_mstd (void *dummy);
void unur_urng_mstd_seed (void *dummy, unsigned long seed);
void unur_urng_mstd_reset (void *dummy);
UNUR_URNG *unur_urng_builtin( void );
UNUR_URNG *unur_urng_builtin_aux( void );
#endif  
