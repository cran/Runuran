/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_H_SEEN
#define URNG_H_SEEN
UNUR_URNG *unur_get_default_urng( void );
UNUR_URNG *unur_set_default_urng( UNUR_URNG *urng_new );
UNUR_URNG *unur_set_default_urng_aux( UNUR_URNG *urng_new );
UNUR_URNG *unur_get_default_urng_aux( void );
int unur_set_urng( UNUR_PAR *parameters, UNUR_URNG *urng );
UNUR_URNG *unur_chg_urng( UNUR_GEN *generator, UNUR_URNG *urng );
UNUR_URNG *unur_get_urng( UNUR_GEN *generator );
int unur_set_urng_aux( UNUR_PAR *parameters, UNUR_URNG *urng_aux );
int unur_use_urng_aux_default( UNUR_PAR *parameters );
int unur_chgto_urng_aux_default( UNUR_GEN *generator );
UNUR_URNG *unur_chg_urng_aux( UNUR_GEN *generator, UNUR_URNG *urng_aux );
UNUR_URNG *unur_get_urng_aux( UNUR_GEN *generator );
double unur_urng_sample (UNUR_URNG *urng);
double unur_sample_urng (UNUR_GEN *gen);
int unur_urng_sample_array (UNUR_URNG *urng, double *X, int dim);
int unur_urng_reset (UNUR_URNG *urng);
#ifdef UNUR_URNG_UNURAN
int unur_urng_sync (UNUR_URNG *urng);
int unur_urng_seed (UNUR_URNG *urng, unsigned long seed);
int unur_urng_anti (UNUR_URNG *urng, int anti);
int unur_urng_nextsub (UNUR_URNG *urng);
int unur_urng_resetsub (UNUR_URNG *urng);
int unur_gen_sync (UNUR_GEN *generator);
int unur_gen_seed (UNUR_GEN *generator, unsigned long seed);
int unur_gen_anti (UNUR_GEN *generator, int anti);
int unur_gen_reset (UNUR_GEN *generator);
int unur_gen_nextsub (UNUR_GEN *generator);
int unur_gen_resetsub (UNUR_GEN *generator);
UNUR_URNG *unur_urng_new( double (*sampleunif)(void *state), void *state );
void unur_urng_free (UNUR_URNG *urng);
int unur_urng_set_sample_array( UNUR_URNG *urng, unsigned int (*samplearray)(void *state, double *X, int dim) );
int unur_urng_set_sync( UNUR_URNG *urng, void (*sync)(void *state) );
int unur_urng_set_seed( UNUR_URNG *urng, void (*setseed)(void *state, unsigned long seed) );
int unur_urng_set_anti( UNUR_URNG *urng, void (*setanti)(void *state, int anti) );
int unur_urng_set_reset( UNUR_URNG *urng, void (*reset)(void *state) );
int unur_urng_set_nextsub( UNUR_URNG *urng, void (*nextsub)(void *state) );
int unur_urng_set_resetsub( UNUR_URNG *urng, void (*resetsub)(void *state) );
int unur_urng_set_delete( UNUR_URNG *urng, void (*fpdelete)(void *state) );
#endif   
#endif  
