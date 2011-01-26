/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_hitro_new( const UNUR_DISTR *distribution );
int unur_hitro_set_variant_coordinate( UNUR_PAR *parameters );
int unur_hitro_set_variant_random_direction( UNUR_PAR *parameters );
int unur_hitro_set_use_adaptiveline( UNUR_PAR *parameters, int adaptive );
int unur_hitro_set_use_boundingrectangle( UNUR_PAR *parameters, int rectangle );
int unur_hitro_set_use_adaptiverectangle( UNUR_PAR *parameters, int adaptive );
int unur_hitro_set_r( UNUR_PAR *parameters, double r );
int unur_hitro_set_v( UNUR_PAR *parameters, double vmax );
int unur_hitro_set_u( UNUR_PAR *parameters, const double *umin, const double *umax );
int unur_hitro_set_adaptive_multiplier( UNUR_PAR *parameters, double factor );
int unur_hitro_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
int unur_hitro_set_thinning( UNUR_PAR *parameters, int thinning );
int unur_hitro_set_burnin( UNUR_PAR *parameters, int burnin );
const double *unur_hitro_get_state( UNUR_GEN *generator );
int unur_hitro_chg_state( UNUR_GEN *generator, const double *state );
int unur_hitro_reset_state( UNUR_GEN *generator );
