/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_gibbs_new( const UNUR_DISTR *distribution );
int unur_gibbs_set_variant_coordinate( UNUR_PAR *parameters );
int unur_gibbs_set_variant_random_direction( UNUR_PAR *parameters );
int unur_gibbs_set_c( UNUR_PAR *parameters, double c );
int unur_gibbs_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
int unur_gibbs_set_thinning( UNUR_PAR *parameters, int thinning );
int unur_gibbs_set_burnin( UNUR_PAR *parameters, int burnin );
const double *unur_gibbs_get_state( UNUR_GEN *generator );
int unur_gibbs_chg_state( UNUR_GEN *generator, const double *state );
int unur_gibbs_reset_state( UNUR_GEN *generator );
