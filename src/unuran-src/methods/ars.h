/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_ars_new( const UNUR_DISTR* distribution );
int unur_ars_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_ars_set_cpoints( UNUR_PAR *parameters, int n_cpoints, const double *cpoints );
int unur_ars_set_reinit_percentiles( UNUR_PAR *parameters, int n_percentiles, const double *percentiles );
int unur_ars_chg_reinit_percentiles( UNUR_GEN *generator, int n_percentiles, const double *percentiles );
int unur_ars_set_reinit_ncpoints( UNUR_PAR *parameters, int ncpoints );
int unur_ars_chg_reinit_ncpoints( UNUR_GEN *generator, int ncpoints );
int unur_ars_set_max_iter( UNUR_PAR *parameters, int max_iter );
int unur_ars_set_verify( UNUR_PAR *parameters, int verify );
int unur_ars_chg_verify( UNUR_GEN *generator, int verify );
int unur_ars_set_pedantic( UNUR_PAR *parameters, int pedantic );
double unur_ars_get_loghatarea( const UNUR_GEN *generator );
double unur_ars_eval_invcdfhat( const UNUR_GEN *generator, double u );
