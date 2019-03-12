/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_hinv_new( const UNUR_DISTR *distribution );
int unur_hinv_set_order( UNUR_PAR *parameters, int order);
int unur_hinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
int unur_hinv_set_cpoints( UNUR_PAR *parameters, const double *stp, int n_stp );
int unur_hinv_set_boundary( UNUR_PAR *parameters, double left, double right );
int unur_hinv_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_hinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_hinv_get_n_intervals( const UNUR_GEN *generator );
double unur_hinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
int unur_hinv_chg_truncated( UNUR_GEN *generator, double left, double right );
int unur_hinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
