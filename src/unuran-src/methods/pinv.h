/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_pinv_new( const UNUR_DISTR *distribution );
int unur_pinv_set_order( UNUR_PAR *parameters, int order);
int unur_pinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
int unur_pinv_set_usepdf( UNUR_PAR *parameters );
int unur_pinv_set_usecdf( UNUR_PAR *parameters );
int unur_pinv_set_boundary( UNUR_PAR *parameters, double left, double right );
int unur_pinv_set_searchboundary( UNUR_PAR *parameters, int left, int right );
int unur_pinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_pinv_get_n_intervals( const UNUR_GEN *generator ); 
double unur_pinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
int unur_pinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
