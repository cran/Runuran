/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_pinv_new( const UNUR_DISTR *distribution );
int unur_pinv_set_order( UNUR_PAR *parameters, int order);
int unur_pinv_set_smoothness( UNUR_PAR *parameters, int smoothness);
int unur_pinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
int unur_pinv_set_extra_testpoints( UNUR_PAR *parameters, int n_points);
int unur_pinv_set_use_upoints( UNUR_PAR *parameters, int use_upoints );
int unur_pinv_set_usepdf( UNUR_PAR *parameters );
int unur_pinv_set_usecdf( UNUR_PAR *parameters );
int unur_pinv_set_boundary( UNUR_PAR *parameters, double left, double right );
int unur_pinv_set_searchboundary( UNUR_PAR *parameters, int left, int right );
int unur_pinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_pinv_get_n_intervals( const UNUR_GEN *generator ); 
int unur_pinv_set_keepcdf( UNUR_PAR *parameters, int keepcdf);
double unur_pinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
double unur_pinv_eval_approxcdf( const UNUR_GEN *generator, double x );
int unur_pinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
