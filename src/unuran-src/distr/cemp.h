/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_cemp_new( void );
int unur_distr_cemp_set_data( UNUR_DISTR *distribution, const double *sample, int n_sample );
int unur_distr_cemp_read_data( UNUR_DISTR *distribution, const char *filename );
int unur_distr_cemp_get_data( const UNUR_DISTR *distribution, const double **sample );
int unur_distr_cemp_set_hist( UNUR_DISTR *distribution, const double *prob, int n_prob, double xmin, double xmax );
int unur_distr_cemp_set_hist_prob( UNUR_DISTR *distribution, const double *prob, int n_prob );
int unur_distr_cemp_set_hist_domain( UNUR_DISTR *distribution, double xmin, double xmax );
int unur_distr_cemp_set_hist_bins( UNUR_DISTR *distribution, const double *bins, int n_bins );
