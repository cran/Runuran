/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_cvemp_new( int dim ); 
int unur_distr_cvemp_set_data( UNUR_DISTR *distribution, const double *sample, int n_sample );
int unur_distr_cvemp_read_data( UNUR_DISTR *distribution, const char *filename );
int unur_distr_cvemp_get_data( const UNUR_DISTR *distribution, const double **sample );
