/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_condi_new( const UNUR_DISTR *distribution, const double *pos, const double *dir, int k );
int unur_distr_condi_set_condition( struct unur_distr *distribution, const double *pos, const double *dir, int k );
int unur_distr_condi_get_condition( struct unur_distr *distribution, const double **pos, const double **dir, int *k );
const UNUR_DISTR *unur_distr_condi_get_distribution( const UNUR_DISTR *distribution );
