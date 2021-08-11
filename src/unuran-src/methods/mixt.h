/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_mixt_new( int n, const double *prob, UNUR_GEN **comp );
int unur_mixt_set_useinversion( UNUR_PAR *parameters, int useinv );
double unur_mixt_eval_invcdf( const UNUR_GEN *generator, double u );
