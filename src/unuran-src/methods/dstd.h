/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_dstd_new( const UNUR_DISTR *distribution );
int unur_dstd_set_variant( UNUR_PAR *parameters, unsigned variant );
int unur_dstd_chg_truncated( UNUR_GEN *generator, int left, int right );
int unur_dstd_eval_invcdf( const UNUR_GEN *generator, double u );
