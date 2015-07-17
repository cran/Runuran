/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_dgt_new( const UNUR_DISTR *distribution );
int unur_dgt_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_dgt_set_variant( UNUR_PAR *parameters, unsigned variant );
int unur_dgt_eval_invcdf_recycle( const UNUR_GEN *generator, double u, double *recycle );
int unur_dgt_eval_invcdf( const UNUR_GEN *generator, double u );
