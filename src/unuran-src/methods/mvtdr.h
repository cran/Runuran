/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_mvtdr_new( const UNUR_DISTR *distribution );
int unur_mvtdr_set_stepsmin( UNUR_PAR *parameters, int stepsmin );
int unur_mvtdr_set_boundsplitting( UNUR_PAR *parameters, double boundsplitting );
int unur_mvtdr_set_maxcones( UNUR_PAR *parameters, int maxcones );
int unur_mvtdr_get_ncones( const UNUR_GEN *generator );
double unur_mvtdr_get_hatvol( const UNUR_GEN *generator );
int unur_mvtdr_set_verify( UNUR_PAR *parameters, int verify );
int unur_mvtdr_chg_verify( UNUR_GEN *generator, int verify );
