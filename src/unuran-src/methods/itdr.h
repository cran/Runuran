/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_itdr_new( const UNUR_DISTR *distribution );
int unur_itdr_set_xi( UNUR_PAR *parameters, double xi );
int unur_itdr_set_cp( UNUR_PAR *parameters, double cp );
int unur_itdr_set_ct( UNUR_PAR *parameters, double ct );
double unur_itdr_get_xi( UNUR_GEN *generator );
double unur_itdr_get_cp( UNUR_GEN *generator );
double unur_itdr_get_ct( UNUR_GEN *generator );
double unur_itdr_get_area( UNUR_GEN *generator );
int unur_itdr_set_verify( UNUR_PAR *parameters, int verify );
int unur_itdr_chg_verify( UNUR_GEN *generator, int verify );
