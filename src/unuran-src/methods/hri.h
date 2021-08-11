/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_hri_new( const UNUR_DISTR *distribution );
int unur_hri_set_p0( UNUR_PAR *parameters, double p0 );
int unur_hri_set_verify( UNUR_PAR *parameters, int verify );
int unur_hri_chg_verify( UNUR_GEN *generator, int verify );
