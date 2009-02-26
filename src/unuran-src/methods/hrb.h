/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_hrb_new( const UNUR_DISTR *distribution );
int unur_hrb_set_upperbound( UNUR_PAR *parameters, double upperbound );
int unur_hrb_set_verify( UNUR_PAR *parameters, int verify );
int unur_hrb_chg_verify( UNUR_GEN *generator, int verify );
