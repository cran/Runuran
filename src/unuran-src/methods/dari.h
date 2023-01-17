/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_dari_new( const UNUR_DISTR *distribution );
int unur_dari_set_squeeze( UNUR_PAR *parameters, int squeeze );
int unur_dari_set_tablesize( UNUR_PAR *parameters, int size );
int unur_dari_set_cpfactor( UNUR_PAR *parameters, double cp_factor );
int unur_dari_set_verify( UNUR_PAR *parameters, int verify );
int unur_dari_chg_verify( UNUR_GEN *generator, int verify );
