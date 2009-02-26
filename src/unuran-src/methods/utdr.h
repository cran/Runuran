/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_utdr_new( const UNUR_DISTR *distribution );
int unur_utdr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
int unur_utdr_set_cpfactor( UNUR_PAR *parameters, double cp_factor );
int unur_utdr_set_deltafactor( UNUR_PAR *parameters, double delta );
int unur_utdr_set_verify( UNUR_PAR *parameters, int verify );
int unur_utdr_chg_verify( UNUR_GEN *generator, int verify );
int unur_utdr_chg_pdfatmode( UNUR_GEN *generator, double fmode );
