/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_dsrou_new( const UNUR_DISTR *distribution );
int unur_dsrou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
int unur_dsrou_set_verify( UNUR_PAR *parameters, int verify );
int unur_dsrou_chg_verify( UNUR_GEN *generator, int verify );
int unur_dsrou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
