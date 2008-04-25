/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_mcorr_new( const UNUR_DISTR *distribution );
int unur_mcorr_set_eigenvalues( UNUR_PAR *par, const double *eigenvalues );
int unur_mcorr_chg_eigenvalues( UNUR_GEN *gen, const double *eigenvalues );
