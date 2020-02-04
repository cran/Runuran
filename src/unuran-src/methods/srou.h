/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_srou_new( const UNUR_DISTR *distribution );
int unur_srou_set_r( UNUR_PAR *parameters, double r );
int unur_srou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
int unur_srou_set_pdfatmode( UNUR_PAR *parameters, double fmode );
int unur_srou_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
int unur_srou_set_usemirror( UNUR_PAR *parameters, int usemirror );
int unur_srou_set_verify( UNUR_PAR *parameters, int verify );
int unur_srou_chg_verify( UNUR_GEN *generator, int verify );
int unur_srou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
int unur_srou_chg_pdfatmode( UNUR_GEN *generator, double fmode );
