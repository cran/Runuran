/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_tabl_new( const UNUR_DISTR* distribution );
int unur_tabl_set_variant_ia( UNUR_PAR *parameters, int use_ia );
int unur_tabl_set_cpoints( UNUR_PAR *parameters, int n_cpoints, const double *cpoints );
int unur_tabl_set_nstp( UNUR_PAR *parameters, int n_stp );
int unur_tabl_set_useear( UNUR_PAR *parameters, int useear );
int unur_tabl_set_areafraction( UNUR_PAR *parameters, double fraction );
int unur_tabl_set_usedars( UNUR_PAR *parameters, int usedars );
int unur_tabl_set_darsfactor( UNUR_PAR *parameters, double factor );
int unur_tabl_set_variant_splitmode( UNUR_PAR *parameters, unsigned splitmode );
int unur_tabl_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
double unur_tabl_get_sqhratio( const UNUR_GEN *generator );
double unur_tabl_get_hatarea( const UNUR_GEN *generator );
double unur_tabl_get_squeezearea( const UNUR_GEN *generator );
int unur_tabl_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_tabl_get_n_intervals( const UNUR_GEN *generator );
int unur_tabl_set_slopes( UNUR_PAR *parameters, const double *slopes, int n_slopes );
int unur_tabl_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_tabl_set_boundary( UNUR_PAR *parameters, double left, double right );
int unur_tabl_chg_truncated(UNUR_GEN *gen, double left, double right);
int unur_tabl_set_verify( UNUR_PAR *parameters, int verify );
int unur_tabl_chg_verify( UNUR_GEN *generator, int verify );
int unur_tabl_set_pedantic( UNUR_PAR *parameters, int pedantic );
