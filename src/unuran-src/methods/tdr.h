/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_tdr_new( const UNUR_DISTR* distribution );
int unur_tdr_set_c( UNUR_PAR *parameters, double c );
int unur_tdr_set_variant_gw( UNUR_PAR *parameters );
int unur_tdr_set_variant_ps( UNUR_PAR *parameters );
int unur_tdr_set_variant_ia( UNUR_PAR *parameters );
int unur_tdr_set_usedars( UNUR_PAR *parameters, int usedars );
int unur_tdr_set_darsfactor( UNUR_PAR *parameters, double factor );
int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
int unur_tdr_set_reinit_percentiles( UNUR_PAR *parameters, int n_percentiles, const double *percentiles );
int unur_tdr_chg_reinit_percentiles( UNUR_GEN *generator, int n_percentiles, const double *percentiles );
int unur_tdr_set_reinit_ncpoints( UNUR_PAR *parameters, int ncpoints );
int unur_tdr_chg_reinit_ncpoints( UNUR_GEN *generator, int ncpoints );
int unur_tdr_chg_truncated(UNUR_GEN *gen, double left, double right);
int unur_tdr_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
double unur_tdr_get_sqhratio( const UNUR_GEN *generator );
double unur_tdr_get_hatarea( const UNUR_GEN *generator );
double unur_tdr_get_squeezearea( const UNUR_GEN *generator );
int unur_tdr_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_tdr_set_verify( UNUR_PAR *parameters, int verify );
int unur_tdr_chg_verify( UNUR_GEN *generator, int verify );
int unur_tdr_set_pedantic( UNUR_PAR *parameters, int pedantic );
double unur_tdr_eval_invcdfhat( const UNUR_GEN *generator, double u, 
				double *hx, double *fx, double *sqx );
int _unur_tdr_is_ARS_running( const UNUR_GEN *generator );
