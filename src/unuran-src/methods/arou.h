/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_arou_new( const UNUR_DISTR *distribution );
int unur_arou_set_usedars( UNUR_PAR *parameters, int usedars );
int unur_arou_set_darsfactor( UNUR_PAR *parameters, double factor );
int unur_arou_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
double unur_arou_get_sqhratio( const UNUR_GEN *generator );
double unur_arou_get_hatarea( const UNUR_GEN *generator );
double unur_arou_get_squeezearea( const UNUR_GEN *generator );
int unur_arou_set_max_segments( UNUR_PAR *parameters, int max_segs );
int unur_arou_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
int unur_arou_set_usecenter( UNUR_PAR *parameters, int usecenter );
int unur_arou_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_arou_set_verify( UNUR_PAR *parameters, int verify );
int unur_arou_chg_verify( UNUR_GEN *generator, int verify );
int unur_arou_set_pedantic( UNUR_PAR *parameters, int pedantic );
