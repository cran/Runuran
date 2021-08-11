/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_DEPRECATED_METHODS_H_SEEN
#define UNUR_DEPRECATED_METHODS_H_SEEN
int unur_cstd_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_dari_reinit( UNUR_GEN *generator );
int unur_dari_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_dari_chg_domain( UNUR_GEN *generator, int left, int right );
int unur_dari_chg_mode( UNUR_GEN *generator, int mode );
int unur_dari_upd_mode( UNUR_GEN *generator );
int unur_dari_chg_pmfsum( UNUR_GEN *generator, double sum );
int unur_dari_upd_pmfsum( UNUR_GEN *generator );
int unur_dsrou_reinit( UNUR_GEN *generator );
int unur_dsrou_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_dsrou_chg_domain( UNUR_GEN *generator, int left, int right );
int unur_dsrou_chg_mode( UNUR_GEN *generator, int mode );
int unur_dsrou_upd_mode( UNUR_GEN *generator );
int unur_dsrou_chg_pmfsum( UNUR_GEN *generator, double sum );
int unur_dsrou_upd_pmfsum( UNUR_GEN *generator );
int unur_dstd_chg_pmfparams( UNUR_GEN *gen, double *params, int n_params );
int unur_ninv_chg_pdfparams(UNUR_GEN *generator, double *params, int n_params);
int unur_srou_reinit( UNUR_GEN *generator );
int unur_srou_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
int unur_srou_upd_mode( UNUR_GEN *generator );
int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
int unur_srou_upd_pdfarea( UNUR_GEN *generator );
int unur_ssr_reinit( UNUR_GEN *generator );
int unur_ssr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_ssr_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_ssr_chg_mode( UNUR_GEN *generator, double mode );
int unur_ssr_upd_mode( UNUR_GEN *generator );
int unur_ssr_chg_pdfarea( UNUR_GEN *generator, double area );
int unur_ssr_upd_pdfarea( UNUR_GEN *generator );
int unur_tdr_reinit( UNUR_GEN *generator );
int unur_tdrgw_reinit( UNUR_GEN *generator );
int unur_utdr_reinit( UNUR_GEN *generator );
int unur_utdr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_utdr_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_utdr_chg_mode( UNUR_GEN *generator, double mode );
int unur_utdr_upd_mode( UNUR_GEN *generator );
int unur_utdr_chg_pdfarea( UNUR_GEN *generator, double area );
int unur_utdr_upd_pdfarea( UNUR_GEN *generator );
#endif  
