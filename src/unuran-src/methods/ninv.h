/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_ninv_new( const UNUR_DISTR *distribution );
int unur_ninv_set_useregula( UNUR_PAR *parameters );
int unur_ninv_set_usenewton( UNUR_PAR *parameters );
int unur_ninv_set_usebisect( UNUR_PAR *parameters );
int unur_ninv_set_max_iter( UNUR_PAR *parameters, int max_iter );
int unur_ninv_chg_max_iter(UNUR_GEN *generator, int max_iter);
int unur_ninv_set_x_resolution( UNUR_PAR *parameters, double x_resolution);
int unur_ninv_chg_x_resolution(UNUR_GEN *generator, double x_resolution);
int unur_ninv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
int unur_ninv_chg_u_resolution(UNUR_GEN *generator, double u_resolution);
int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
int unur_ninv_chg_start(UNUR_GEN *gen, double left, double right);
int unur_ninv_set_table(UNUR_PAR *parameters, int no_of_points);
int unur_ninv_chg_table(UNUR_GEN *gen, int no_of_points);
int unur_ninv_chg_truncated(UNUR_GEN *gen, double left, double right);
double unur_ninv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
