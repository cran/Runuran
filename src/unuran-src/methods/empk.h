/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_empk_new( const UNUR_DISTR *distribution );
int unur_empk_set_kernel( UNUR_PAR *parameters, unsigned kernel);
int unur_empk_set_kernelgen( UNUR_PAR *parameters, const UNUR_GEN *kernelgen, double alpha, double kernelvar );
int unur_empk_set_beta( UNUR_PAR *parameters, double beta );
int unur_empk_set_smoothing( UNUR_PAR *parameters, double smoothing );
int unur_empk_chg_smoothing( UNUR_GEN *generator, double smoothing );
int unur_empk_set_varcor( UNUR_PAR *parameters, int varcor );
int unur_empk_chg_varcor( UNUR_GEN *generator, int varcor );
int unur_empk_set_positive( UNUR_PAR *parameters, int positive );
