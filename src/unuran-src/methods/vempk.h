/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_vempk_new( const UNUR_DISTR *distribution );
int unur_vempk_set_smoothing( UNUR_PAR *parameters, double smoothing );
int unur_vempk_chg_smoothing( UNUR_GEN *generator, double smoothing );
int unur_vempk_set_varcor( UNUR_PAR *parameters, int varcor );
int unur_vempk_chg_varcor( UNUR_GEN *generator, int varcor );
