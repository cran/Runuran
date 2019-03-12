/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_nrou_new( const UNUR_DISTR *distribution );
int unur_nrou_set_u( UNUR_PAR *parameters, double umin, double umax );
int unur_nrou_set_v( UNUR_PAR *parameters, double vmax );
int unur_nrou_set_r( UNUR_PAR *parameters, double r );
int unur_nrou_set_center( UNUR_PAR *parameters, double center );
int unur_nrou_set_verify( UNUR_PAR *parameters, int verify );
int unur_nrou_chg_verify( UNUR_GEN *generator, int verify );
