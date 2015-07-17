/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_vnrou_new( const UNUR_DISTR *distribution );
int unur_vnrou_set_u( UNUR_PAR *parameters, double *umin, double *umax );
int unur_vnrou_chg_u( UNUR_GEN *generator, double *umin, double *umax );
int unur_vnrou_set_v( UNUR_PAR *parameters, double vmax );
int unur_vnrou_chg_v( UNUR_GEN *generator, double vmax );
int unur_vnrou_set_r( UNUR_PAR *parameters, double r );
int unur_vnrou_set_verify( UNUR_PAR *parameters, int verify );
int unur_vnrou_chg_verify( UNUR_GEN *generator, int verify );
double unur_vnrou_get_volumehat( const UNUR_GEN *generator );
