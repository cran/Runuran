/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_PAR *unur_dext_new( const UNUR_DISTR *distribution );
int unur_dext_set_init( UNUR_PAR *parameters, int (*init)(UNUR_GEN *gen) );
int unur_dext_set_sample( UNUR_PAR *parameters, int (*sample)(UNUR_GEN *gen) );
void *unur_dext_get_params( UNUR_GEN *generator, size_t size );
double *unur_dext_get_distrparams( UNUR_GEN *generator );
int unur_dext_get_ndistrparams( UNUR_GEN *generator );
