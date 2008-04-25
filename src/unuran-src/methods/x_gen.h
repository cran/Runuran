/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_GEN *unur_init( UNUR_PAR *parameters );
int unur_reinit( UNUR_GEN *generator );
int    unur_sample_discr(UNUR_GEN *generator);
double unur_sample_cont(UNUR_GEN *generator);
int    unur_sample_vec(UNUR_GEN *generator, double *vector);
int    unur_sample_matr(UNUR_GEN *generator, double *matrix);
void  unur_free( UNUR_GEN *generator );
const char *unur_gen_info( UNUR_GEN *generator, int help );
int unur_get_dimension( const UNUR_GEN *generator );
const char *unur_get_genid( const UNUR_GEN *generator );
UNUR_DISTR *unur_get_distr( const UNUR_GEN *generator );
int unur_set_use_distr_privatecopy( UNUR_PAR *parameters, int use_privatecopy );
UNUR_GEN *unur_gen_clone( const UNUR_GEN *gen );
void unur_par_free( UNUR_PAR *par);
