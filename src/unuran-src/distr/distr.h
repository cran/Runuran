/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

enum {
  UNUR_DISTR_CONT  = 0x010u,      
  UNUR_DISTR_CEMP  = 0x011u,      
  UNUR_DISTR_CVEC  = 0x110u,      
  UNUR_DISTR_CVEMP = 0x111u,      
  UNUR_DISTR_MATR  = 0x210u,      
  UNUR_DISTR_DISCR = 0x020u       
};
void unur_distr_free( UNUR_DISTR *distribution );
int unur_distr_set_name( UNUR_DISTR *distribution, const char *name );
const char *unur_distr_get_name( const UNUR_DISTR *distribution );
int unur_distr_get_dim( const UNUR_DISTR *distribution );
unsigned int unur_distr_get_type( const UNUR_DISTR *distribution );
int unur_distr_is_cont( const UNUR_DISTR *distribution );
int unur_distr_is_cvec( const UNUR_DISTR *distribution );
int unur_distr_is_cemp( const UNUR_DISTR *distribution );
int unur_distr_is_cvemp( const UNUR_DISTR *distribution );
int unur_distr_is_discr( const UNUR_DISTR *distribution );
int unur_distr_is_matr( const UNUR_DISTR *distribution );
int unur_distr_set_extobj( UNUR_DISTR *distribution, const void *extobj );
const void *unur_distr_get_extobj( const UNUR_DISTR *distribution );
UNUR_DISTR *unur_distr_clone( const UNUR_DISTR *distr );
