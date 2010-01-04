/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr.h"
#include "distr_source.h"
struct unur_distr *
_unur_distr_generic_new( void )
{
  register struct unur_distr *distr;
  distr = _unur_xmalloc( sizeof(struct unur_distr) );
  if (!distr) return NULL;
  distr->type = UNUR_DISTR_GENERIC;
  distr->id = UNUR_DISTR_GENERIC;
  distr->dim = 1;   
  distr->name = "unknown";
  distr->name_str = NULL;
  distr->base = NULL;
  distr->destroy = NULL;     
  distr->clone   = NULL;     
  distr->extobj  = NULL;     
  distr->set     = 0u;       
  return distr;
} 
void
unur_distr_free( struct unur_distr *distr )
{
  if (distr) _unur_distr_free( distr );
} 
int 
unur_distr_set_name( struct unur_distr *distr, const char *name )
{
  size_t len;
  char *name_str;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  len = strlen(name) + 1;
  name_str = _unur_xrealloc(distr->name_str,len);
  memcpy( name_str, name, len );
  distr->name_str = name_str;
  distr->name = name_str;
  return UNUR_SUCCESS;
} 
const char *
unur_distr_get_name( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  return distr->name;
} 
int
unur_distr_get_dim( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, 0 );
  return distr->dim;
} 
unsigned int 
unur_distr_get_type( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, 0u );
  return (distr->type);
} 
int 
unur_distr_is_cont( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  return ((distr->type == UNUR_DISTR_CONT) ? TRUE : FALSE);
} 
int 
unur_distr_is_cvec( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  return ((distr->type == UNUR_DISTR_CVEC) ? TRUE : FALSE);
} 
int 
unur_distr_is_cvemp( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  return ((distr->type == UNUR_DISTR_CVEMP) ? TRUE : FALSE);
} 
int 
unur_distr_is_matr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  return ((distr->type == UNUR_DISTR_MATR) ? TRUE : FALSE);
} 
int 
unur_distr_is_discr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  return ((distr->type == UNUR_DISTR_DISCR) ? TRUE : FALSE);
} 
int 
unur_distr_is_cemp( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, FALSE );
  return ((distr->type == UNUR_DISTR_CEMP) ? TRUE : FALSE);
} 
struct unur_distr *
unur_distr_clone( const struct unur_distr *distr )
{
  _unur_check_NULL( "Clone", distr, NULL );
  _unur_check_NULL( "Clone", distr->clone, NULL );
  return (distr->clone(distr));
} 
int
unur_distr_set_extobj( struct unur_distr *distr, const void *extobj )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  distr->extobj = extobj;
  return UNUR_SUCCESS;
} 
const void *
unur_distr_get_extobj( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  return distr->extobj;
} 
