/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <utils/matrix_source.h>
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "matr.h"
#define DISTR distr->data.matr
static void _unur_distr_matr_free( struct unur_distr *distr );
struct unur_distr *
unur_distr_matr_new( int n_rows, int n_cols )
{
  register struct unur_distr *distr;
  if (n_rows < 1 || n_cols < 1) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"n_rows or n_cols < 1");
    return NULL;
  }
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;
  COOKIE_SET(distr,CK_DISTR_MATR);
  distr->type = UNUR_DISTR_MATR;
  distr->id = UNUR_DISTR_GENERIC;
  DISTR.n_rows = n_rows;
  DISTR.n_cols = n_cols;
  distr->dim = n_rows * n_cols;
  distr->destroy = _unur_distr_matr_free;
  distr->clone = _unur_distr_matr_clone;
  DISTR.init      = NULL;   
  return distr;
} 
struct unur_distr *
_unur_distr_matr_clone( const struct unur_distr *distr )
{
#define CLONE clone->data.matr
  struct unur_distr *clone;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, MATR, NULL );
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  memcpy( clone, distr, sizeof( struct unur_distr ) );
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }
  return clone;
#undef CLONE
} 
void
_unur_distr_matr_free( struct unur_distr *distr )
{
  if( distr == NULL ) 
    return;
  COOKIE_CHECK(distr,CK_DISTR_MATR,RETURN_VOID);
  if (distr->name_str) free(distr->name_str);
  COOKIE_CLEAR(distr);
  free( distr );
} 
int
unur_distr_matr_get_dim( const struct unur_distr *distr, int *n_rows, int *n_cols )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, MATR, 0 );
  CHECK_NULL( n_rows, 0 );
  CHECK_NULL( n_cols, 0 );
  *n_rows = DISTR.n_rows;
  *n_cols = DISTR.n_cols;
  return distr->dim;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_matr_debug( const struct unur_distr *distr, const char *genid )
{
  FILE *LOG;
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_MATR,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = matrix distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);
  fprintf(LOG,"%s:\tdimension = %d x %d   (= %d)\n",genid,DISTR.n_rows,DISTR.n_cols,distr->dim);
  fprintf(LOG,"%s:\n",genid);
} 
#endif    
#undef DISTR
