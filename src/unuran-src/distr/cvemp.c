/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "cvemp.h"
#define DISTR distr->data.cvemp
static void _unur_distr_cvemp_free( struct unur_distr *distr );
struct unur_distr *
unur_distr_cvemp_new( int dim )
{
  register struct unur_distr *distr;
  if (dim < 2) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"dimension < 2");
    return NULL;
  }
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;
  COOKIE_SET(distr,CK_DISTR_CVEMP);
  distr->type = UNUR_DISTR_CVEMP;
  distr->id = UNUR_DISTR_GENERIC;
  distr->dim = dim;   
  distr->name = "(empirical)";
  distr->name_str = NULL;
  distr->destroy = _unur_distr_cvemp_free;
  distr->clone = _unur_distr_cvemp_clone;
  DISTR.sample    = NULL;    
  DISTR.n_sample  = 0;       
  return distr;
} 
struct unur_distr *
_unur_distr_cvemp_clone( const struct unur_distr *distr )
{
#define CLONE clone->data.cvemp
  struct unur_distr *clone;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEMP, NULL );
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  memcpy( clone, distr, sizeof( struct unur_distr ) );
  if (DISTR.sample) {
    CLONE.sample = _unur_xmalloc( DISTR.n_sample * distr->dim * sizeof(double) );
    memcpy( CLONE.sample, DISTR.sample, DISTR.n_sample * distr->dim * sizeof(double) );
  }
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
_unur_distr_cvemp_free( struct unur_distr *distr )
{
  if( distr == NULL ) 
    return;
  COOKIE_CHECK(distr,CK_DISTR_CVEMP,RETURN_VOID);
  if (DISTR.sample) free( DISTR.sample );
  if (distr->name_str) free(distr->name_str);
  COOKIE_CLEAR(distr);
  free( distr );
} 
int
unur_distr_cvemp_set_data( struct unur_distr *distr, const double *sample, int n_sample )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, sample, UNUR_ERR_NULL );
  if (n_sample <= 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"sample size");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.sample = _unur_xmalloc( n_sample * distr->dim * sizeof(double) );
  if (!DISTR.sample) return UNUR_ERR_MALLOC;
  memcpy( DISTR.sample, sample, n_sample * distr->dim * sizeof(double) );
  DISTR.n_sample = n_sample;
  return UNUR_SUCCESS;
} 
int
unur_distr_cvemp_read_data( struct unur_distr *distr, const char *filename )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEMP, UNUR_ERR_DISTR_INVALID );
  DISTR.n_sample = _unur_read_data( filename, distr->dim, &(DISTR.sample) );
  return (DISTR.n_sample > 0) ? UNUR_SUCCESS : UNUR_ERR_DISTR_DATA;
} 
int 
unur_distr_cvemp_get_data( const struct unur_distr *distr, const double **sample )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEMP, 0 );
  *sample = (DISTR.sample) ? DISTR.sample : NULL;
  return DISTR.n_sample;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_cvemp_debug( const struct unur_distr *distr, const char *genid, unsigned printvector )
{
#define idx(k,l)  (k * distr->dim + l)
  FILE *LOG;
  int i,j;
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CVEMP,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous multivariate distribution (ie. a sample)\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);
  fprintf(LOG,"%s:\tdimension = %d\n",genid,distr->dim);
  if (DISTR.n_sample>0) {
    fprintf(LOG,"%s:\tsample size = %d\n",genid,DISTR.n_sample);
    if (printvector) {
      for (i=0; i<DISTR.n_sample; i++) {
	fprintf(LOG,"%s:\t( %.5f",genid,DISTR.sample[idx(i,0)]);
	for (j=1; j<distr->dim; j++) 
	  fprintf(LOG,", %.5f",DISTR.sample[idx(i,j)]);
	fprintf(LOG,")\n");
      }
    }
  }
  fprintf(LOG,"%s:\n",genid);
#undef idx
} 
#endif    
