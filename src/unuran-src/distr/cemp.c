/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "cemp.h"
#define DISTR distr->data.cemp
static void _unur_distr_cemp_free( struct unur_distr *distr );
struct unur_distr *
unur_distr_cemp_new( void )
{
  register struct unur_distr *distr;
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;
  COOKIE_SET(distr,CK_DISTR_CEMP);
  distr->type = UNUR_DISTR_CEMP;
  distr->id = UNUR_DISTR_GENERIC;
  distr->dim = 1;   
  distr->name = "(empirical)";
  distr->name_str = NULL;
  distr->destroy = _unur_distr_cemp_free;
  distr->clone = _unur_distr_cemp_clone;
  DISTR.sample    = NULL;   
  DISTR.n_sample  = 0;      
  DISTR.n_hist    = 0;          
  DISTR.hist_prob = NULL;       
  DISTR.hist_bins = NULL;       
  DISTR.hmin      = -INFINITY;  
  DISTR.hmax      = INFINITY;          
  DISTR.hist_bins = NULL;       
  return distr;
} 
struct unur_distr *
_unur_distr_cemp_clone( const struct unur_distr *distr )
{
#define CLONE clone->data.cemp
  struct unur_distr *clone;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CEMP, NULL );
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  memcpy( clone, distr, sizeof( struct unur_distr ) );
  if (DISTR.sample) {
    CLONE.sample = _unur_xmalloc( DISTR.n_sample * sizeof(double) );
    memcpy( CLONE.sample, DISTR.sample, DISTR.n_sample * sizeof(double) );
  }
  if (DISTR.hist_prob) {
    CLONE.hist_prob = _unur_xmalloc( DISTR.n_hist * sizeof(double) );
    memcpy( CLONE.hist_prob, DISTR.hist_prob, DISTR.n_hist * sizeof(double) );
  }
  if (DISTR.hist_bins) {
    CLONE.hist_bins = _unur_xmalloc( (DISTR.n_hist+1) * sizeof(double) );
    memcpy( CLONE.hist_bins, DISTR.hist_bins, (DISTR.n_hist+1) * sizeof(double) );
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
_unur_distr_cemp_free( struct unur_distr *distr )
{
  if( distr == NULL ) 
    return;
  COOKIE_CHECK(distr,CK_DISTR_CEMP,RETURN_VOID);
  if (DISTR.sample)    free( DISTR.sample );
  if (DISTR.hist_prob) free( DISTR.hist_prob );
  if (DISTR.hist_bins) free( DISTR.hist_bins );
  if (distr->name_str) free(distr->name_str);
  COOKIE_CLEAR(distr);
  free( distr );
} 
int
unur_distr_cemp_set_data( struct unur_distr *distr, const double *sample, int n_sample )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, sample, UNUR_ERR_NULL );
  if (n_sample <= 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"sample size");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.sample = _unur_xmalloc( n_sample * sizeof(double) );
  if (!DISTR.sample) return UNUR_ERR_MALLOC;
  memcpy( DISTR.sample, sample, n_sample * sizeof(double) );
  DISTR.n_sample = n_sample;
  return UNUR_SUCCESS;
} 
int
unur_distr_cemp_read_data( struct unur_distr *distr, const char *filename )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  DISTR.n_sample = _unur_read_data( filename, 1, &(DISTR.sample) );
  return (DISTR.n_sample > 0) ? UNUR_SUCCESS : UNUR_ERR_DISTR_DATA;
} 
int 
unur_distr_cemp_get_data( const struct unur_distr *distr, const double **sample )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CEMP, 0 );
  *sample = (DISTR.sample) ? DISTR.sample : NULL;
  return DISTR.n_sample;
} 
int
unur_distr_cemp_set_hist( struct unur_distr *distr, const double *prob, 
			  int n_prob, double xmin, double xmax )
{
  int rcode;
  if ((rcode = unur_distr_cemp_set_hist_domain(distr, xmin, xmax)) != UNUR_SUCCESS)
    return rcode;
  if ((rcode = unur_distr_cemp_set_hist_prob(distr, prob, n_prob)) != UNUR_SUCCESS) {
    distr->set &= ~UNUR_DISTR_SET_DOMAIN;
    return rcode;
  }
  return UNUR_SUCCESS;
} 
int
unur_distr_cemp_set_hist_prob( struct unur_distr *distr, const double *prob, int n_prob )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, prob, UNUR_ERR_NULL );
  if (n_prob <= 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"histogram size");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.hist_prob = _unur_xmalloc( n_prob * sizeof(double) );
  if (!DISTR.hist_prob) return UNUR_ERR_MALLOC;
  memcpy( DISTR.hist_prob, prob, n_prob * sizeof(double) );
  DISTR.n_hist = n_prob;
  return UNUR_SUCCESS;
} 
int
unur_distr_cemp_set_hist_domain( struct unur_distr *distr, double xmin, double xmax )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  if (xmin >= xmax) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"histogram, min >= max");
    return UNUR_ERR_DISTR_SET;
  }
  if (!_unur_isfinite(xmin) || !_unur_isfinite(xmax)) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"histogram, unbounded domain");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.hmin = xmin;
  DISTR.hmax = xmax;
  distr->set |= UNUR_DISTR_SET_DOMAIN;
  return UNUR_SUCCESS;
} 
int
unur_distr_cemp_set_hist_bins( struct unur_distr *distr, const double *bins, int n_bins )
{
  int i;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, bins, UNUR_ERR_NULL );
  if (DISTR.hist_prob == NULL) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"probabilities of histogram not set");
    return UNUR_ERR_DISTR_SET;
  }
  if (n_bins != (DISTR.n_hist+1)) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"histogram size");
    return UNUR_ERR_DISTR_SET;
  }
  for( i=1; i<n_bins; i++ )
    if (bins[i] <= bins[i-1]) {
      _unur_error(distr->name,UNUR_ERR_DISTR_SET,"bins not strictly increasing");
      return UNUR_ERR_DISTR_SET;
    }
  if (unur_distr_cemp_set_hist_domain(distr, bins[0], bins[n_bins-1]) != UNUR_SUCCESS)
    return UNUR_ERR_DISTR_SET;
  DISTR.hist_bins = _unur_xmalloc( n_bins * sizeof(double) );
  if (!DISTR.hist_bins) return UNUR_ERR_MALLOC;
  memcpy( DISTR.hist_bins, bins, n_bins * sizeof(double) );
  distr->set |= UNUR_DISTR_SET_DOMAIN;
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_cemp_debug( const struct unur_distr *distr, const char *genid, unsigned printvector )
{
  FILE *LOG;
  int i;
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CEMP,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous univariate distribution (ie. a sample)\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);
  if (DISTR.n_sample>0) {
    fprintf(LOG,"%s:\tsample size = %d",genid,DISTR.n_sample);
    if (printvector) {
      for (i=0; i<DISTR.n_sample; i++) {
	if (i%10 == 0)
	  fprintf(LOG,"\n%s:\t",genid);
	fprintf(LOG,"  %.5f",DISTR.sample[i]);
      }
    }
    fprintf(LOG,"\n%s:\n",genid);
  }
  if (DISTR.n_hist>0) {
    fprintf(LOG,"%s:\thistogram: #bins = %d, ",genid,DISTR.n_hist);
    fprintf(LOG,"min = %g, max = %g", DISTR.hmin, DISTR.hmax);
    if (DISTR.hist_bins) {
      fprintf(LOG," (bins with different width)\n");
      if (printvector) {
	fprintf(LOG,"%s:\t> bins (breaks) = ",genid);
	for (i=0; i<=DISTR.n_hist; i++) {
	  if (i%10 == 0)
	    fprintf(LOG,"\n%s:\t",genid);
	  fprintf(LOG,"  %.5f",DISTR.hist_bins[i]);
	}
	fprintf(LOG,"\n%s:\n",genid);
      }
    }
    else {
      fprintf(LOG,", width = %g (equally spaced)\n", (DISTR.hmax-DISTR.hmin)/DISTR.n_hist);
    }
    if (printvector) {
      fprintf(LOG,"%s:\t> bin probabilities = ",genid);
      for (i=0; i<DISTR.n_hist; i++) {
	if (i%10 == 0)
	  fprintf(LOG,"\n%s:\t",genid);
	fprintf(LOG,"  %.5f",DISTR.hist_prob[i]);
      }
      fprintf(LOG,"\n%s:\n",genid);
    }
  }
} 
#endif    
