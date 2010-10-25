/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include <methods/cstd.h>
#include <methods/unif.h>
#include <distr/distr.h>
#include <distributions/unur_distributions.h>
#include <parser/parser.h>
#include <urng/urng.h>
#include "unuran_tests.h"
#if defined(HAVE_GETTIMEOFDAY) && defined(HAVE_SYS_TIME_H)
#include <sys/time.h>
static struct timeval tv;
#define _unur_get_time() ( gettimeofday(&tv, NULL), ((tv).tv_sec * 1.e6 + (tv).tv_usec) )
#else
#include <time.h>
#define _unur_get_time() ( (1.e6 * clock()) / CLOCKS_PER_SEC )
#endif
static char test_name[] = "Timing";
static double unur_test_timing_total_run( const struct unur_par *par, int samplesize, int repeat );
inline static int
compare_doubles (const void *a, const void *b)
{ 
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}
struct unur_gen*
unur_test_timing( struct unur_par *par, 
		  int log10_samplesize, 
		  double *time_setup,
		  double *time_sample,
		  int verbosity,
		  FILE *out )
{
  struct unur_gen *gen;
  int k;
  double x;
  double *vec = NULL;
  double time_uniform, time_exponential;
  double time_start, *time_gen;
  long samples, samplesize, log10_samples;
  _unur_check_NULL(test_name,par,NULL);
  if (log10_samplesize < 2) log10_samplesize = 2;
  time_gen = _unur_xmalloc((log10_samplesize+1) * sizeof(double));
  time_uniform = unur_test_timing_uniform( par,log10_samplesize );
  time_exponential = unur_test_timing_exponential( par,log10_samplesize );
  if (par->distr && _unur_gen_is_vec(par))
    vec = _unur_xmalloc( par->distr->dim * sizeof(double) );
  time_start = _unur_get_time();
  gen = _unur_init(par);
  *time_setup = _unur_get_time();
  if (!gen) {
    free (time_gen);
    if (vec) free(vec);
    return NULL;
  }
  samplesize = 10;
  samples = 0;
  for( log10_samples=1; log10_samples<=log10_samplesize; log10_samples++ ) {
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( ; samples < samplesize; samples++ )
	k = unur_sample_discr(gen);
      break;
    case UNUR_METH_CONT:
    case UNUR_METH_CEMP:
      for( ; samples < samplesize; samples++ )
	x = unur_sample_cont(gen);
      break;
    case UNUR_METH_VEC:
      for( ; samples < samplesize; samples++ )
	unur_sample_vec(gen,vec);
      break;
    default: 
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      free(time_gen);
      if (vec) free(vec);
      return NULL;
    }
    time_gen[log10_samples] = _unur_get_time();
    samplesize *= 10;
  }
  *time_sample = (time_gen[log10_samplesize] - time_gen[log10_samplesize-1]) / (0.09 * samplesize);
  samplesize = 1;
  for( log10_samples=1; log10_samples<=log10_samplesize; log10_samples++ ) {
    samplesize *= 10;
    time_gen[log10_samples] = (time_gen[log10_samples] - time_start) / samplesize;
  }
  *time_setup -= time_start;
  if (verbosity) {
    fprintf(out,"\nTIMING:\t\t    usec \t relative to \t relative to\n");
    fprintf(out,"\t\t\t\t uniform\t exponential\n\n");
    fprintf(out,"   setup time:\t    %#g \t %#g \t %#g\n",
	    (*time_setup),
	    (*time_setup)/time_uniform,
	    (*time_setup)/time_exponential);
    fprintf(out,"   generation time: %#g \t %#g \t %#g\n",
	    (*time_sample),
	    (*time_sample)/time_uniform,
	    (*time_sample)/time_exponential);
    fprintf(out,"\n   average generation time for samplesize:\n");
    for( log10_samples=1; log10_samples<=log10_samplesize; log10_samples++ )
      fprintf(out,"\t10^%ld:\t    %#g \t %#g \t %#g\n",log10_samples,
	      time_gen[log10_samples],
	      time_gen[log10_samples]/time_uniform,
	      time_gen[log10_samples]/time_exponential);
  }
  free(time_gen);
  if (vec) free(vec);
  return gen;
} 
double 
unur_test_timing_R( struct unur_par *par, const char *distrstr, const char *methodstr,
		    double log10_samplesize, double *time_setup, double *time_marginal )
{
  const int n_steps = 2;    
  const int n_reps = 10;    
  struct unur_distr *distr_tmp = NULL; 
  struct unur_par   *par_tmp   = NULL; 
  struct unur_gen   *gen_tmp   = NULL; 
  struct unur_slist *mlist     = NULL; 
  int k;
  double x;
  double *vec = NULL;
  double time_start, *time_gen = NULL;
  int sample, rep;
  long samplesize, n;
  double sx = 0.;         
  double sy = 0.;
  double sxx = 0.;
  double syy = 0.;
  double sxy = 0.;
  double Rsq = -100.;     
  *time_setup = -100.;
  *time_marginal = -100.;
  Rsq = -100.;
  if (log10_samplesize < 2) log10_samplesize = 2;
  if (par == NULL) {
    distr_tmp = unur_str2distr(distrstr);
    if (distr_tmp == NULL) goto error;
    par = _unur_str2par( distr_tmp, methodstr, &mlist );
    if (par == NULL) goto error;
  }
  time_gen = _unur_xmalloc( n_reps * sizeof(double) );
  if (par->distr && _unur_gen_is_vec(par))
    vec = _unur_xmalloc( par->distr->dim * sizeof(double) );
  for (sample=0; sample<n_steps; sample++) {
    samplesize = (long) exp(M_LN10 * (1. + sample * (log10_samplesize - 1.) / (n_steps - 1.)));
    for (rep=0; rep<n_reps; rep++) {
      par_tmp = _unur_par_clone(par);
      time_start = _unur_get_time();
      gen_tmp = _unur_init(par_tmp);
      if (!gen_tmp) goto error;  
      switch (gen_tmp->method & UNUR_MASK_TYPE) {
      case UNUR_METH_DISCR:
	for( n=0; n<samplesize; n++ )
	  k = unur_sample_discr(gen_tmp);
	break;
      case UNUR_METH_CONT:
	for( n=0; n<samplesize; n++ )
	  x = unur_sample_cont(gen_tmp);
	break;
      case UNUR_METH_VEC:
	for( n=0; n<samplesize; n++ )
	  unur_sample_vec(gen_tmp,vec);
	break;
      default: 
	_unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      }
      time_gen[rep]= _unur_get_time() - time_start;
      unur_free(gen_tmp);
    }
    qsort( time_gen, (size_t)n_reps, sizeof(double), compare_doubles);
    for (rep=2; rep<n_reps-3; rep++) {
      sx += samplesize;
      sy += time_gen[rep];
      sxx += ((double)samplesize) * ((double)samplesize);
      syy += time_gen[rep] * time_gen[rep];
      sxy += ((double)samplesize) * time_gen[rep];
    }
  }
  n = n_steps * (n_reps-5);
  *time_marginal = ( n*sxy - sx*sy ) / ( n*sxx - sx*sx );
  *time_setup = sy/n - *time_marginal * sx/n;
  Rsq = ( n*sxy - sx*sy ) / sqrt( (n*sxx - sx*sx) * (n*syy - sy*sy) );
 error:
  if (distr_tmp) unur_distr_free(distr_tmp);
  if (par)       _unur_par_free(par);
  if (mlist)     _unur_slist_free(mlist);
  if (time_gen)  free(time_gen);
  if (vec)       free(vec);
  return Rsq;
} 
double 
unur_test_timing_total( const UNUR_PAR *par, int samplesize, double avg_duration )
{
  double time_pilot, time_result;
  int size_pilot, size_result;
  int repeat_pilot, repeat_result;
  double time_2nd, d, k;
  _unur_check_NULL(test_name,par,-1.);
  if (samplesize < 0) return -1.;
  avg_duration = (avg_duration < 1.e-3) ? 1000. : 1.e6 * avg_duration;
  repeat_pilot = 11 - (int)(log((double)samplesize)/M_LN2);
  if (repeat_pilot<1) repeat_pilot = 1;
  size_pilot = _unur_min(samplesize,1000);
  time_pilot = unur_test_timing_total_run(par, size_pilot, repeat_pilot);
  if (time_pilot < 0) return -1.;   
  if (samplesize > 1000) {
    time_2nd = unur_test_timing_total_run(par, 2*size_pilot, repeat_pilot);
    if (time_2nd < 0) return -1.;
    d = 2*time_pilot - time_2nd;
    if (d<0.) d=0.;
    k = (time_2nd - time_pilot)/size_pilot;
    if (k<=0.) k = time_pilot/size_pilot;
    time_pilot = d + samplesize * k;
  }
  else {
    d = 0;
    k = time_pilot / size_pilot;
  }
  repeat_result = (int) (avg_duration / time_pilot);
  if (repeat_result > 1000) repeat_result = 1000;
  size_result = samplesize;
  if (repeat_result >= 1) {
    repeat_result = _unur_max(4,repeat_result);
    if (repeat_result <= repeat_pilot && size_result == size_pilot) {
      time_result = time_pilot;
    }
    else {
      time_result =  unur_test_timing_total_run(par,size_result,repeat_result);
    }
  }
  else {
    repeat_result = 4;
    size_result = (int) ((avg_duration - d)/k);
    size_result /= 2;
    time_result =  unur_test_timing_total_run(par,size_result,repeat_result);
    time_2nd =  unur_test_timing_total_run(par,2*size_result,repeat_result);
    d = 2*time_result - time_2nd;
    if (d<0.) d=0.;
    k = (time_2nd - time_result)/size_result;
    if (k<=0.) k = time_result/size_result;
    time_result = d + samplesize * k;
  }      
  return time_result;
} 
double unur_test_timing_total_run( const struct unur_par *par, int samplesize, int n_repeat )
{
  struct unur_par *par_tmp;    
  struct unur_gen *gen_tmp;    
  double *time;
  int n, rep;
  double time_total;           
  double time_start;
  int i,k;
  double x;
  double *vec = NULL;
  _unur_check_NULL(test_name,par,-1.);
  if (samplesize < 0 || n_repeat < 1) 
    return -1.;
  time = _unur_xmalloc( n_repeat * sizeof(double) );
  if (par->distr && _unur_gen_is_vec(par))
    vec = _unur_xmalloc( par->distr->dim * sizeof(double) );
  for (rep = 0; rep < n_repeat; rep++) {
    par_tmp = _unur_par_clone(par);
    time_start = _unur_get_time();
    gen_tmp = _unur_init(par_tmp);
    if (!gen_tmp) {  
      if (vec) free(vec);
      free(time);
      for (x=0,i=0; i<100000; i++) x+=i;  
      return -1.;
    }
    switch (gen_tmp->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( n=0; n<samplesize; n++ )
	k = unur_sample_discr(gen_tmp);
      break;
    case UNUR_METH_CONT:
      for( n=0; n<samplesize; n++ )
	x = unur_sample_cont(gen_tmp);
      break;
    case UNUR_METH_VEC:
      for( n=0; n<samplesize; n++ )
	unur_sample_vec(gen_tmp,vec);
      break;
    default: 
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }
    time[rep]= _unur_get_time() - time_start;
    unur_free(gen_tmp);
  }
  qsort( time, (size_t)n_repeat, sizeof(double), compare_doubles);
  time_total = time[n_repeat/2];
  if (vec) free(vec);
  free(time);
  return time_total;
} 
double
unur_test_timing_uniform( const struct unur_par *par, int log10_samplesize )
{
#define TIMING_REPETITIONS (21)
  struct unur_gen *gen_urng;
  static double uniform_time = -1.;
  double time[TIMING_REPETITIONS];
  double x;
  int j,n;
  if (uniform_time <= 0.) {  
    int samplesize = 1;
    for( j=0; j<log10_samplesize; j++ )
      samplesize *= 10;
    gen_urng = unur_init( unur_unif_new(NULL) );
    _unur_check_NULL( test_name,gen_urng,-1. );
    unur_chg_urng(gen_urng,par->urng);
    for( n=0; n<TIMING_REPETITIONS; n++ ) {
      time[n] = _unur_get_time();
      for( j=0; j<samplesize; j++ )
	x = unur_sample_cont(gen_urng);
      time[n] = (_unur_get_time() - time[n])/samplesize;
    }
    qsort( time, (size_t)TIMING_REPETITIONS, sizeof(double), compare_doubles);
    uniform_time = time[TIMING_REPETITIONS/2];
    _unur_free(gen_urng);
  }
  return uniform_time;
#undef TIMING_REPETITIONS
} 
double
unur_test_timing_exponential( const struct unur_par *par, int log10_samplesize )
{
#define TIMING_REPETITIONS (21)
  struct unur_distr *unit_distr;
  struct unur_par   *unit_par;
  struct unur_gen   *unit_gen;
  static double exponential_time = -1.;
  double time[TIMING_REPETITIONS];
  double x;
  int j,n;
  if (exponential_time <= 0.) {  
    int samplesize = 1;
    for( j=0; j<log10_samplesize; j++ )
      samplesize *= 10;
    unit_distr = unur_distr_exponential(NULL,0);
    unit_par = unur_cstd_new(unit_distr);
    unur_cstd_set_variant(unit_par,UNUR_STDGEN_INVERSION);
    unit_gen = unur_init(unit_par); 
    _unur_check_NULL( test_name,unit_gen,-1. );
    unur_chg_urng(unit_gen,par->urng);
    for( n=0; n<TIMING_REPETITIONS; n++ ) {
      time[n] = _unur_get_time();
      for( j=0; j<samplesize; j++ )
	x = unur_sample_cont(unit_gen);
      time[n] = (_unur_get_time() - time[n])/samplesize;
    }
    qsort( time, (size_t)TIMING_REPETITIONS, sizeof(double), compare_doubles);
    exponential_time = time[TIMING_REPETITIONS/2];
    unur_distr_free(unit_distr);
    unur_free(unit_gen);
  }
  return exponential_time;
#undef TIMING_REPETITIONS
} 
