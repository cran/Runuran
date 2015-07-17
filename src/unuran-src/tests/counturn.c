/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"
static char test_name[] = "CountURN";
static long urng_counter = 0;                
#ifdef UNUR_URNG_UNURAN
static double (*urng_to_use)(void*);         
static double
_urng_with_counter(void *params)
{
  ++urng_counter;
  return urng_to_use(params);
} 
#endif  
int
unur_test_count_urn( struct unur_gen *gen, int samplesize, int verbosity, FILE *out )
{
  long j;
  UNUR_URNG *urng_aux;
  _unur_check_NULL(test_name,gen,-1);
  urng_counter = 0;
  urng_aux = gen->urng_aux;
#ifdef UNUR_URNG_UNURAN
  urng_to_use = gen->urng->sampleunif;
  gen->urng->sampleunif = _urng_with_counter;
  if (gen->urng_aux) gen->urng_aux = gen->urng;
#else
  if (verbosity)
    fprintf(out,"\nCOUNT: ---  (cannot count URNs)\n");
  return -1;
#endif
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    for( j=0; j<samplesize; j++ )
      _unur_sample_discr(gen);
    break;
  case UNUR_METH_CONT:
  case UNUR_METH_CEMP:
    for( j=0; j<samplesize; j++ )
      _unur_sample_cont(gen);
    break;
  case UNUR_METH_VEC: 
    { 
      double *vec;
      int dim;
      dim = unur_get_dimension(gen);
      vec = _unur_xmalloc( dim * sizeof(double) );
      for( j=0; j<samplesize; j++ )
	_unur_sample_vec(gen,vec);
      free(vec);
    }
    break;
  default: 
    _unur_error(test_name,UNUR_ERR_GENERIC,"method unknown!");
    return -1;
  }
#ifdef UNUR_URNG_UNURAN
  gen->urng->sampleunif = urng_to_use;
#endif
  gen->urng_aux = urng_aux;
  if (verbosity) {
    fprintf(out,"\nCOUNT: %g urng per generated number (total = %ld)\n",
	   ((double)urng_counter)/((double) samplesize),urng_counter);
  }
  return urng_counter;
} 
