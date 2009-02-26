/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"
static char test_name[] = "Moments";
int
unur_test_moments( UNUR_GEN *gen, double *moments, int n_moments, int samplesize,
		   int verbosity, FILE *out )
{
#define idx(d,n) ((d)*(n_moments+1)+(n))
  double an, an1, dx, dx2;
  int n, mom, d, dim;
  double *x;
  _unur_check_NULL(test_name, gen, UNUR_ERR_NULL);
  if (! ( ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_CONT)  ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_VEC) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute moments for distribution");
    return UNUR_ERR_GENERIC;
  }
  CHECK_NULL(moments,0);
  if (n_moments <= 0 || n_moments > 4) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"number of moments < 1 or > 4");
    return UNUR_ERR_GENERIC;
  }
  if (samplesize < 10) 
    samplesize = 10;
  dim = 1;
  if ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_VEC) dim = gen->distr->dim;
  x = _unur_xmalloc(dim * sizeof(double));
  for (d=0; d<dim; d++) {
    moments[idx(d,0)] = 1.; 
    for (mom = 1; mom <= n_moments; mom++ )
      moments[idx(d,mom)] = 0.;
  }
  for (n=1; n<=samplesize; n++) {
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x[0] = _unur_sample_discr(gen); break;
    case UNUR_METH_CONT:
      x[0] = _unur_sample_cont(gen); break;
    case UNUR_METH_VEC:
      _unur_sample_vec(gen, x); break;
    }
    for (d=0; d<dim; d++) {
      an = (double)n;
      an1 = an-1.;
      dx = (x[d] - moments[idx(d,1)]) / an;
      dx2 = dx * dx;
      switch (n_moments) {
      case 4:
        moments[idx(d,4)] -= dx * (4.*moments[idx(d,3)] - dx * (6.*moments[idx(d,2)] + an1*(1. + an1*an1*an1)*dx2));
      case 3:
        moments[idx(d,3)] -= dx * (3.*moments[idx(d,2)] - an*an1*(an-2.)*dx2);
      case 2:
        moments[idx(d,2)] += an * an1 * dx2;
      case 1:
        moments[idx(d,1)] += dx;
      }
    }
  }
  for (d=0; d<dim; d++) {
    for (mom = 2; mom <= n_moments; mom++ )
      moments[idx(d,mom)] /= samplesize;
    if (verbosity) {
      if (dim==1) fprintf(out,"\nCentral MOMENTS:\n");
      else        fprintf(out,"\nCentral MOMENTS for dimension #%d:\n", d);
      for (mom = 1; mom <= n_moments; mom++ )
        fprintf(out,"\t[%d] =\t%g\n",mom,moments[idx(d,mom)]);
      fprintf(out,"\n");
    }
  }
  free(x);
  return UNUR_SUCCESS;
#undef idx  
} 
