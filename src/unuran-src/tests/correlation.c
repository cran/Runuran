/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/cont.h>
#include <utils/matrix_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"
static char test_name[] = "Correlation";
#define CORR_DEFAULT_SAMPLESIZE 10000
#define CORR_MAX_SAMPLESIZE 10000000
double
unur_test_correlation( UNUR_GEN *genx, UNUR_GEN *geny, int samplesize, int verbosity, FILE *out )
{
  double x  =0., y =0.;   
  double mx =0., my=0.;  
  double dx =0., dy=0.;  
  double sx =0., sy=0.;  
  double sxy=0.;         
  double factor;
  int n;
  _unur_check_NULL(test_name,genx,-3.);
  _unur_check_NULL(test_name,geny,-3.);
  if (! ( ((genx->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((genx->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
         "dont know how to compute correlation coefficient for distribution");
    return -2.;
  }
  if (! ( ((geny->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((geny->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
         "dont know how to compute correlation coefficient for distribution");
    return -2.;
  }
  if( samplesize <= 0 ) samplesize = CORR_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CORR_MAX_SAMPLESIZE );
  for (n=1; n<=samplesize; n++) {
    switch (genx->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = _unur_sample_discr(genx); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(genx); break;
    }
    switch (geny->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      y = _unur_sample_discr(geny); break;
    case UNUR_METH_CONT:
      y = _unur_sample_cont(geny); break;
    }
    factor = (double) ( n*(n-1) );
    dx = (x - mx) / n;
    dy = (y - my) / n;
    mx += dx;
    my += dy;
    sx  += factor * dx*dx;
    sy  += factor * dy*dy;
    sxy += factor * dx*dy;
  }
  if (verbosity) {
    fprintf(out,"\nCorrelation coefficient: %g\n\n", sxy/(sqrt(sx*sy)) );
  }
  return ( sxy/(sqrt(sx*sy)) );
} 
int
unur_test_cvec_rankcorr( double *rc, struct unur_gen *gen, int samplesize, int verbose, FILE *out )
{
#define DISTR   gen->distr->data.cvec
#define idx(a,b) ((a)*dim+(b))
  double *mx;   
  double *dx;   
  double factor;
  int dim;         
  UNUR_DISTR **marginals;  
  UNUR_FUNCT_CONT **marginal_cdf;  
  double *X;             
  double *U;             
  int i, j, n;     
  CHECK_NULL(gen,UNUR_ERR_NULL);
  if (verbose >= 1)
    fprintf(out,"\nRank correlations of random vector:\n");
  if( samplesize <= 0 ) samplesize = CORR_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CORR_MAX_SAMPLESIZE );
  dim = gen->distr->dim;
  if (dim < 1) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"distribution dimension < 1 ?");
    return UNUR_ERR_GENERIC;
  }
  if ( (gen->method & UNUR_MASK_TYPE) != UNUR_METH_VEC ) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
		"rank correlation coefficients cannot be computed");
    return UNUR_ERR_GENERIC;
  }
  if (DISTR.marginals==NULL) {
    _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"marginal distributions");
    return UNUR_ERR_DISTR_REQUIRED; }
  marginals = _unur_xmalloc(dim * sizeof(UNUR_DISTR *));
  marginal_cdf = _unur_xmalloc(dim * sizeof(UNUR_FUNCT_CONT *));
  for (i=0; i<dim; i++) {
    marginals[i] = DISTR.marginals[i];
    marginal_cdf[i] = unur_distr_cont_get_cdf(DISTR.marginals[i]);
    if (marginals[i]==NULL || marginal_cdf[i]==NULL) {
      _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"CDF of continuous marginal");
      free (marginals);  free (marginal_cdf);
      return UNUR_ERR_DISTR_REQUIRED; }
  }
  X   = _unur_xmalloc( dim * sizeof(double));
  U   = _unur_xmalloc( dim * sizeof(double));
  mx  = _unur_xmalloc( dim * sizeof(double));
  dx  = _unur_xmalloc( dim * sizeof(double));
  for (i=0; i<dim; i++)
    mx[i] = dx[i] = 0.;
  for (i=0; i<dim*dim; i++) 
    rc[i] = 0.;
  for (n=1; n<=samplesize; n++) {
    factor = ((double)n)*(n-1.);
    _unur_sample_vec(gen, X);
    for (i=0; i<dim; i++) {
      U[i] = marginal_cdf[i](X[i],marginals[i]);
      dx[i] = (U[i] - mx[i]) / n;
      mx[i] += dx[i];
    }
    for (i=0; i<dim; i++)
      for (j=i; j<dim; j++)
	rc[idx(i,j)] += factor * dx[i] * dx[j];
  }
  for (i=0; i<dim; i++) {
    for (j=0; j<i; j++)
      rc[idx(i,j)] = rc[idx(j,i)];
    for (j=i+1; j<dim; j++)
      rc[idx(i,j)] /= sqrt(rc[idx(i,i)] * rc[idx(j,j)]);
    rc[idx(i,i)] = 1.;
  }
  if (verbose >= 1)
    _unur_matrix_print_matrix( dim, rc, "rank correlation =", 
			       out, "", "\t");
  if (X)    free(X);
  if (U)    free(U);
  if (mx)   free(mx);
  if (dx)   free(dx);
  if (marginals)  free (marginals);
  if (marginal_cdf)  free (marginal_cdf);
  return UNUR_SUCCESS;
#undef DISTR
#undef idx
} 
