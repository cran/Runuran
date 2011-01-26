/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <limits.h>
#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include <distr/cont.h>
#include <distr/cvec.h>
#include <distr/discr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <specfunct/unur_specfunct_source.h>
#include <utils/matrix_source.h>
#include "unuran_tests.h"
#define CHI2_SAMPLEFAC  40
#define CHI2_CLASSMIN_DEFAULT  20
#define CHI2_INTERVALS_DEFAULT 50
#define CHI2_DEFAULT_SAMPLESIZE 10000
#define CHI2_MAX_SAMPLESIZE 1000000
#define CHI2_MAX_TOTALINTERVALS 1000000
static char test_name[] = "Chi^2-Test";
static double _unur_test_chi2_discr( struct unur_gen *gen, int samplesize, int classmin,
				     int verbose, FILE *out );
static double _unur_test_chi2_cont( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				    int verbose, FILE *out );
static double _unur_test_chi2_cemp( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				    int verbose, FILE *out );
static double _unur_test_chi2_vec( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				   int verbose, FILE *out );
static double _unur_test_chi2_cvemp( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				     int verbose, FILE *out );
static double _unur_test_chi2test( double *prob, int *observed, int len, int classmin,
				   int verbose, FILE *out );
double
unur_test_chi2( struct unur_gen *gen,
    int n_intervals,
    int samplesize,
    int classmin,
    int verbose,
    FILE *out )
{
  _unur_check_NULL(test_name,gen,-1.);
  if (verbose >= 1)
    fprintf(out,"\nGOODNESS-OF-FIT TESTS:\n");
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    return _unur_test_chi2_discr(gen, samplesize, classmin, verbose, out);
  case UNUR_METH_CONT:
    return _unur_test_chi2_cont(gen, n_intervals, samplesize, classmin, verbose, out);
  case UNUR_METH_CEMP:
    return _unur_test_chi2_cemp(gen, n_intervals, samplesize, classmin, verbose, out);
  case UNUR_METH_VEC:
    return _unur_test_chi2_vec(gen, n_intervals, samplesize, classmin, verbose, out);
  case UNUR_METH_CVEMP:
    return _unur_test_chi2_cvemp(gen, n_intervals, samplesize, classmin, verbose, out);
  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"Not implemented for such distributions!");
    return -1.;
  }
} 
double
_unur_test_chi2_discr( struct unur_gen *gen,
           int samplesize,
           int classmin,
           int verbose,
           FILE *out)
{
#define DISTR   gen->distr->data.discr
  double *pv;           
  int n_pv;             
  int *observed;        
  double pval;          
  int had_PV;           
  int i,j;
  CHECK_NULL(gen,-1.);
  if (DISTR.pv == NULL) {
    had_PV = FALSE;
    if (!unur_distr_discr_make_pv( gen->distr )) {
      return -2.;
    }
  }
  else
    had_PV = TRUE;
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;
  observed = _unur_xmalloc( n_pv * sizeof(int));
  for( i=0; i<n_pv; i++ )
    observed[i] = 0;
  if( samplesize <= 0 ) {
    samplesize = (INT_MAX/n_pv > n_pv) ? n_pv * n_pv : 1000000;
    samplesize = _unur_max(samplesize,1000000);
  }
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );
  for( i=0; i<samplesize; i++ ) {
    j = _unur_sample_discr(gen);
    if (verbose >= 3) fprintf(out,"i = %d\n",j);
    j -= DISTR.domain[0];
    if (j >= 0 && j < n_pv) 
      ++observed[j];
  }
  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for discrete distribution with given probability vector:");
    fprintf(out,"\n  length     = %d\n",n_pv);
  }
  pval = _unur_test_chi2test(pv, observed, n_pv, classmin, verbose, out);
  free(observed);
  if (!had_PV) {
    free (DISTR.pv);
    DISTR.pv = NULL;
    DISTR.n_pv = 0;
  }
  return pval;
#undef DISTR
} 
double
_unur_test_chi2_cont( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
{
#define DISTR   gen->distr->data.cont
  double F, Fl, Fr, Fdelta;  
  UNUR_FUNCT_CONT *cdf;      
  int *observed;             
  double pval;               
  int i,j;
  CHECK_NULL(gen,-1.);
  cdf = DISTR.cdf;
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required for continuous random variates!");
    return -2.;
  }
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;
  observed = _unur_xmalloc( n_intervals * sizeof(int));
  for( i=0; i<n_intervals; i++ )
    observed[i] = 0;
  if( samplesize <= 0 )
    samplesize = (INT_MAX/n_intervals > n_intervals) ? n_intervals*n_intervals : INT_MAX;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    Fl = (DISTR.trunc[0] <= -INFINITY) ? 0. : cdf(DISTR.trunc[0], gen->distr);
    Fr = (DISTR.trunc[1] >=  INFINITY) ? 1. : cdf(DISTR.trunc[1], gen->distr);
  }
  else {
    Fl = (DISTR.domain[0] <= -INFINITY) ? 0. : cdf(DISTR.domain[0], gen->distr);
    Fr = (DISTR.domain[1] >=  INFINITY) ? 1. : cdf(DISTR.domain[1], gen->distr);
  }
  Fdelta = Fr - Fl;
  if (Fdelta <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
    free (observed);
    return -1.;
  }
  for( i=0; i<samplesize; i++ ) {
    if (verbose >= 3) {
      double x = _unur_sample_cont(gen);
      F = cdf( x, gen->distr );
      fprintf(out,"x = %g\n",x);
    }
    else {
      F = cdf( _unur_sample_cont(gen), gen->distr );
    }
    F = (F-Fl)/Fdelta;
    j = (int)(n_intervals * F);
    if (j > n_intervals) {
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) > Fmax (out of domain).");
      j = n_intervals-1;
    }
    if (j >= n_intervals)    
      j = n_intervals-1;
    if (j < 0 ) {           
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) < 0 (out of domain).");
      j = 0;
    }
    ++observed[j];
  }
  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for continuous distribution:");
    fprintf(out,"\n  intervals  = %d\n",n_intervals);
  }
  pval = _unur_test_chi2test(NULL, observed, n_intervals, classmin, verbose, out );
  free(observed);
  return pval;
#undef DISTR
} 
double
_unur_test_chi2_cemp( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
{
  UNUR_DISTR *distr_normal;  
  UNUR_FUNCT_CONT *cdf;      
  double F;                  
  int *observed;             
  double pval;               
  int i,j;
  CHECK_NULL(gen,-1.);
  distr_normal = unur_distr_normal( NULL, 0 );
  cdf = distr_normal->data.cont.cdf;
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;
  observed = _unur_xmalloc( n_intervals * sizeof(int));
  for( i=0; i<n_intervals; i++ )
    observed[i] = 0;
  if( samplesize <= 0 )
    samplesize = (INT_MAX/n_intervals > n_intervals) ? n_intervals*n_intervals : INT_MAX;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );
  for( i=0; i<samplesize; i++ ) {
    F = cdf( _unur_sample_cont(gen), distr_normal );
    j = (int)(n_intervals * F);
    if (j > n_intervals) {
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) > Fmax (out of domain).");
      j = n_intervals-1;
    }
    if (j >= n_intervals)    
      j = n_intervals-1;
    if (j < 0 ) {           
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) < 0 (out of domain).");
      j = 0;
    }
    ++observed[j];
  }
  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for continuous empirical distribution:");
    fprintf(out,"\n(Assumes standard normal distribution!)");
    fprintf(out,"\n  intervals  = %d\n",n_intervals);
  }
  pval = _unur_test_chi2test(NULL, observed, n_intervals, classmin, verbose, out );
  _unur_distr_free(distr_normal);
  free(observed);
  return pval;
} 
double
_unur_test_chi2_vec ( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
{
#define DISTR     gen->distr->data.cvec
#define idx(i,j)  ((i)*dim+(j))
  int dim;                       
  double *Fl = NULL;             
  double *Fr = NULL;
  double *Fdelta = NULL;
  UNUR_DISTR **marginals = NULL; 
  UNUR_FUNCT_CONT **marginal_cdf = NULL;  
  const double *mean = NULL;     
  const double *L = NULL;        
  double *Linv = NULL;           
  double Linv_det;               
  double *X = NULL;              
  double *U = NULL;              
  int *bm = NULL;                
  double pval, pval_min;         
  int i, j, k;                   
  CHECK_NULL(gen,-1.);
  if (DISTR.marginals==NULL) {
    _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"marginals");
    return -2.; 
  }
  if ((gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
      !(gen->distr->set & UNUR_DISTR_SET_COVAR_IDENT) )
    _unur_warning(test_name,UNUR_ERR_GENERIC,"correlated and domain truncated --> test might fail");
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;
  if( samplesize <= 0 ) samplesize = CHI2_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );
  dim = gen->distr->dim;  
  if (dim < 1) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"distribution dimension < 1 ?");
    return -1.;
  }
  mean = unur_distr_cvec_get_mean(gen->distr);
  L = unur_distr_cvec_get_cholesky(gen->distr);
  marginals = _unur_xmalloc(dim * sizeof(UNUR_DISTR *));
  marginal_cdf = _unur_xmalloc(dim * sizeof(UNUR_FUNCT_CONT *));
  for (i=0; i<dim; i++) {
    marginals[i] = DISTR.marginals[i];
    marginal_cdf[i] = unur_distr_cont_get_cdf(DISTR.marginals[i]);
    if (marginals[i]==NULL || marginal_cdf[i]==NULL) {
      _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"CDF of continuous standardized marginal");
      pval_min = -2.; goto free_memory;
    }
  }
  Fl  = _unur_xmalloc(dim * sizeof(double));
  Fr  = _unur_xmalloc(dim * sizeof(double));
  Fdelta = _unur_xmalloc(dim * sizeof(double));
  if (gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) {
    for (i=0; i<dim; i++) {
      Fl[i] = ( (!_unur_isfinite(DISTR.domainrect[2*i])) ? 0. 
		: marginal_cdf[i](DISTR.domainrect[2*i],marginals[i]) );
      Fr[i] = ( (!_unur_isfinite(DISTR.domainrect[2*i+1])) ? 1. 
		: marginal_cdf[i](DISTR.domainrect[2*i+1],marginals[i]) );
      Fdelta[i] = Fr[i] - Fl[i];
      if (Fdelta[i] <= 0.) {
	_unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
	pval_min = -1.; goto free_memory;
      }
    }
  }
  else {
    for (i=0; i<dim; i++) {
      Fl[i] = 0.;
      Fr[i] = 1.;
      Fdelta[i] = 1.;
    }
  }
  X   = _unur_xmalloc( dim * sizeof(double));
  U   = _unur_xmalloc( dim * sizeof(double));
  bm  = _unur_xmalloc( dim * n_intervals * sizeof(int)); 
  memset(bm , 0, dim * n_intervals  * sizeof(int));
  if (L != NULL) {
    Linv  = _unur_xmalloc( dim * dim * sizeof(double));
    if (_unur_matrix_invert_matrix (dim, L, Linv, &Linv_det) != UNUR_SUCCESS) {
      _unur_error(test_name,UNUR_ERR_DISTR_DATA,"cannot compute inverse of Cholesky factor");
      pval_min = -2.; goto free_memory;
    }
  }
  for( i=0; i<samplesize; i++ ) {
    _unur_sample_vec(gen, X);
    if (mean) { for (j=0; j<dim; j++)  X[j] -= mean[j]; }
    for (j=0; j<dim; j++) {
      double Z=0;
      if (Linv)
	for (k=0; k<=j; k++)  Z += Linv[idx(j,k)] * X[k]; 
      else
	Z = X[j];
      U[j] = (marginal_cdf[j](Z,marginals[j]) - Fl[j]) / Fdelta[j];
    }
    for (j=0; j<dim; j++) {
      int iv;
      iv = (int)( n_intervals * U[j] );
      if (iv>=n_intervals) iv = n_intervals-1;
      if (iv < 0) iv = 0;
      bm[j*n_intervals + iv] += 1;
    }
  }
  pval_min = 1.;
  for (j=0; j<dim; j++) {
    if (verbose >= 1) {
      fprintf(out,"\nChi^2-Test for marginal distribution [%d]\n",j);
    }
    pval = _unur_test_chi2test(NULL, bm+(j*n_intervals), n_intervals, classmin, verbose, out );
    pval_min = _unur_min(pval_min,pval);
  }
  if (verbose >= 1) {
    fprintf(out,"\nSummary:\n");
    fprintf(out,"  Minimal p-value * number_of_tests = %g:\n\n",pval_min * dim);
  }
free_memory:
  if (X)    free(X);
  if (U)    free(U);
  if (bm)   free(bm);
  if (Linv) free(Linv);
  if (marginals)  free (marginals);
  if (marginal_cdf)  free (marginal_cdf);
  if (Fl)   free(Fl);
  if (Fr)   free(Fr);
  if (Fdelta)   free(Fdelta);
  return pval_min * dim;
#undef idx
#undef DISTR
} 
double
_unur_test_chi2_cvemp ( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
{
#define DISTR     gen->distr->data.cvec
#define idx(i,j)  ((i)*dim+(j))
  int dim;                       
  UNUR_DISTR *distr_normal;      
  UNUR_FUNCT_CONT *cdf;          
  double *X = NULL;              
  double *U = NULL;              
  int *bm = NULL;                
  double pval, pval_min;         
  int i, j;                      
  CHECK_NULL(gen,-1.);
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;
  if( samplesize <= 0 ) samplesize = CHI2_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );
  dim = gen->distr->dim;  
  if (dim < 1) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"distribution dimension < 1 ?");
    return -1.;
  }
  distr_normal = unur_distr_normal( NULL, 0 );
  cdf = distr_normal->data.cont.cdf;
  X   = _unur_xmalloc( dim * sizeof(double));
  U   = _unur_xmalloc( dim * sizeof(double));
  bm  = _unur_xmalloc( dim * n_intervals * sizeof(int)); 
  memset(bm , 0, dim * n_intervals  * sizeof(int));
  for( i=0; i<samplesize; i++ ) {
    _unur_sample_vec(gen, X);
    for (j=0; j<dim; j++)
      U[j] = cdf(X[j], distr_normal);
    for (j=0; j<dim; j++) {
      int iv;
      iv = (int)( n_intervals * U[j] );
      if (iv>=n_intervals) iv = n_intervals-1;
      if (iv < 0) iv = 0;
      bm[j*n_intervals + iv] += 1;
    }
  }
  pval_min = 1.;
  for (j=0; j<dim; j++) {
    if (verbose >= 1) {
      fprintf(out,"\nChi^2-Test for marginal distribution [%d]\n",j);
    }
    pval = _unur_test_chi2test(NULL, bm+(j*n_intervals), n_intervals, classmin, verbose, out );
    pval_min = _unur_min(pval_min,pval);
  }
  if (verbose >= 1) {
    fprintf(out,"\nSummary:\n");
    fprintf(out,"  Minimal p-value * number_of_tests = %g:\n\n",pval_min * dim);
  }
  _unur_distr_free(distr_normal);
  free(X);
  free(U);
  free(bm);
  return pval_min * dim;
#undef idx
#undef DISTR
} 
double
_unur_test_chi2test( double *prob,
         int *observed,
         int len,
         int classmin,
         int verbose,
         FILE *out )
{
  double chi2 = 0.;     
  double df;            
  double pval;          
  double clexpd = 0.;   
  int clobsd = 0;       
  int classes = 0;      
  double probsum = 0.;  
  int samplesize = 0;
  double factor;        
  int i;
  UNUR_DISTR *distr_chisquare = NULL; 
  CHECK_NULL(observed,-1.);
  classmin = (classmin > 0) ? classmin : CHI2_CLASSMIN_DEFAULT;
  for( samplesize=0, i=0; i<len; i++) 
    samplesize += observed[i];
  if (prob != NULL) {
    for( probsum=0., i=0; i<len; i++ )
      probsum += prob[i];
    factor = samplesize/probsum;
  }
  else   
    factor = ((double)samplesize)/len;
  clexpd = 0.;
  clobsd = 0;
  classes = 0;
  for( chi2=0., i=0; i<len; i++ ) {
    clexpd += (prob) ? prob[i]*factor : factor;  
    clobsd += observed[i];                       
    if (clexpd >= classmin || i == len-1) {
      if (clobsd <= 0 && clexpd <= 0.) break;
      chi2 += (clobsd-clexpd)*(clobsd-clexpd)/clexpd;
      if (verbose >= 2)
	fprintf(out,"Class #%d:\tobserved %d\texpected %.2f\n",classes,clobsd,clexpd);
      clexpd = 0.;
      clobsd = 0;
      classes++;
    }
  }
  if (classes < 2) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"too few classes!");
    if (verbose >= 1)
      fprintf(out,"\nCannot run chi^2-Test: too few classes\n");
    return -1.;
  }
  df = (double)(classes-1);
  distr_chisquare = unur_distr_chisquare( &df, 1 );
  if (distr_chisquare->data.cont.cdf) {
    pval = 1. - _unur_cont_CDF( chi2, distr_chisquare );
  }
  else {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF for CHI^2 distribution required");
    pval = -2.;
  }
  _unur_distr_free(distr_chisquare);
  if (verbose >= 1 && pval >= 0.) {
    fprintf(out,"\nResult of chi^2-Test:\n  samplesize = %d\n",samplesize);
    fprintf(out,"  classes    = %d\t (minimum per class = %d)\n", classes, classmin);
    fprintf(out,"  chi2-value = %g\n  p-value    = %g\n\n", chi2, pval);
  }
  return pval;
} 
