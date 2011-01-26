/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distr/cont.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "cstd.h"
#include "hinv.h"
#include "ninv.h"
#include "pinv.h"
#include "mvstd.h"
#include "norta.h"
#include "norta_struct.h"
#define UNUR_NORTA_MIN_EIGENVALUE  (1.e-10)
#define NORTA_DEBUG_SIGMA_Y     0x00000010u   
#define GENTYPE "NORTA"          
static struct unur_gen *_unur_norta_init( struct unur_par *par );
static struct unur_gen *_unur_norta_create( struct unur_par *par );
static struct unur_gen *_unur_norta_clone( const struct unur_gen *gen );
static void _unur_norta_free( struct unur_gen *gen);
static int _unur_norta_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_norta_nortu_setup( struct unur_gen *gen );
static int _unur_norta_make_correlationmatrix( int dim, double *M);
static struct unur_gen *_unur_norta_make_marginalgen( const struct unur_gen *gen,
						      const struct unur_distr *marginal );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_norta_debug_init( const struct unur_gen *gen );
static void _unur_norta_debug_sigma_y( const struct unur_gen *gen, 
				       const double *sigma_y, 
				       const char *comment );
static void _unur_norta_debug_eigensystem( const struct unur_gen *gen,
					   const double *eigenvalues,
					   const double *eigenvectors );
static void _unur_norta_debug_nmgenerator( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_norta_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cvec      
#define PAR       ((struct unur_norta_par*)par->datap) 
#define GEN       ((struct unur_norta_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec 
#define SAMPLE    gen->sample.cvec           
#define MNORMAL   gen->gen_aux          
#define _unur_norta_getSAMPLE(gen)   (_unur_norta_sample_cvec)
struct unur_par *
unur_norta_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);
  if (!(distr->set & UNUR_DISTR_SET_RANKCORR)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"rank correlation matrix");
    return NULL; }
  if (!(distr->set & UNUR_DISTR_SET_MARGINAL)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"marginals");
    return NULL; }
  par = _unur_par_new( sizeof(struct unur_norta_par) );
  COOKIE_SET(par,CK_NORTA_PAR);
  par->distr    = distr;            
  par->method   = UNUR_METH_NORTA ;   
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_norta_init;
  return par;
} 
struct unur_gen *
_unur_norta_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_NORTA ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NORTA_PAR,NULL);
  gen = _unur_norta_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if ( gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED ) {
    if ( DISTR.domainrect == NULL ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot handle non-rectangular domain");
      _unur_norta_free(gen); return NULL;
    }
    else { 
      if (_unur_distr_cvec_marginals_are_equal(DISTR.marginals, GEN->dim)) {
	if ( _unur_distr_cvec_duplicate_firstmarginal(gen->distr) != UNUR_SUCCESS ) {
	  _unur_norta_free(gen); return NULL;
	}
      }
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_norta_debug_init(gen);
#endif
  if (_unur_norta_nortu_setup(gen) != UNUR_SUCCESS) {
    _unur_norta_free(gen); return NULL;
  }
  GEN->normaldistr = unur_distr_normal(NULL,0);
  if (gen->distr->id != UNUR_DISTR_COPULA) {
    if (_unur_distr_cvec_marginals_are_equal(DISTR.marginals, GEN->dim)) {
      struct unur_gen *marginalgen = _unur_norta_make_marginalgen( gen, DISTR.marginals[0] );
      if (marginalgen)
	GEN->marginalgen_list = _unur_gen_list_set(marginalgen,GEN->dim);
    }
    else {
      int i,j;
      int failed = FALSE;
      struct unur_gen **marginalgens = _unur_xmalloc( GEN->dim * sizeof(struct unur_gen*) );
      if ( gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED ) {
	for (i=0; i<GEN->dim && !failed; i++) {
	  if ( (unur_distr_cont_set_domain(DISTR.marginals[i],
					   DISTR.domainrect[2*i], DISTR.domainrect[2*i+1]))
	       != UNUR_SUCCESS) {
	    failed = TRUE; break;
	  }
	}
      }
      for (i=0; i<GEN->dim && !failed; i++) {
	marginalgens[i] = _unur_norta_make_marginalgen( gen, DISTR.marginals[i] );
	if (marginalgens[i]==NULL) {
	  failed=TRUE; break;
	}
      }
      if (failed) {
	for (j=0; j<i; j++) _unur_free(marginalgens[j]);
	free (marginalgens);
      }
      else
	GEN->marginalgen_list = marginalgens;
    }
    if (GEN->marginalgen_list == NULL) {
      _unur_error(gen->genid,UNUR_ERR_GENERIC,"init of marginal generators failed");
      _unur_norta_free(gen);
      return NULL;
    }
  }
  return gen;
} 
static struct unur_gen *
_unur_norta_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NORTA_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_norta_gen) );
  COOKIE_SET(gen,CK_NORTA_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_norta_getSAMPLE(gen);
  gen->destroy = _unur_norta_free;
  gen->clone = _unur_norta_clone;
  GEN->dim = gen->distr->dim;
  GEN->copula = _unur_xmalloc(sizeof(double)*GEN->dim);
  MNORMAL = NULL;
  GEN->normaldistr = NULL;
  GEN->marginalgen_list = NULL;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_norta_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_norta_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_norta_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NORTA_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->copula = _unur_xmalloc(sizeof(double)*GEN->dim);
  CLONE->normaldistr = _unur_distr_clone(GEN->normaldistr);
  if (GEN->marginalgen_list)
    CLONE->marginalgen_list = _unur_gen_list_clone( GEN->marginalgen_list, GEN->dim );
  return clone;
#undef CLONE
} 
void
_unur_norta_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_NORTA ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);
  if (GEN->copula) free (GEN->copula);
  if (GEN->normaldistr) _unur_distr_free (GEN->normaldistr);
  if (GEN->marginalgen_list)
    _unur_gen_list_free( GEN->marginalgen_list, GEN->dim);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
int
_unur_norta_sample_cvec( struct unur_gen *gen, double *vec )
{
#define idx(a,b) ((a)*GEN->dim+(b))
  int j;
  double *u;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_NORTA_GEN,UNUR_ERR_COOKIE);
  u = GEN->copula;
  _unur_sample_vec(MNORMAL,u);
  for (j=0; j<GEN->dim; j++)
    vec[j] = unur_distr_cont_eval_cdf( u[j], GEN->normaldistr );
  if (gen->distr->id == UNUR_DISTR_COPULA)
    return UNUR_SUCCESS;
  for (j=0; j<GEN->dim; j++) {
    vec[j] = unur_quantile(GEN->marginalgen_list[j], vec[j]);
  }
  return UNUR_SUCCESS;
#undef idx
} 
int
_unur_norta_nortu_setup( struct unur_gen *gen )
{
#define idx(a,b) ((a)*dim+(b))
  int dim = GEN->dim;    
  double *sigma_y;      
  double *eigenvalues;  
  double *eigenvectors; 
  int eigenvalues_positive; 
  struct unur_distr *mn_distr;  
  struct unur_gen   *mn_gen;    
  int i,j;
  sigma_y = _unur_xmalloc(dim * dim * sizeof(double));
  for(i=0; i<dim; i++) {
    for(j=0; j<i; j++)
      sigma_y[idx(i,j)] = sigma_y[idx(j,i)];
    sigma_y[idx(i,i)] = 1.;
    for(j=i+1; j<dim; j++)
      sigma_y[idx(i,j)] = 2.*sin(DISTR.rankcorr[idx(i,j)]*(M_PI/6.));  
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
    _unur_norta_debug_sigma_y( gen, sigma_y, "NORTU setup:" );
#endif
  eigenvalues = _unur_xmalloc(dim * sizeof(double));
  eigenvectors = _unur_xmalloc(dim * dim * sizeof(double));
  if (_unur_matrix_eigensystem(dim, sigma_y, eigenvalues, eigenvectors) != UNUR_SUCCESS) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"cannot compute eigenvalues for given sigma_y");
    free(sigma_y); free(eigenvalues); free(eigenvectors);
    return UNUR_ERR_GEN_DATA;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
    _unur_norta_debug_eigensystem( gen, eigenvalues, eigenvectors );
#endif
  eigenvalues_positive = TRUE;
  for(i=0; i<dim; i++)
    if(eigenvalues[i] < UNUR_NORTA_MIN_EIGENVALUE) {
      eigenvalues[i] = UNUR_NORTA_MIN_EIGENVALUE;
      eigenvalues_positive = FALSE;
    }
  if (!eigenvalues_positive) {
    _unur_matrix_transform_diagonal(dim,eigenvectors,eigenvalues,sigma_y);
    _unur_norta_make_correlationmatrix(dim,sigma_y);
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,
		  "sigma_y not positive definite -> corrected matrix");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
      _unur_norta_debug_sigma_y( gen, sigma_y, "\tEigenvalues < 0 --> correction required" );
#endif
  }
  free(eigenvalues);
  free(eigenvectors);
  mn_distr = unur_distr_multinormal(dim, NULL, sigma_y);
  mn_gen = NULL;
  if (mn_distr) {
    mn_gen = unur_init(unur_mvstd_new(mn_distr));
    _unur_distr_free(mn_distr);
  }
  if (mn_gen == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"(corrected) sigma_y not positive definit");
    free(sigma_y);
    return UNUR_ERR_GEN_DATA;
  }
  MNORMAL = mn_gen;
  MNORMAL->urng = gen->urng;
  MNORMAL->debug = gen->debug;
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
      _unur_norta_debug_nmgenerator( gen );
#endif
  free(sigma_y);
  return UNUR_SUCCESS;
#undef idx
} 
int
_unur_norta_make_correlationmatrix( int dim, double *M)
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;
  for (i=0; i<dim; i++)
    M[idx(i,i)] = sqrt(M[idx(i,i)]);
  for (i=0; i<dim; i++)
    for (j=i; j<dim; j++)
      if(i!=j) {
	M[idx(i,j)] /= M[idx(i,i)] * M[idx(j,j)];
	M[idx(j,i)] = M[idx(i,j)];
      }
  for (i=0; i<dim; i++) 
    M[idx(i,i)] = 1.;
  return UNUR_SUCCESS;
#undef idx
} 
struct unur_gen *
_unur_norta_make_marginalgen( const struct unur_gen *gen,
			      const struct unur_distr *marginal )
{
  struct unur_gen *marginalgen;
  struct unur_par *par;
  CHECK_NULL(gen,NULL);      COOKIE_CHECK(gen,CK_NORTA_GEN,NULL);
  CHECK_NULL(marginal,NULL);
  if (marginal->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(marginal,CK_DISTR_CONT,NULL);
  do {
    par = unur_pinv_new( marginal );
    if ( (marginalgen = _unur_init(par)) != NULL )
      break;
    par = unur_cstd_new( marginal );
    if (unur_cstd_set_variant( par, UNUR_STDGEN_INVERSION)==UNUR_SUCCESS) {
      marginalgen = _unur_init(par);
      break;
    }
    else {
      _unur_par_free(par);
    }
    par = unur_hinv_new( marginal );
    if ( (marginalgen = _unur_init(par)) != NULL )
      break;
    par = unur_ninv_new( marginal );
    unur_ninv_set_table( par, 100 );
    if ( (marginalgen = _unur_init(par)) != NULL ) 
      break;
    _unur_error(gen->genid,UNUR_ERR_DISTR_REQUIRED,"data for (numerical) inversion of marginal missing");
    return NULL;
  } while (1);
  marginalgen->debug = gen->debug;
  return marginalgen;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_norta_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = NORTA (Vector Matrix Transformation)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cvec_debug( gen->distr, gen->genid );
} 
void
_unur_norta_debug_sigma_y( const struct unur_gen *gen, 
			   const double *sigma_y, 
			   const char *comment )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: %s\n",gen->genid,comment);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_matrix_print_matrix( GEN->dim, sigma_y, "\tsigma_y =", 
			     LOG, gen->genid, "\t   ");
} 
void
_unur_norta_debug_eigensystem( const struct unur_gen *gen,
			       const double *eigenvalues,
			       const double *eigenvectors )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  _unur_matrix_print_vector( GEN->dim, eigenvalues, 
			     "\teigenvalues of sigma_y =", 
			     LOG, gen->genid, "\t   ");
  _unur_matrix_print_matrix( GEN->dim, eigenvectors, 
			     "\teigenvectors of sigma_y [rows] =", 
			     LOG, gen->genid, "\t   ");
} 
void
_unur_norta_debug_nmgenerator( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: generator for multinormal auxiliary distribution = %s\n", gen->genid,
	  MNORMAL->genid );
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_norta_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int i;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = MARGINAL distributions\n");
  _unur_string_append(info,"   marginals =");
  for (i=0; i<distr->dim; i++)
    _unur_string_append(info," %s", distr->data.cvec.marginals[i]->name);
  _unur_string_append(info,"\n\n");
  _unur_string_append(info,"method: NORTA (NORmal To Anything)\n");
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
