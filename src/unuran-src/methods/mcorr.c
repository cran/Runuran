/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/matr.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include <utils/unur_fp_source.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "mcorr.h"
#include "mcorr_struct.h"
#define MCORR_DEBUG_REINIT    0x00000010u  
#define MCORR_SET_EIGENVALUES  0x001u   
#define GENTYPE "MCORR"        
static struct unur_gen *_unur_mcorr_init( struct unur_par *par );
static int _unur_mcorr_init_HH( struct unur_gen *gen );
static int _unur_mcorr_init_eigen( struct unur_gen *gen );
static int _unur_mcorr_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_mcorr_create( struct unur_par *par );
static struct unur_gen *_unur_mcorr_clone( const struct unur_gen *gen );
static void _unur_mcorr_free( struct unur_gen *gen);
static int _unur_mcorr_sample_matr_HH( struct unur_gen *gen, double *mat );
static int _unur_mcorr_sample_matr_eigen( struct unur_gen *gen, double *mat );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_mcorr_debug_init( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_mcorr_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.matr      
#define PAR       ((struct unur_mcorr_par*)par->datap) 
#define GEN       ((struct unur_mcorr_gen*)gen->datap) 
#define DISTR     gen->distr->data.matr 
#define SAMPLE    gen->sample.matr      
#define NORMAL    gen->gen_aux        
#define _unur_mcorr_getSAMPLE(gen) \
   ( ((gen)->set & MCORR_SET_EIGENVALUES) \
     ? _unur_mcorr_sample_matr_eigen : _unur_mcorr_sample_matr_HH )
struct unur_par *
unur_mcorr_new( const struct unur_distr *distr )
{
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if ( !(distr->type == UNUR_DISTR_MATR &&
	 distr->id == UNUR_DISTR_MCORRELATION) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_MATR,NULL);
  par = _unur_par_new( sizeof(struct unur_mcorr_par) );
  COOKIE_SET(par,CK_MCORR_PAR);
  par->distr    = distr;      
  par->method   = UNUR_METH_MCORR;    
  par->variant  = 0u;                 
  par->set      = 0u;                 
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  PAR->dim = distr->data.matr.n_rows;
  PAR->eigenvalues = NULL; 
  par->init = _unur_mcorr_init;
  return par;
} 
int
unur_mcorr_set_eigenvalues( UNUR_PAR *par, const double *eigenvalues )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCORR );
  _unur_check_NULL( GENTYPE, eigenvalues, UNUR_ERR_NULL );
  for (i=0; i<PAR->dim; i++)
    if (eigenvalues[i] <= 0.) {
      _unur_error(GENTYPE, UNUR_ERR_PAR_SET,"eigenvalue <= 0");
      return UNUR_ERR_PAR_SET;
    }
  PAR->eigenvalues = eigenvalues;
  par->set |= MCORR_SET_EIGENVALUES;
  return UNUR_SUCCESS;
} 
int
unur_mcorr_chg_eigenvalues( UNUR_GEN *gen, const double *eigenvalues )
{
  int i;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, MCORR, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( GENTYPE, eigenvalues, UNUR_ERR_NULL );
  for (i=0; i<GEN->dim; i++)
    if (eigenvalues[i] <= 0.) {
      _unur_error(GENTYPE, UNUR_ERR_PAR_SET,"eigenvalue <= 0");
      return UNUR_ERR_PAR_SET;
    }
  if (GEN->eigenvalues == NULL)
    GEN->eigenvalues = _unur_xmalloc(GEN->dim * sizeof(double));
  memcpy(GEN->eigenvalues, eigenvalues, GEN->dim * sizeof(double));
  gen->set |= MCORR_SET_EIGENVALUES;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_mcorr_init( struct unur_par *par )
{
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_MCORR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MCORR_PAR,NULL);
  gen = _unur_mcorr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL; 
  if (gen->set && MCORR_SET_EIGENVALUES) {
    if (_unur_mcorr_init_eigen(gen) != UNUR_SUCCESS) {
      _unur_mcorr_free(gen); return NULL;
    }
  }
  else {
    if (_unur_mcorr_init_HH(gen) != UNUR_SUCCESS) {
      _unur_mcorr_free(gen); return NULL;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_mcorr_debug_init(gen);
#endif
  return gen;
} 
int
_unur_mcorr_init_HH( struct unur_gen *gen )
{
  if (NORMAL==NULL) {
    struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
    struct unur_par   *normalpar = unur_arou_new( normaldistr );
    unur_arou_set_usedars( normalpar, TRUE );
    NORMAL = unur_init( normalpar );
    _unur_distr_free( normaldistr );
    if (NORMAL == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"Cannot create aux Gaussian generator");
      return UNUR_FAILURE;
    }
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
  }
  return UNUR_SUCCESS;
} 
int
_unur_mcorr_init_eigen( struct unur_gen *gen )
{
  int i;
  double sum_eigenvalues = 0.;
  GEN->M = _unur_xrealloc(GEN->M, (5*GEN->dim + 2*GEN->dim*GEN->dim)*sizeof(double));
  for (i=0; i<GEN->dim; i++) {
    if (GEN->eigenvalues[i] <= 0.) {
      _unur_error(GENTYPE, UNUR_ERR_SHOULD_NOT_HAPPEN,"eigenvalue <= 0");
      return UNUR_FAILURE;
    }
    sum_eigenvalues += GEN->eigenvalues[i];
  }
  if (!_unur_FP_equal(sum_eigenvalues, (double) GEN->dim))
    _unur_warning(GENTYPE, UNUR_ERR_GENERIC,"scaling sum(eigenvalues) -> dim");
  for (i=0; i<GEN->dim; i++)
    GEN->eigenvalues[i] *= GEN->dim / sum_eigenvalues;
  return UNUR_SUCCESS;
} 
int
_unur_mcorr_reinit( struct unur_gen *gen )
{
  SAMPLE = _unur_mcorr_getSAMPLE(gen);
  if (gen->set && MCORR_SET_EIGENVALUES)
    return _unur_mcorr_init_eigen(gen);
  else
    return _unur_mcorr_init_HH(gen);
} 
struct unur_gen *
_unur_mcorr_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MCORR_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_mcorr_gen) );
  COOKIE_SET(gen,CK_MCORR_GEN);
  GEN->dim = DISTR.n_rows;
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_mcorr_getSAMPLE(gen);
  gen->destroy = _unur_mcorr_free;
  gen->clone = _unur_mcorr_clone;
  gen->reinit = _unur_mcorr_reinit;
  GEN->M = NULL;
  GEN->H = NULL;
  GEN->eigenvalues = NULL;
  if (gen->set && MCORR_SET_EIGENVALUES) {
    GEN->eigenvalues = _unur_xmalloc(GEN->dim * sizeof(double));
    memcpy(GEN->eigenvalues, PAR->eigenvalues, GEN->dim * sizeof(double));
  }
  if (gen->set && MCORR_SET_EIGENVALUES) {
    GEN->M = _unur_xmalloc((5*GEN->dim + 2*GEN->dim*GEN->dim) * sizeof(double));
  }
  else {
    GEN->H = _unur_xmalloc(GEN->dim * GEN->dim * sizeof(double));
  }
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_mcorr_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_mcorr_clone( const struct unur_gen *gen )
{
#define CLONE  ((struct unur_mcorr_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MCORR_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  if (GEN->M) 
    CLONE->M = _unur_xmalloc((5*GEN->dim + 2*GEN->dim*GEN->dim) * sizeof(double));
  if (GEN->H)
    CLONE->H = _unur_xmalloc(GEN->dim * GEN->dim * sizeof(double));
  if (GEN->eigenvalues) {
    CLONE->eigenvalues = _unur_xmalloc(GEN->dim * sizeof(double));
    memcpy(CLONE->eigenvalues, GEN->eigenvalues, GEN->dim * sizeof(double));
  }
  return clone;
#undef CLONE
} 
void
_unur_mcorr_free( struct unur_gen *gen )
{
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_MCORR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->eigenvalues) free(GEN->eigenvalues);
  if (GEN->H)           free(GEN->H);
  if (GEN->M)           free(GEN->M);
  _unur_generic_free(gen);
} 
int
_unur_mcorr_sample_matr_HH( struct unur_gen *gen, double *mat )
{
#define idx(a,b) ((a)*(GEN->dim)+(b))
  int i,j,k;
  double sum, norm, x;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MCORR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(mat,UNUR_ERR_NULL);
  for (i=0; i<GEN->dim; i++) {
    sum=0.;
    for (j=0; j<GEN->dim; j++) {
      x = _unur_sample_cont(NORMAL);
      GEN->H[idx(i,j)] = x;
      sum += x * x;
    }
    norm = sqrt(sum);
    for (j=0; j<GEN->dim; j++) GEN->H[idx(i,j)] /= norm;
  }
  for (i=0; i<GEN->dim; i++)
    for (j=0; j<GEN->dim; j++) {
      if (j<i)
	mat[idx(i,j)] = mat[idx(j,i)];
      else if(j==i)
	mat[idx(i,j)] = 1.;
      else {
	sum=0.;
	for (k=0; k<GEN->dim; k++)
	  sum += GEN->H[idx(i,k)]*GEN->H[idx(j,k)];
	mat[idx(i,j)] = sum;
      }
    }
  return UNUR_SUCCESS;
#undef idx
} 
int
_unur_mcorr_sample_matr_eigen( struct unur_gen *gen, double *mat )
{
#define idx(a,b) ((a)*dim+(b))
  int i,j,k, dim;
  double *E, *P;
  double *x, *y, *z, *w, *r; 
  double a, b, c, e, e2;
  int s; 
  CHECK_NULL(gen, UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MCORR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(mat, UNUR_ERR_NULL);
  dim = GEN->dim; 
  if (dim<1) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"dimension < 1");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  x = GEN->M + (0*dim);
  y = GEN->M + (1*dim);
  z = GEN->M + (2*dim);
  w = GEN->M + (3*dim);
  r = GEN->M + (4*dim);
  E = GEN->M + (5*dim);
  P = GEN->M + (5*dim+dim*dim);
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++)
      E[idx(i,j)] = (i==j) ? 1 : 0;
  for (k=0; k<dim-1; k++) {
    for (i=0; i<dim; i++) w[i] = _unur_call_urng(gen->urng);
    for (i=0; i<dim; i++) {
      x[i]=0;
      for (j=0; j<dim; j++) {
        x[i] += E[idx(i,j)]*w[j];
      }
    }
    a=0;
    for (i=0; i<dim; i++)
      a += (1-GEN->eigenvalues[i])*x[i]*x[i];
    if (fabs(a)<DBL_EPSILON) {
      for (i=0; i<dim; i++) {
	for (j=0; j<dim; j++) {
	  mat[idx(i,j)] = (i==j) ? 1: 0;
	}}
      _unur_warning(gen->genid, UNUR_ERR_GEN_CONDITION,"all eigenvalues are ~1 -> identity matrix");
      return UNUR_ERR_GEN_CONDITION;
    }
    do {
      for (i=0; i<dim; i++) z[i] = _unur_call_urng(gen->urng);
      for (i=0; i<dim; i++) {
        y[i]=0;
        for (j=0; j<dim; j++) {
          y[i] += E[idx(i,j)]*z[j];
        }
      }
      b=0; c=0;
      for (i=0; i<dim; i++) {
        b += (1-GEN->eigenvalues[i])*x[i]*y[i];
        c += (1-GEN->eigenvalues[i])*y[i]*y[i];
      }
      e2 = b*b - a*c;
    } while (e2<0);
    e=sqrt(e2);
    s = ( _unur_call_urng(gen->urng) >.5) ? 1: -1 ;
    for (i=0; i<dim; i++) r[i] = x[i]*(b+s*e)/a - y[i];
    s = ( _unur_call_urng(gen->urng) >.5) ? 1: -1 ;
    _unur_vector_normalize(dim, r);
    for (i=0; i<dim; i++) P[idx(k,i)] = s * r[i];
    for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        E[idx(i,j)] -= r[i]*r[j];
      }
    }
  } 
  for (i=0; i<dim; i++) w[i] = _unur_call_urng(gen->urng);
  for (i=0; i<dim; i++) {
    x[i]=0;
    for (j=0; j<dim; j++) {
      x[i] += E[idx(i,j)]*w[j];
    }
  }
  _unur_vector_normalize(dim, x);
  for (i=0; i<dim; i++) {
    P[idx(dim-1,i)] = x[i];
  }
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      mat[idx(i,j)] = 0;
      for (k=0; k<dim; k++) {
        mat[idx(i,j)] += P[idx(i,k)] * GEN->eigenvalues[k] * P[idx(j,k)];
      }
    }
  }
  for (i=0; i<dim; i++) {
    for (j=(i+1); j<dim; j++) {
      mat[idx(i,j)] = (mat[idx(i,j)]+mat[idx(j,i)])/2.;
      mat[idx(j,i)] = mat[idx(i,j)];
    }
  }  
  return UNUR_SUCCESS;
#undef idx
} 
#ifdef UNUR_ENABLE_LOGGING
static void
_unur_mcorr_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = random matrix\n",gen->genid);
  fprintf(LOG,"%s: method  = MCORR (Matrix - CORRELATION matrix)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_matr_debug( gen->distr, gen->genid );
  if (gen->set && MCORR_SET_EIGENVALUES)
    _unur_matrix_print_vector( GEN->dim, GEN->eigenvalues, "eigenvalues =", LOG, gen->genid, "\t   ");
  if (gen->set && MCORR_SET_EIGENVALUES)
    fprintf(LOG,"%s: sampling routine = _unur_mcorr_sample_matr_eigen()\n",gen->genid);
  else
    fprintf(LOG,"%s: sampling routine = _unur_mcorr_sample_matr_HH()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_mcorr_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d x %d   (= %d)\n", 
		      DISTR.n_rows, DISTR.n_cols, gen->distr->dim);
  if (gen->set && MCORR_SET_EIGENVALUES) {
    _unur_string_append(info,"   eigenvalues = ");
    _unur_distr_info_vector( gen, GEN->eigenvalues, GEN->dim);
    _unur_string_append(info,"\n");
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: MCORR (Random CORRelation matrix)\n");
  if (gen->set && MCORR_SET_EIGENVALUES)
    _unur_string_append(info,"   generate correlation matrix with given eigenvalues\n");
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: \n");
    if (gen->set && MCORR_SET_EIGENVALUES) {
      _unur_string_append(info,"   eigenvalues = ");
      _unur_distr_info_vector( gen, GEN->eigenvalues, GEN->dim);
      _unur_string_append(info,"\n");
    }
    _unur_string_append(info,"\n");
  }
} 
#endif   
