/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "mvstd.h"
#include "vempk.h"
#include "vempk_struct.h"
#define VEMPK_VARFLAG_VARCOR    0x001u   
#define VEMPK_DEBUG_PRINTDATA   0x00000100u
#define VEMPK_SET_SMOOTHING      0x008u    
#define GENTYPE "VEMPK"         
static struct unur_gen *_unur_vempk_init( struct unur_par *par );
static struct unur_gen *_unur_vempk_create( struct unur_par *par );
static struct unur_gen *_unur_vempk_clone( const struct unur_gen *gen );
static void _unur_vempk_free( struct unur_gen *gen);
static int _unur_vempk_sample_cvec( struct unur_gen *gen, double *result );
static int compute_mean_covar( double *data, int n_data, int dim, double *xbar, double *S );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_vempk_debug_init( const struct unur_par *par, const struct unur_gen *gen );
#endif
#define DISTR_IN  distr->data.cvemp      
#define PAR       ((struct unur_vempk_par*)par->datap) 
#define GEN       ((struct unur_vempk_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvemp 
#define SAMPLE    gen->sample.cvec           
#define _unur_vempk_getSAMPLE(gen)  ( _unur_vempk_sample_cvec )
struct unur_par *
unur_vempk_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CVEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEMP,NULL);
  if (DISTR_IN.sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"observed sample"); return NULL; }
  if (DISTR_IN.n_sample < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"size of observed sample"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_vempk_par) );
  COOKIE_SET(par,CK_VEMPK_PAR);
  par->distr    = distr;         
  PAR->smoothing = 1.;            
  par->method   = UNUR_METH_VEMPK; 
  par->variant  = 0u;              
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init     = _unur_vempk_init;
  return par;
} 
int
unur_vempk_set_smoothing( struct unur_par *par, double smoothing )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VEMPK );
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->smoothing = smoothing;
  par->set |= VEMPK_SET_SMOOTHING;
  return UNUR_SUCCESS;
} 
int
unur_vempk_chg_smoothing( struct unur_gen *gen, double smoothing )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VEMPK, UNUR_ERR_GEN_INVALID );
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  GEN->smoothing = smoothing;
  GEN->hact = GEN->hopt * GEN->smoothing;
  GEN->corfac = 1./sqrt( 1. + GEN->hact * GEN->hact);
  gen->set |= VEMPK_SET_SMOOTHING;
  return UNUR_SUCCESS;
} 
int
unur_vempk_set_varcor( struct unur_par *par, int varcor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VEMPK );
  par->variant = (varcor) 
    ? (par->variant | VEMPK_VARFLAG_VARCOR) 
    : (par->variant & (~VEMPK_VARFLAG_VARCOR));
  return UNUR_SUCCESS;
} 
int
unur_vempk_chg_varcor( struct unur_gen *gen, int varcor )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VEMPK, UNUR_ERR_GEN_INVALID );
  gen->variant = (varcor) 
    ? (gen->variant | VEMPK_VARFLAG_VARCOR) 
    : (gen->variant & (~VEMPK_VARFLAG_VARCOR));
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_vempk_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  double *S;                  
  UNUR_DISTR *kernel_distr;   
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_VEMPK ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VEMPK_PAR,NULL);
  gen = _unur_vempk_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }
  GEN->xbar = _unur_xmalloc( GEN->dim * sizeof(double) );
  S  = _unur_xmalloc( GEN->dim * GEN->dim * sizeof(double) );
  compute_mean_covar( DISTR.sample, DISTR.n_sample, GEN->dim, GEN->xbar, S );
  kernel_distr = unur_distr_multinormal( GEN->dim, NULL, S );
  GEN->kerngen = unur_init( unur_mvstd_new( kernel_distr ) );
  if (GEN->kerngen==NULL) {
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_par_free (par); free (S); _unur_vempk_free(gen);
    return NULL;
  }
  GEN->kerngen->urng = par->urng;
  GEN->kerngen->debug = par->debug;
  gen->gen_aux = GEN->kerngen;
  GEN->hopt = (exp((1./(GEN->dim+4.))*log(4./(GEN->dim+2.)))*
	      exp(-1./(GEN->dim+4.) *log((double)GEN->n_observ)));
  GEN->hact = GEN->hopt * GEN->smoothing;
  GEN->corfac = 1./sqrt(1. + GEN->hact * GEN->hact);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_vempk_debug_init(par,gen);
#endif
  _unur_par_free(par);
  free(S);
  unur_distr_free(kernel_distr);
  return gen;
} 
static struct unur_gen *
_unur_vempk_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VEMPK_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_vempk_gen) );
  COOKIE_SET(gen,CK_VEMPK_GEN);
  GEN->dim = gen->distr->dim; 
  GEN->observ   = DISTR.sample;          
  GEN->n_observ = DISTR.n_sample;        
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_vempk_getSAMPLE(gen);
  gen->destroy = _unur_vempk_free;
  gen->clone = _unur_vempk_clone;
  GEN->smoothing = PAR->smoothing;    
  GEN->kerngen = NULL;               
  GEN->xbar = NULL;                  
  return gen;
} 
struct unur_gen *
_unur_vempk_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_vempk_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VEMPK_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->observ = clone->distr->data.cvemp.sample;   
  if (GEN->xbar) {
    CLONE->xbar = _unur_xmalloc( GEN->dim * sizeof(double) );
    memcpy( CLONE->xbar, GEN->xbar, GEN->dim * sizeof(double) );
  }
  CLONE->kerngen = clone->gen_aux;
  return clone;
#undef CLONE
} 
void
_unur_vempk_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_VEMPK ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VEMPK_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->xbar)   free( GEN->xbar );
  _unur_generic_free(gen);
} 
int
_unur_vempk_sample_cvec( struct unur_gen *gen, double *result )
{ 
  #define idx(a,b) (a*GEN->dim+b)
  double U;
  int j,k;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_VEMPK_GEN,UNUR_ERR_COOKIE);
  U = _unur_call_urng(gen->urng) * GEN->n_observ;
  j = (int) (U);
  unur_sample_vec( GEN->kerngen, result );
  if (gen->variant & VEMPK_VARFLAG_VARCOR)
    for (k=0; k<GEN->dim; k++) 
      result[k] = GEN->xbar[k] + (GEN->observ[idx(j,k)] - GEN->xbar[k] + result[k]*GEN->hact) * GEN->corfac;
  else
    for (k=0; k<GEN->dim; k++) 
      result[k] = GEN->hact * result[k] + GEN->observ[idx(j,k)];
  return UNUR_SUCCESS;
#undef idx
} 
int
compute_mean_covar( double *data, int n_data, int dim, 
		    double *xbar, double *S ) 
{
#define idx(a,b) (a*dim+b)
  int i,j,k;
  double *x;
  x = malloc(dim*sizeof(double));
  for(j=0; j<dim; j++) {
    xbar[j] = 0.;
    for(k=0; k<dim; k++)
      S[idx(j,k)] = 0.;
  }
  for (i=0; i<n_data; i++)
    for(j=0; j<dim; j++)
      xbar[j] += data[idx(i,j)];
  for(j=0; j<dim; j++)
    xbar[j] /= n_data;
  for (i=0; i<n_data; i++) {
    for(j=0; j<dim; j++)
      x[j] = data[idx(i,j)] - xbar[j];
    for(j=0; j<dim; j++)
      for(k=0; k<=j; k++) 
	S[idx(j,k)] += x[j] * x[k];
  }
  for (j=dim-1; j>=0; j--)
    for (k=0; k<=j; k++) {
      S[idx(j,k)] /= (n_data-1);
      if (k!=j)
	S[idx(k,j)] = S[idx(j,k)];   
    }
  free(x);
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
static void
_unur_vempk_debug_init( const struct unur_par *par, const struct unur_gen *gen )
{
  FILE *log;
  int i;
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_VEMPK_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VEMPK_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = VEMPK ((Vector) EMPirical distribution with Kernel smoothing)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_distr_cvemp_debug( gen->distr, gen->genid, (gen->debug & VEMPK_DEBUG_PRINTDATA));
  fprintf(log,"%s:\tmean vector =\n",gen->genid);
  fprintf(log,"%s:\t   ( %g",gen->genid,GEN->xbar[0]);
  for (i=1; i<GEN->dim; i++) 
    fprintf(log,", %g",GEN->xbar[i]);
  fprintf(log,")\n%s:\n",gen->genid);
  fprintf(log,"%s:\tcovariance matrix = [see %s]\n",gen->genid, GEN->kerngen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: sampling routine = _unur_vempk_sample_cvec()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: smoothing factor = %g",gen->genid, PAR->smoothing);
  _unur_print_if_default(par,VEMPK_SET_SMOOTHING); fprintf(log,"\n");
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: bandwith hopt = %g\n",gen->genid, GEN->hopt);
  fprintf(log,"%s: (used)   hact = %g\n",gen->genid, GEN->hact);
  fprintf(log,"%s:\n",gen->genid);
  if (gen->variant & VEMPK_VARFLAG_VARCOR) {
    fprintf(log,"%s: use variance correction\n",gen->genid);
    fprintf(log,"%s:\tcorrection factor = %g\n",gen->genid, GEN->corfac);
  }
  else
    fprintf(log,"%s: no variance correction\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
} 
#endif   
