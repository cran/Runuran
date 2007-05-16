/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <distributions/unur_stddistr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "cstd.h"
#include "empk.h"
#include "empk_struct.h"
#define EMPK_VARFLAG_VARCOR     0x001u   
#define EMPK_VARFLAG_POSITIVE   0x002u   
#define EMPK_DEBUG_PRINTDATA   0x00000100u
#define EMPK_SET_KERNEL         0x010u    
#define EMPK_SET_KERNGEN        0x020u    
#define EMPK_SET_KERNELVAR      0x001u    
#define EMPK_SET_ALPHA          0x002u    
#define EMPK_SET_BETA           0x004u    
#define EMPK_SET_SMOOTHING      0x008u    
#define GENTYPE "EMPK"         
static struct unur_gen *_unur_empk_init( struct unur_par *par );
static struct unur_gen *_unur_empk_create( struct unur_par *par );
static struct unur_gen *_unur_empk_clone( const struct unur_gen *gen );
static void _unur_empk_free( struct unur_gen *gen);
static double _unur_empk_sample( struct unur_gen *gen );
inline static int _unur_empk_comp_stddev( double *data, int n_data,
					  double *mean, double *stddev);
inline static double _unur_empk_comp_iqrtrange( double *data, int n_data );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_empk_debug_init( const struct unur_par *par, const struct unur_gen *gen );
#endif
#define DISTR_IN  distr->data.cemp      
#define PAR       ((struct unur_empk_par*)par->datap) 
#define GEN       ((struct unur_empk_gen*)gen->datap) 
#define DISTR     gen->distr->data.cemp 
#define SAMPLE    gen->sample.cont           
#define SQU(a) ((a)*(a))
inline static int 
compare_doubles (const void *a, const void *b)
{ 
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}
#define _unur_empk_getSAMPLE(gen)   (_unur_empk_sample)
struct unur_par *
unur_empk_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CEMP,NULL);
  if (DISTR_IN.sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"observed sample"); return NULL; }
  if (DISTR_IN.n_sample < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"number of observed sample"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_empk_par) );
  COOKIE_SET(par,CK_EMPK_PAR);
  par->distr    = distr;          
  PAR->kernvar   = 1.;            
  PAR->alpha     = 0.7763884;     
  PAR->beta      = 1.3637439;     
  PAR->smoothing = 1.;            
  PAR->kerngen   = NULL;          
  PAR->kernel    = NULL;          
  par->method   = UNUR_METH_EMPK; 
  par->variant  = 0u;             
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init     = _unur_empk_init;
  return par;
} 
int 
unur_empk_set_kernel( struct unur_par *par, unsigned kernel)
{
  UNUR_DISTR *kerndist;
  double fpar[4];
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );
  if (par->set & EMPK_SET_KERNEL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"Cannot overwrite kernel");
    return UNUR_ERR_PAR_SET;
  }
  switch (kernel) {
  case UNUR_DISTR_EPANECHNIKOV:
    fpar[0]     = 2.;
    fpar[1]     = 2.;
    fpar[2]     = -1.;
    fpar[3]     = 1.;
    kerndist    = unur_distr_beta( fpar, 4 );
    PAR->kernel  = unur_init( unur_arou_new ( kerndist ) );
    PAR->alpha   = 1.718771928;
    PAR->kernvar = 0.2;
    unur_distr_free( kerndist );
    break;
  case UNUR_DISTR_GAUSSIAN:
    kerndist    = unur_distr_normal( NULL, 0 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 0.7763884;
    PAR->kernvar = 1.;
    unur_distr_free( kerndist );
    break;
  case UNUR_DISTR_BOXCAR:
    fpar[0]     = -1.;
    fpar[1]     = 1.;
    kerndist    = unur_distr_uniform( fpar, 2 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 1.351;
    PAR->kernvar = 1./3.;
    unur_distr_free( kerndist );
    break;
  case UNUR_DISTR_STUDENT:
    fpar[0] = 3.;
    kerndist    = unur_distr_student( fpar, 1 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 0.48263;
    PAR->kernvar = 3.;
    unur_distr_free( kerndist );
    break;
  case UNUR_DISTR_LOGISTIC:
    kerndist    = unur_distr_logistic( NULL, 0 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 0.434;
    PAR->kernvar = 3.289868133696; 
    unur_distr_free( kerndist );
    break;
  default:
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"Unknown kernel. make it manually");
    return UNUR_ERR_PAR_SET;
  }
  if (PAR->kernel == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"Could not initialize kernel generator");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  par->set &= ~EMPK_SET_KERNGEN;   
  par->set |= EMPK_SET_KERNEL | EMPK_SET_ALPHA | EMPK_SET_KERNELVAR;
  return UNUR_SUCCESS;
} 
int
unur_empk_set_kernelgen( struct unur_par *par, const struct unur_gen *kernelgen,
			 double alpha, double kernelvar )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, kernelgen, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );
  if (par->set & EMPK_SET_KERNEL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"Cannot overwrite kernel");
    return UNUR_ERR_PAR_SET;
  }
  if ( (kernelgen->method & UNUR_MASK_TYPE) != UNUR_METH_CONT ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return UNUR_ERR_DISTR_INVALID;
  }
  if (alpha <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"alpha <= 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->kerngen = kernelgen;
  PAR->alpha = alpha;
  par->set |= EMPK_SET_KERNGEN | EMPK_SET_ALPHA;
  PAR->kernvar = kernelvar;
  if (kernelvar > 0.)
    par->set |= EMPK_SET_KERNELVAR;
  else
    par->set &= ~EMPK_SET_KERNELVAR;
  return UNUR_SUCCESS;
} 
int
unur_empk_set_beta( struct unur_par *par, double beta )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );
  if (beta <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"beta <= 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->beta = beta;
  par->set |= EMPK_SET_BETA;
  return UNUR_SUCCESS;
} 
int
unur_empk_set_smoothing( struct unur_par *par, double smoothing )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->smoothing = smoothing;
  par->set |= EMPK_SET_SMOOTHING;
  return UNUR_SUCCESS;
} 
int
unur_empk_chg_smoothing( struct unur_gen *gen, double smoothing )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, EMPK, UNUR_ERR_GEN_INVALID );
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  GEN->bwidth = smoothing * GEN->bwidth_opt;
  GEN->sconst = 1./sqrt(1. + GEN->kernvar * SQU( GEN->bwidth/GEN->stddev_observ ) );
  GEN->smoothing = smoothing;
  gen->set |= EMPK_SET_SMOOTHING;
  return UNUR_SUCCESS;
} 
int
unur_empk_set_varcor( struct unur_par *par, int varcor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );
  par->variant = (varcor) 
    ? (par->variant | EMPK_VARFLAG_VARCOR) 
    : (par->variant & (~EMPK_VARFLAG_VARCOR));
  return UNUR_SUCCESS;
} 
int
unur_empk_chg_varcor( struct unur_gen *gen, int varcor )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, EMPK, UNUR_ERR_GEN_INVALID );
  if (! (gen->set & EMPK_SET_KERNELVAR) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"variance correction disabled");
    return UNUR_ERR_PAR_SET;
  }
  gen->variant = (varcor) 
    ? (gen->variant | EMPK_VARFLAG_VARCOR) 
    : (gen->variant & (~EMPK_VARFLAG_VARCOR));
  return UNUR_SUCCESS;
} 
int
unur_empk_set_positive( struct unur_par *par, int positive )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );
  par->variant = (positive) 
    ? (par->variant | EMPK_VARFLAG_POSITIVE) 
    : (par->variant & (~EMPK_VARFLAG_POSITIVE));
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_empk_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  double iqrtrange;       
  double sigma;           
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_EMPK ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_EMPK_PAR,NULL);
  if ( PAR->kerngen==NULL && PAR->kernel == NULL) {
    if ( unur_empk_set_kernel( par, UNUR_DISTR_GAUSSIAN )!=UNUR_SUCCESS ) {
      _unur_par_free(par); return NULL; 
    }
  }
  gen = _unur_empk_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }
  if( (gen->variant & EMPK_VARFLAG_VARCOR) &&
      !( (gen->set & EMPK_SET_KERNELVAR) && GEN->kernvar > 0. )) {
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"variance correction disabled");
    gen->variant &= ~EMPK_SET_KERNELVAR;
  }
  GEN->kerngen->urng = par->urng;
  GEN->kerngen->debug = par->debug;
  qsort( GEN->observ, (size_t)GEN->n_observ, sizeof(double), compare_doubles);
  _unur_empk_comp_stddev( GEN->observ, GEN->n_observ, &(GEN->mean_observ), &(GEN->stddev_observ) );
  iqrtrange = _unur_empk_comp_iqrtrange( GEN->observ, GEN->n_observ );
  sigma = iqrtrange / 1.34;
  if (GEN->stddev_observ < sigma) sigma = GEN->stddev_observ;
  GEN->bwidth_opt = PAR->alpha * PAR->beta * sigma / exp(0.2 * log((double)GEN->n_observ));
  GEN->bwidth = PAR->smoothing * GEN->bwidth_opt;
  GEN->sconst = 1./sqrt(1. + PAR->kernvar * SQU( GEN->bwidth/GEN->stddev_observ ) );
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_empk_debug_init(par,gen);
#endif
  _unur_par_free(par);
  return gen;
} 
static struct unur_gen *
_unur_empk_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_EMPK_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_empk_gen) );
  COOKIE_SET(gen,CK_EMPK_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_empk_getSAMPLE(gen);
  gen->destroy = _unur_empk_free;
  gen->clone = _unur_empk_clone;
  GEN->observ   = DISTR.sample;          
  GEN->n_observ = DISTR.n_sample;        
  GEN->smoothing = PAR->smoothing;    
  if (PAR->kerngen)
    GEN->kerngen = _unur_gen_clone(PAR->kerngen);
  else
    GEN->kerngen = PAR->kernel;
  GEN->kernvar = PAR->kernvar;
  gen->gen_aux = GEN->kerngen;
  return gen;
} 
struct unur_gen *
_unur_empk_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_empk_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_EMPK_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->observ = clone->distr->data.cemp.sample;   
  CLONE->kerngen = clone->gen_aux;
  return clone;
#undef CLONE
} 
void
_unur_empk_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_EMPK ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_EMPK_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_empk_sample( struct unur_gen *gen )
{ 
  double U,K,X;
  int j;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_EMPK_GEN,INFINITY);
  U = _unur_call_urng(gen->urng) * GEN->n_observ;
  j = (int) (U);
  K = unur_sample_cont( GEN->kerngen );
  if (gen->variant & EMPK_VARFLAG_VARCOR)
    X = GEN->mean_observ + (GEN->observ[j] - GEN->mean_observ + GEN->bwidth * K) * GEN->sconst;
  else
    X = GEN->observ[j] + GEN->bwidth * K;
  if (gen->variant & EMPK_VARFLAG_POSITIVE)
    X = (X<0.) ? -X : X;
  return X;
} 
int
_unur_empk_comp_stddev( double *data, int n_data, double *mean, double *stddev)
{
  double xsqu_sum;   
  double dx;
  int n;
  if (n_data < 2)
    return 0;
  *mean = 0.;
  xsqu_sum = 0.;
  for (n=1; n <= n_data; n++) {
    dx = (data[n-1] - *mean) / n;
    xsqu_sum += n * (n - 1.) * dx * dx;
    *mean += dx;
  }
  *stddev = sqrt( xsqu_sum / (n_data - 1.));
  return UNUR_SUCCESS;
} 
double
_unur_empk_comp_iqrtrange( double *data, int n )
{
  double lowerqrt,upperqrt;  
  int j;
  j = n/2;
  if (j % 2) {
    lowerqrt = data[(j+1)/2-1];
    upperqrt = data[n-(j+1)/2];
  }
  else {
    lowerqrt = (data[j/2-1] + data[j/2+1-1])/2.;
    upperqrt = (data[n-j/2] + data[n-j/2-1])/2.;
  }
  return (upperqrt - lowerqrt);
} 
#ifdef UNUR_ENABLE_LOGGING
static void
_unur_empk_debug_init( const struct unur_par *par, const struct unur_gen *gen )
{
  FILE *log;
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_EMPK_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_EMPK_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = EMPK (EMPirical distribution with Kernel smoothing)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_distr_cemp_debug( gen->distr, gen->genid, (gen->debug & EMPK_DEBUG_PRINTDATA));
  fprintf(log,"%s: sampling routine = _unur_empk_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: smoothing factor = %g",gen->genid, PAR->smoothing);
  _unur_print_if_default(par,EMPK_SET_SMOOTHING); fprintf(log,"\n");
  if (gen->variant & EMPK_VARFLAG_POSITIVE)
    fprintf(log,"%s: positive random variable only; use mirroring \n",gen->genid);
  if (gen->variant & EMPK_VARFLAG_VARCOR)
    fprintf(log,"%s: use variance correction\n",gen->genid);
  else
    fprintf(log,"%s: no variance correction\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: Kernel:\n",gen->genid);
  fprintf(log,"%s:    type = %s  ",gen->genid,GEN->kerngen->distr->name);
  if (gen->set & EMPK_SET_KERNGEN)
    fprintf(log,"[kernel generator set]\n");
  else if (gen->set & EMPK_SET_KERNEL)
    fprintf(log,"[standard kernel]\n");
  else 
    fprintf(log,"[default kernel]\n");
  fprintf(log,"%s:    window width = %g\t(opt = %g)\n",gen->genid, GEN->bwidth, GEN->bwidth_opt);
  fprintf(log,"%s:    alpha = %g",gen->genid, PAR->alpha);
  _unur_print_if_default(par,EMPK_SET_ALPHA); fprintf(log,"\n");
  if (gen->variant & EMPK_VARFLAG_VARCOR) {
    fprintf(log,"%s:    kernel variance = %g",gen->genid, PAR->kernvar);
    _unur_print_if_default(par,EMPK_SET_KERNELVAR); fprintf(log,"\n");
    fprintf(log,"%s:    variance correction factor = %g\n",gen->genid, GEN->sconst);
  }
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: Data:\n",gen->genid);
  fprintf(log,"%s:    beta  = %g",gen->genid, PAR->beta);
  _unur_print_if_default(par,EMPK_SET_BETA); fprintf(log,"\n");
  fprintf(log,"%s:    mean (data) = %g\n",gen->genid, GEN->mean_observ);
  fprintf(log,"%s:    stddev (data) = %g\n",gen->genid, GEN->stddev_observ);
  fprintf(log,"%s:\n",gen->genid);
} 
#endif   
