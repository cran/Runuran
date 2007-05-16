/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "hri.h"
#include "hri_struct.h"
#define HRI_VARFLAG_VERIFY     0x01u    
#define HRB_DEBUG_REINIT    0x00000010u   
#define HRI_DEBUG_SAMPLE    0x01000000u    
#define HRI_SET_P0             0x001u
#define GENTYPE "HRI"          
static struct unur_gen *_unur_hri_init( struct unur_par *par );
static int _unur_hri_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_hri_create( struct unur_par *par );
static int _unur_hri_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_hri_clone( const struct unur_gen *gen );
static void _unur_hri_free( struct unur_gen *gen );
static double _unur_hri_sample( struct unur_gen *gen );
static double _unur_hri_sample_check( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_hri_debug_init( const struct unur_gen *gen );
static void _unur_hri_debug_sample( const struct unur_gen *gen,
				    double x, double p1, int i0, int i1 );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_hri_par*)par->datap) 
#define GEN       ((struct unur_hri_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define SAMPLE    gen->sample.cont      
#define HR(x)     _unur_cont_HR((x),(gen->distr))   
#define _unur_hri_getSAMPLE(gen) \
   ( ((gen)->variant & HRI_VARFLAG_VERIFY) \
     ? _unur_hri_sample_check : _unur_hri_sample )
struct unur_par *
unur_hri_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.hr == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"HR"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_hri_par) );
  COOKIE_SET(par,CK_HRI_PAR);
  par->distr   = distr;         
  PAR->p0        = 1.;                     
  par->method   = UNUR_METH_HRI;           
  par->variant  = 0u;                      
  par->set      = 0u;                      
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_hri_init;
  return par;
} 
int
unur_hri_set_p0( struct unur_par *par, double p0 )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HRI );
  if (p0 <= par->distr->data.cont.domain[0]) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"p0 <= left boundary");
    return UNUR_ERR_PAR_SET;
  }
  PAR->p0 = p0;
  par->set |= HRI_SET_P0;
  return UNUR_SUCCESS;
} 
int
unur_hri_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HRI );
  par->variant = (verify) ? (par->variant | HRI_VARFLAG_VERIFY) : (par->variant & (~HRI_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_hri_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, HRI, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  gen->variant = (verify) 
    ? (gen->variant | HRI_VARFLAG_VERIFY) 
    : (gen->variant & (~HRI_VARFLAG_VERIFY));
  SAMPLE = _unur_hri_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hri_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_HRI ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HRI_PAR,NULL);
  gen = _unur_hri_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_hri_check_par(gen) != UNUR_SUCCESS) {
    _unur_hri_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_hri_debug_init(gen);
#endif
  return gen;
} 
int
_unur_hri_reinit( struct unur_gen *gen )
{
  int rcode;
  SAMPLE = _unur_hri_getSAMPLE(gen);
  if ( (rcode = _unur_hri_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hri_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HRI_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_hri_gen) );
  COOKIE_SET(gen,CK_HRI_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_hri_getSAMPLE(gen);
  gen->destroy = _unur_hri_free;
  gen->clone = _unur_hri_clone;
  gen->reinit = _unur_hri_reinit;
  GEN->p0 = PAR->p0;                  
  GEN->left_border = 0.;             
  GEN->hrp0 = 0.;                    
  GEN->left_border = 0.;             
  return gen;
} 
int
_unur_hri_check_par( struct unur_gen *gen )
{
  if (DISTR.domain[0] < 0.)       DISTR.domain[0] = 0.;
  if (DISTR.domain[1] < INFINITY) DISTR.domain[1] = INFINITY;
  GEN->left_border = DISTR.domain[0];
  if (gen->set & HRI_SET_P0) {
    if (GEN->p0 <= GEN->left_border) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"p0 <= left boundary");
      GEN->p0 = GEN->left_border + 1.;
    }
  }
  else {
    GEN->p0 = GEN->left_border + 1.;
  }
  GEN->hrp0 = HR(GEN->p0);
  if (GEN->hrp0 <= 0. || _unur_FP_is_infinity(GEN->hrp0)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"design point p0 not valid");
    return UNUR_ERR_GEN_CONDITION;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_hri_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_hri_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HRI_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_hri_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_HRI ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HRI_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_hri_sample( struct unur_gen *gen )
{ 
  double U, V, E, X, hrx1;
  double lambda0, p1, lambda1;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HRI_GEN,INFINITY);
  lambda0 = GEN->hrp0;
  X = GEN->left_border;
  for(;;) {
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );
    E = -log(U) / lambda0;
    X += E;
    hrx1 = HR(X);
    V =  lambda0 * _unur_call_urng(gen->urng);
    if( V <= hrx1 )
      break;
  }
  if (X <= GEN->p0)
    return X;
  lambda1 = hrx1 - lambda0;
  if (lambda1 <= 0.)
    return X;
  p1 = X;
  X = GEN->p0;
  for(;;) {
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );
    E = -log(U) / lambda1;
    X += E;
    V = lambda0 + lambda1 * _unur_call_urng(gen->urng);
    if (V <= GEN->hrp0)
      break;
    if (V <= HR(X))
      break;
  }
  return ((X <= p1) ? X : p1);
} 
double
_unur_hri_sample_check( struct unur_gen *gen )
{ 
  double U, V, E, X, hrx, hrx1;
  double lambda0, p1, lambda1;
  int i0 = 0;
  int i1 = 0;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HRI_GEN,INFINITY);
  lambda0 = GEN->hrp0;
  X = GEN->left_border;
  for(i0=1;;i0++) {
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );
    E = -log(U) / lambda0;
    X += E;
    hrx1 = HR(X);
    V =  lambda0 * _unur_call_urng(gen->urng);
    if ( (X <= GEN->p0 && (1.+UNUR_EPSILON) * lambda0 < hrx1 ) || 
	 (X >= GEN->p0 && (1.-UNUR_EPSILON) * lambda0 > hrx1 ) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not increasing");
    if( V <= hrx1 )
      break;
  }
  if (X <= GEN->p0) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & HRI_DEBUG_SAMPLE)
      _unur_hri_debug_sample( gen, X, X, i0, 0 );
#endif
    return X;
  }
  lambda1 = hrx1 - lambda0;
  if (lambda1 <= 0.) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & HRI_DEBUG_SAMPLE)
      _unur_hri_debug_sample( gen, X, X, i0, 0 );
#endif
    return X;
  }
  p1 = X;
  X = GEN->p0;
  for(i1=1;;i1++) {
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );
    E = -log(U) / lambda1;
    X += E;
    V = lambda0 + lambda1 * _unur_call_urng(gen->urng);
    hrx = HR(X);
    if ( (X <= p1 && (1.+UNUR_EPSILON) * (lambda0+lambda1) < hrx ) || 
	 (X >= p1 && (1.-UNUR_EPSILON) * (lambda0+lambda1) > hrx ) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not increasing");
    if (V <= GEN->hrp0)
      break;
    if (V <= hrx)
      break;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & HRI_DEBUG_SAMPLE)
    _unur_hri_debug_sample( gen, X, p1, i0, i1 );
#endif
  return ((X <= p1) ? X : p1);
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_hri_debug_init( const struct unur_gen *gen )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HRI_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = HRI (Hazard Rate Increasing)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(log,"%s: sampling routine = _unur_hri_sample",gen->genid);
  if (gen->variant & HRI_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: design point p0 = %g  (HR(p0)=%g)",gen->genid,GEN->p0,GEN->hrp0);
  _unur_print_if_default(gen,HRI_SET_P0);
  fprintf(log,"\n%s: left boundary = %g\n",gen->genid,GEN->left_border);
  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} 
void
_unur_hri_debug_sample( const struct unur_gen *gen, 
			double x, double p1, int i0, int i1 )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HRI_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: X = %g\t(p1=%g)\t#iterations = %d + %d = %d",gen->genid,
	  x, p1, i0, i1, i0+i1);
  if (i1) 
    fprintf(log,"   2nd loop\n");
  else
    fprintf(log,"\n");
} 
#endif   
