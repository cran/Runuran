/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "arou.h"
#include "arou_struct.h"
#define AROU_VARFLAG_VERIFY     0x01u   
#define AROU_VARFLAG_USECENTER  0x02u   
#define AROU_VARFLAG_PEDANTIC   0x04u   
#define AROU_VARFLAG_USEDARS    0x10u   
#define AROU_DEBUG_SEGMENTS     0x00000010u   
#define AROU_DEBUG_SPLIT        0x00010000u   
#define AROU_DEBUG_DARS         0x00020000u
#define AROU_SET_CENTER         0x001u
#define AROU_SET_STP            0x002u
#define AROU_SET_N_STP          0x004u
#define AROU_SET_GUIDEFACTOR    0x010u
#define AROU_SET_MAX_SQHRATIO   0x020u
#define AROU_SET_MAX_SEGS       0x040u
#define AROU_SET_USE_DARS       0x100u
#define AROU_SET_DARS_FACTOR    0x200u
#define GENTYPE "AROU"         
static struct unur_gen *_unur_arou_init( struct unur_par *par );
static struct unur_gen *_unur_arou_create( struct unur_par *par );
static struct unur_gen *_unur_arou_clone( const struct unur_gen *gen );
static void _unur_arou_free( struct unur_gen *gen);
static double _unur_arou_sample( struct unur_gen *gen );
static double _unur_arou_sample_check( struct unur_gen *gen );
static int _unur_arou_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen );
static int _unur_arou_get_starting_segments( struct unur_gen *gen );
static double _unur_arou_compute_x( double v, double u );
static int _unur_arou_run_dars( struct unur_gen *gen );
static struct unur_arou_segment *_unur_arou_segment_new( struct unur_gen *gen, double x, double fx );
static int _unur_arou_segment_parameter( struct unur_gen *gen, struct unur_arou_segment *seg );
static int _unur_arou_segment_split( struct unur_gen *gen, struct unur_arou_segment *seg_old, double x, double fx );
static int _unur_arou_make_guide_table( struct unur_gen *gen );
static double _unur_arou_segment_arcmean( struct unur_arou_segment *seg );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_arou_debug_init( const struct unur_par *par, const struct unur_gen *gen );
static void _unur_arou_debug_dars_start( const struct unur_gen *gen );
static void _unur_arou_debug_dars( const struct unur_gen *gen );
static void _unur_arou_debug_free( const struct unur_gen *gen );
static void _unur_arou_debug_segments( const struct unur_gen *gen );
static void _unur_arou_debug_split_start( const struct unur_gen *gen, 
					  const struct unur_arou_segment *seg,
					  double x, double fx );
static void _unur_arou_debug_split_stop( const struct unur_gen *gen, 
					 const struct unur_arou_segment *seg_left,
					 const struct unur_arou_segment *seg_right );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_arou_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_arou_par*)par->datap) 
#define GEN       ((struct unur_arou_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))  
#define dPDF(x)   _unur_cont_dPDF((x),(gen->distr)) 
#define _unur_arou_getSAMPLE(gen) \
   ( ((gen)->variant & AROU_VARFLAG_VERIFY) \
     ? _unur_arou_sample_check : _unur_arou_sample )
struct unur_par *
unur_arou_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL(GENTYPE,distr,NULL);
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL;
  }
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of PDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_arou_par) );
  COOKIE_SET(par,CK_AROU_PAR);
  par->distr              = distr;  
  PAR->guide_factor        = 2.;     
  PAR->darsfactor          = 0.99;   
  PAR->starting_cpoints    = NULL;   
  PAR->n_starting_cpoints  = 30;     
  PAR->max_segs            = 100;    
  PAR->max_ratio           = 0.99;   
  par->method   = UNUR_METH_AROU;             
  par->variant  = ( AROU_VARFLAG_USECENTER |
		    AROU_VARFLAG_USEDARS );   
  par->set      = 0u;                          
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = par->urng;               
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_arou_init;
  return par;
} 
int
unur_arou_set_usedars( struct unur_par *par, int usedars )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  par->variant = (usedars) ? (par->variant | AROU_VARFLAG_USEDARS) : (par->variant & (~AROU_VARFLAG_USEDARS));
  par->set |= AROU_SET_USE_DARS;
  return UNUR_SUCCESS;
} 
int
unur_arou_set_darsfactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->darsfactor = factor;
  par->set |= AROU_SET_DARS_FACTOR;
  return UNUR_SUCCESS;
} 
int
unur_arou_set_cpoints( struct unur_par *par, int n_stp, const double *stp )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }
  if (stp) 
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,
		      "starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
  PAR->starting_cpoints = stp;
  PAR->n_starting_cpoints = n_stp;
  par->set |= AROU_SET_N_STP | ((stp) ? AROU_SET_STP : 0);
  return UNUR_SUCCESS;
} 
int
unur_arou_set_guidefactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->guide_factor = factor;
  par->set |= AROU_SET_GUIDEFACTOR;
  return UNUR_SUCCESS;
} 
int
unur_arou_set_max_sqhratio( struct unur_par *par, double max_ratio )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ratio = max_ratio;
  par->set |= AROU_SET_MAX_SQHRATIO;
  return UNUR_SUCCESS;
} 
double
unur_arou_get_sqhratio( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, AROU, UNUR_INFINITY );
  return (GEN->Asqueeze / GEN->Atotal);
} 
double
unur_arou_get_hatarea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, AROU, UNUR_INFINITY );
  return GEN->Atotal;
} 
double
unur_arou_get_squeezearea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, AROU, UNUR_INFINITY );
  return GEN->Asqueeze;
} 
int
unur_arou_set_max_segments( struct unur_par *par, int max_segs )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  if (max_segs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of segments < 1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_segs = max_segs;
  par->set |= AROU_SET_MAX_SEGS;
  return UNUR_SUCCESS;
} 
int
unur_arou_set_usecenter( struct unur_par *par, int usecenter )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  par->variant = (usecenter) ? (par->variant | AROU_VARFLAG_USECENTER) : (par->variant & (~AROU_VARFLAG_USECENTER));
  return UNUR_SUCCESS;
} 
int
unur_arou_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  par->variant = (verify) ? (par->variant | AROU_VARFLAG_VERIFY) : (par->variant & (~AROU_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_arou_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, AROU, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= AROU_VARFLAG_VERIFY;
  else
    gen->variant &= ~AROU_VARFLAG_VERIFY;
  SAMPLE = _unur_arou_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
int
unur_arou_set_pedantic( struct unur_par *par, int pedantic )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AROU );
  par->variant = (pedantic) ? (par->variant | AROU_VARFLAG_PEDANTIC) : (par->variant & (~AROU_VARFLAG_PEDANTIC));
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_arou_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  int i,k;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_AROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_AROU_PAR,NULL);
  gen = _unur_arou_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }
  if (_unur_arou_get_starting_cpoints(par,gen)!=UNUR_SUCCESS ) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_arou_debug_init(par,gen);
#endif
    _unur_par_free(par); _unur_arou_free(gen);
    return NULL;
  }
  if ( _unur_arou_get_starting_segments(gen)!=UNUR_SUCCESS ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_arou_debug_init(par,gen);
#endif
    _unur_par_free(par); _unur_arou_free(gen);
    return NULL;
  }
  if (GEN->n_segs > GEN->max_segs) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"maximal number of segments too small. increase.");
    GEN->max_segs = GEN->n_segs;
  }
  if (gen->variant & AROU_VARFLAG_USEDARS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & AROU_DEBUG_DARS) {
      _unur_arou_make_guide_table(gen);
      _unur_arou_debug_init(par,gen);
      _unur_arou_debug_dars_start(gen);
    }
#endif
    for (i=0; i<3; i++) {
      if ( _unur_arou_run_dars(gen)!=UNUR_SUCCESS ) {
	_unur_par_free(par); _unur_arou_free(gen);
	return NULL;
      }
      _unur_arou_make_guide_table(gen);
      if (GEN->n_segs < GEN->max_segs) {
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) {
      if (gen->debug & AROU_DEBUG_DARS)
	_unur_arou_debug_dars(gen);
      else 
  	_unur_arou_debug_init(par,gen);
    }
#endif
  }
  else { 
    _unur_arou_make_guide_table(gen);
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_arou_debug_init(par,gen);
#endif
  }
  _unur_par_free(par);
  if (GEN->Atotal <= 0. || !_unur_isfinite(GEN->Atotal)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points");
    _unur_arou_free(gen);
    return NULL;
  }
  return gen;
} 
struct unur_gen *
_unur_arou_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_AROU_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_arou_gen) );
  COOKIE_SET(gen,CK_AROU_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_arou_getSAMPLE(gen);
  gen->destroy = _unur_arou_free;
  gen->clone = _unur_arou_clone;
  GEN->seg         = NULL;
  GEN->n_segs      = 0;
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;
  GEN->guide_factor = PAR->guide_factor; 
  GEN->max_segs = PAR->max_segs;      
#ifdef UNUR_ENABLE_INFO
  GEN->max_segs_info = PAR->max_segs;   
#endif
  GEN->max_ratio = PAR->max_ratio;    
  GEN->darsfactor = PAR->darsfactor;
  if ( (gen->distr->set & UNUR_DISTR_SET_CENTER) ||
       (gen->distr->set & UNUR_DISTR_SET_MODE) ) {
    GEN->center = unur_distr_cont_get_center(gen->distr);
    GEN->center = _unur_max(GEN->center,DISTR.BD_LEFT);
    GEN->center = _unur_min(GEN->center,DISTR.BD_RIGHT);
    gen->set |= AROU_SET_CENTER;
  }
  else {
    GEN->center = 0.;
    gen->variant &= ~AROU_VARFLAG_USECENTER;
  }
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_arou_info;
#endif
  return(gen);
} 
struct unur_gen *
_unur_arou_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_arou_gen*)clone->datap)
  struct unur_gen *clone;
  struct unur_arou_segment *seg,*next, *clone_seg, *clone_prev;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  clone_seg = NULL;
  clone_prev = NULL;
  for (seg = GEN->seg; seg != NULL; seg = next) {
    clone_seg = _unur_xmalloc( sizeof(struct unur_arou_segment) );
    memcpy( clone_seg, seg, sizeof(struct unur_arou_segment) );
    if (clone_prev == NULL) {
      CLONE->seg = clone_seg;
    }
    else {
      clone_prev->next = clone_seg;
      clone_prev->rtp  = clone_seg->ltp;
      clone_prev->drtp = clone_seg->dltp;
    }
    next = seg->next;
    clone_prev = clone_seg;
  }
  if (clone_seg) clone_seg->next = NULL;
  CLONE->guide = NULL;
  _unur_arou_make_guide_table(clone);
  return clone;
#undef CLONE
} 
void
_unur_arou_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_AROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  SAMPLE = NULL;   
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_arou_debug_free(gen);
#endif
  {
    struct unur_arou_segment *seg,*next;
    for (seg = GEN->seg; seg != NULL; seg = next) {
      next = seg->next;
      free(seg);
    }
  }
  if (GEN->guide) free(GEN->guide);
  _unur_generic_free(gen);
} 
double
_unur_arou_sample( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_arou_segment *seg;
  int result_split;
  double R,R1,R2,R3,tmp,x,fx,u;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_INFINITY);
  urng = gen->urng;
  while (1) {
    R = _unur_call_urng(urng);
    seg =  GEN->guide[(int) (R * GEN->guide_size)];
    R *= GEN->Atotal;
    while (seg->Acum < R) {
      seg = seg->next;
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_INFINITY);
    R = seg->Acum - R;
    if (R < seg->Ain) {
      return( ( seg->Ain * seg->rtp[0] + R * (seg->ltp[0] - seg->rtp[0]) ) /
	      ( seg->Ain * seg->rtp[1] + R * (seg->ltp[1] - seg->rtp[1]) ) );
    }
    else {
      urng = gen->urng_aux;
      R1 = (R - seg->Ain) / seg->Aout;  
      R2 = _unur_call_urng(urng);
      if (R1>R2) { tmp = R1; R1=R2; R2=tmp; }  
      R3 = 1.-R2;
      R2 -= R1;
      u = seg->ltp[1]*R1 + seg->rtp[1]*R2 + seg->mid[1]*R3;
      x = (seg->ltp[0]*R1 + seg->rtp[0]*R2 + seg->mid[0]*R3) / u;
      fx = PDF(x);
      if (GEN->n_segs < GEN->max_segs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  result_split = _unur_arou_segment_split(gen,seg,x,fx);
	  if ( !(result_split == UNUR_SUCCESS || result_split == UNUR_ERR_SILENT) ) {
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	    if (gen->variant & AROU_VARFLAG_PEDANTIC) {
	      SAMPLE = _unur_sample_cont_error;
	      return UNUR_INFINITY;
	    }
	  }
	  else {
	    _unur_arou_make_guide_table(gen);
	  }
	}
	else 
	  GEN->max_segs = GEN->n_segs;
      }
      if (u*u <= fx) 
	return x;
    }
  }
} 
double
_unur_arou_sample_check( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_arou_segment *seg;
  int result_split;
  double R,R1,R2,R3,tmp,x,fx,u,sqx,a;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_INFINITY);
  urng = gen->urng;
  while (1) {
    R = _unur_call_urng(urng);
    seg =  GEN->guide[(int) (R * GEN->guide_size)];
    R *= GEN->Atotal;
    while (seg->Acum < R) {
      seg = seg->next;
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_INFINITY);
    R = seg->Acum - R;
    if (R < seg->Ain) {
      x = ( ( seg->Ain * seg->rtp[0] + R * (seg->ltp[0] - seg->rtp[0]) ) /
	    ( seg->Ain * seg->rtp[1] + R * (seg->ltp[1] - seg->rtp[1]) ) );
      fx = PDF(x);
      a = ( (seg->rtp[0] - x * seg->rtp[1]) / 
	    (seg->rtp[0] - seg->ltp[0] + x * (seg->ltp[1] - seg->rtp[1])) );
      sqx = a * seg->ltp[1] + (1.-a) * seg->rtp[1];
      if (sqx*sqx > fx * (1.+UNUR_EPSILON))
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave.");
      return x;
    }
    else {
      urng = gen->urng_aux;
      R1 = (R - seg->Ain) / seg->Aout;  
      R2 = _unur_call_urng(urng);
      if (R1>R2) { tmp = R1; R1=R2; R2=tmp; }  
      R3 = 1.-R2;
      R2 -= R1;
      u = seg->ltp[1]*R1 + seg->rtp[1]*R2 + seg->mid[1]*R3;
      x = (seg->ltp[0]*R1 + seg->rtp[0]*R2 + seg->mid[0]*R3) / u;
      fx = PDF(x);
      a = ( (seg->rtp[0] - x * seg->rtp[1]) / 
	    (seg->rtp[0] - seg->ltp[0] + x * (seg->ltp[1] - seg->rtp[1])) );
      sqx = a * seg->ltp[1] + (1.-a) * seg->rtp[1];
      if (sqx*sqx > fx)
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave.");
      if (GEN->n_segs < GEN->max_segs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  result_split = _unur_arou_segment_split(gen,seg,x,fx);
	  if ( !(result_split == UNUR_SUCCESS || result_split == UNUR_ERR_SILENT) ) {
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	    if (gen->variant & AROU_VARFLAG_PEDANTIC) {
	      SAMPLE = _unur_sample_cont_error;
	      return UNUR_INFINITY;
	    }
	  }
	  else {
	    _unur_arou_make_guide_table(gen);
	  }
	}
	else 
	  GEN->max_segs = GEN->n_segs;
      }
      if (u*u <= fx) 
	return x;
    }
  }
} 
int
_unur_arou_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
{
  struct unur_arou_segment *seg, *seg_new;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, fx, fx_last;
  int i, use_center, is_increasing;  
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_AROU_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  use_center = (gen->variant & AROU_VARFLAG_USECENTER) ? TRUE : FALSE;
  GEN->n_segs = 0;
  if (!PAR->starting_cpoints) {
    left_angle =  ( DISTR.BD_LEFT  <= -UNUR_INFINITY ) ? -M_PI/2. : atan(DISTR.BD_LEFT  - GEN->center);  
    right_angle = ( DISTR.BD_RIGHT >= UNUR_INFINITY )  ? M_PI/2.  : atan(DISTR.BD_RIGHT - GEN->center);
    diff_angle = (right_angle-left_angle) / (PAR->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   
  x = x_last = DISTR.BD_LEFT;
  fx = fx_last = (x <= -UNUR_INFINITY) ? 0. : PDF(x);
  seg = GEN->seg = _unur_arou_segment_new( gen, x, fx );
  if (seg == NULL) return UNUR_ERR_GEN_CONDITION;  
  is_increasing = 1; 
  for( i=0; i<=PAR->n_starting_cpoints; i++ ) {
    if (i < PAR->n_starting_cpoints) {
      if (PAR->starting_cpoints) {   
	x = PAR->starting_cpoints[i];
	if (x <= DISTR.BD_LEFT || x >= DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
	if (x<=x_last) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting points not increasing -> skip");
	  continue;
	}
      }
      else {
	angle += diff_angle;
	x = tan( angle ) + GEN->center;
      }
    }
    else {
      x = DISTR.BD_RIGHT;
    }
    if (use_center && x >= GEN->center) {
      use_center = FALSE;   
      if (x>GEN->center) {
	x = GEN->center;   
	--i;              
	if (!PAR->starting_cpoints)
	  angle -= diff_angle; 
      }
    }
    fx = (x >= UNUR_INFINITY) ? 0. : PDF(x);
    if (!is_increasing && fx > fx_last) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (fx <= 0. && fx_last <= 0.) {
      if (is_increasing) {
	if (i<PAR->n_starting_cpoints) {
	  seg->dltp[0] = -1.;
	  seg->dltp[1] = x;
	  x_last = x;
	  continue;   
	}
      }
      else
	break;
    }
    seg_new = _unur_arou_segment_new( gen, x, fx );
    if (seg_new == NULL) {
      seg->next = NULL;  
      return UNUR_ERR_GEN_CONDITION;
    }
    seg->next =seg_new;
    seg->rtp = seg_new->ltp;
    seg->drtp = seg_new->dltp;
    seg = seg_new;
    if (is_increasing && fx < fx_last)
      is_increasing = 0;
    x_last = x;
    fx_last = fx;
  }
  seg->Ain = seg->Aout = 0.;
  seg->Acum = UNUR_INFINITY;
  seg->next = NULL;         
  --(GEN->n_segs);           
  return UNUR_SUCCESS;
} 
int
_unur_arou_get_starting_segments( struct unur_gen *gen )
{
#define MAX_IT   (1000)      
  struct unur_arou_segment *seg, *seg_new, *seg_tmp; 
  double x,fx;              
  int n_it = 0;             
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  for( seg=GEN->seg; seg->next != NULL; ) {
    switch (_unur_arou_segment_parameter(gen, seg)) {
    case UNUR_SUCCESS:     
      seg = seg->next;
      continue;
    case UNUR_ERR_SILENT:    
      if (seg->next != NULL) {
	seg_tmp = seg->next;
	seg->next = seg->next->next;
	seg->rtp = seg->next->ltp;
	seg->drtp = seg->next->dltp;
	free(seg_tmp);
	--(GEN->n_segs);
      }
      else { 
	seg->Ain = seg->Aout = 0.;
	seg->Acum = UNUR_INFINITY;
      }
      continue;
    case UNUR_ERR_INF:    
      break;
    default:     
      return UNUR_ERR_GEN_CONDITION;
    }
    ++n_it;
    if (n_it > MAX_IT) {
      return UNUR_ERR_GEN_CONDITION;
    }
    x = _unur_arou_segment_arcmean(seg);  
    fx = PDF(x);
    if (fx < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF < 0");
      return UNUR_ERR_GEN_DATA;
    }
    if (GEN->n_segs >= GEN->max_segs) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded envelope!");
      return UNUR_ERR_GEN_CONDITION;
    }
    seg_new = _unur_arou_segment_new( gen, x, fx );
    if (seg_new == NULL) return UNUR_ERR_GEN_CONDITION;  
    seg_new->next = seg->next;
    seg->next = seg_new;
    seg_new->rtp = seg->rtp;
    seg_new->drtp = seg->drtp;
    seg->rtp = seg_new->ltp;
    seg->drtp = seg_new->dltp;
  }
  return UNUR_SUCCESS;
} 
struct unur_arou_segment *
_unur_arou_segment_new( struct unur_gen *gen, double x, double fx )
{
  struct unur_arou_segment *seg;
  double u,v,dfx;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,NULL);
  if (fx<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.");
    return NULL;
  }
  if (_unur_FP_is_infinity(fx)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return NULL;
  }
  seg = _unur_xmalloc( sizeof(struct unur_arou_segment) );
  seg->next = NULL; 
  ++(GEN->n_segs);   
  COOKIE_SET(seg,CK_AROU_SEG);
  seg->Ain = seg->Aout = seg->Acum = 0.;
  seg->mid[0] = seg->mid[1] = 0.;
  if ( _unur_iszero(fx) ) {
    seg->ltp[0] = 0.;   
    seg->ltp[1] = 0.;
    if (x <= -UNUR_INFINITY || x >= UNUR_INFINITY ) {
      seg->dltp[0] = 0.;   
      seg->dltp[1] = 1.;   
      seg->dltp[2] = 0.;   
    }
    else {
      seg->dltp[0] = -1.;  
      seg->dltp[1] = x;    
      seg->dltp[2] = 0.;   
    }
    return seg;
  }
  u = sqrt( fx );
  v = x * u;
  seg->ltp[0] = v;
  seg->ltp[1] = u; 
  dfx = dPDF(x);
  if ( dfx > -UNUR_INFINITY && dfx < UNUR_INFINITY ) {
    seg->dltp[0] = -dfx / u;                 
    seg->dltp[1] = 2 * u + dfx * x / u;      
    seg->dltp[2] = seg->dltp[0] * v + seg->dltp[1] * u;
    return seg;
  }
  seg->dltp[0] = -u;   
  seg->dltp[1] = v;    
  seg->dltp[2] = 0.;
  return seg;
} 
#define MAX_NORM_OF_INTERSECTION_POINT  1.e6
int
_unur_arou_segment_parameter( struct unur_gen *gen, struct unur_arou_segment *seg )
{
  double coeff_det, cramer_det[2];
  double norm_vertex;      
  double diff_tangents;    
  double det_bound;        
  double tmp_a, tmp_b;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(seg,UNUR_ERR_NULL);  COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
  norm_vertex = fabs(seg->ltp[0]) + fabs(seg->ltp[1]) + fabs(seg->rtp[0]) + fabs(seg->rtp[1]);
  seg->Ain = (seg->ltp[1] * seg->rtp[0] - seg->ltp[0] * seg->rtp[1]) / 2.;
  if( seg->Ain < 0. ) {
    if (fabs(seg->Ain) < 1.e-8 * norm_vertex) {
      seg->Ain = seg->Aout = 0.;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }
    return UNUR_ERR_SILENT;
  }
  coeff_det     = seg->dltp[0] * seg->drtp[1] - seg->dltp[1] * seg->drtp[0];
  cramer_det[0] = seg->dltp[2] * seg->drtp[1] - seg->dltp[1] * seg->drtp[2];
  cramer_det[1] = seg->dltp[0] * seg->drtp[2] - seg->dltp[2] * seg->drtp[0];
  det_bound = fabs(coeff_det) * norm_vertex * MAX_NORM_OF_INTERSECTION_POINT;
  diff_tangents = ( fabs(seg->dltp[0] - seg->drtp[0]) + fabs(seg->dltp[1] - seg->drtp[1])
		    + fabs(seg->dltp[2] - seg->drtp[2]) );
  if (!_unur_iszero(coeff_det) && !_unur_iszero(diff_tangents)) {
    if ( fabs(cramer_det[0]) > det_bound || fabs(cramer_det[1]) > det_bound ) {
      seg->Aout = UNUR_INFINITY;
      return UNUR_ERR_INF;
    }
    seg->mid[0] = cramer_det[0] / coeff_det;
    seg->mid[1] = cramer_det[1] / coeff_det;
    seg->Aout = ( (seg->ltp[0] - seg->mid[0]) * (seg->rtp[1] - seg->mid[1])
		  - (seg->ltp[1] - seg->mid[1]) * (seg->rtp[0] - seg->mid[0])) / 2.;
    if( seg->mid[1] < 0. ) {
      seg->Aout = UNUR_INFINITY;
      return UNUR_ERR_INF;
    }
    if ( seg->Aout > 0. ) {
      tmp_a = seg->mid[0] * seg->ltp[1];
      tmp_b = seg->ltp[0] * seg->mid[1];
      if ( ! _unur_FP_less(tmp_a, tmp_b) ) {
	tmp_a = seg->mid[0] * seg->rtp[1];
	tmp_b = seg->rtp[0] * seg->mid[1];
	if ( ! _unur_FP_greater(tmp_a, tmp_b) )
	  return UNUR_SUCCESS;
      }
    }
    if ( !_unur_iszero(seg->ltp[1]) && !_unur_iszero(seg->rtp[1]) ) {
      tmp_a = seg->ltp[0] * seg->rtp[1];
      tmp_b = seg->rtp[0] * seg->ltp[1];
      if ( _unur_FP_equal(tmp_a, tmp_b) ) {
	seg->Ain = seg->Aout = 0.;
	return UNUR_ERR_SILENT;
      }
    }
    if (!(fabs(seg->Aout) < fabs(seg->Ain) * UNUR_EPSILON)) {
      seg->Aout = UNUR_INFINITY;
      return UNUR_ERR_INF;
    }
  }
  seg->mid[0] =  0.5 * (seg->ltp[0] + seg->rtp[0]);
  seg->mid[1] =  0.5 * (seg->ltp[1] + seg->rtp[1]);
  seg->Aout = 0.;
  return UNUR_SUCCESS;
} 
#undef MAX_NORM_INTERSECTION
int
_unur_arou_segment_split( struct unur_gen *gen, struct unur_arou_segment *seg_oldl, double x, double fx )
{
  struct unur_arou_segment *seg_newr;    
  struct unur_arou_segment seg_bak;      
  double Adiff;
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(seg_oldl,UNUR_ERR_NULL); COOKIE_CHECK(seg_oldl,CK_AROU_SEG,UNUR_ERR_COOKIE);
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & AROU_DEBUG_SPLIT) 
      _unur_arou_debug_split_start( gen,seg_oldl,x,fx );
#endif
  if (GEN->n_segs * seg_oldl->Aout / (GEN->Atotal - GEN->Asqueeze) < GEN->darsfactor )
    return UNUR_SUCCESS;
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }
  memcpy(&seg_bak, seg_oldl, sizeof(struct unur_arou_segment));
  if (fx <= 0.) {
    if (seg_oldl->rtp[1] <= 0. && seg_oldl->rtp[0] <= 0. ) {
      seg_oldl->drtp[1] = x;    
    }
    else if (seg_oldl->ltp[1] <= 0. && seg_oldl->ltp[0] <= 0. ) {
      seg_oldl->dltp[1] = x;    
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    if( _unur_arou_segment_parameter(gen,seg_oldl)!=UNUR_SUCCESS ) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot chop segment at given point");
      memcpy(seg_oldl, &seg_bak, sizeof(struct unur_arou_segment));
      return UNUR_ERR_SILENT;
    }
    seg_newr = seg_oldl;
  }
  else {  
    seg_newr = _unur_arou_segment_new(gen,x,fx);
    if (seg_newr == NULL) return UNUR_ERR_GEN_DATA;  
    seg_newr->next = seg_oldl->next;
    seg_oldl->next = seg_newr;
    seg_newr->rtp = seg_oldl->rtp;
    seg_newr->drtp = seg_oldl->drtp;
    seg_oldl->rtp = seg_newr->ltp;
    seg_oldl->drtp = seg_newr->dltp;
    if( _unur_arou_segment_parameter(gen,seg_oldl)!=UNUR_SUCCESS || 
	_unur_arou_segment_parameter(gen,seg_newr)!=UNUR_SUCCESS 
	 ) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split segment at given point.");
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug & AROU_DEBUG_SPLIT)
	_unur_arou_debug_split_stop( gen,seg_oldl,seg_newr );
#endif
      memcpy(seg_oldl, &seg_bak, sizeof(struct unur_arou_segment));
      if (seg_newr) {
	--(GEN->n_segs); 
	free( seg_newr );
      }
      return UNUR_ERR_SILENT;
    }
  }
  Adiff =  - seg_bak.Ain  + seg_oldl->Ain  + ((seg_newr!=seg_oldl) ? seg_newr->Ain : 0. );
  GEN->Asqueeze += Adiff;
  Adiff += - seg_bak.Aout + seg_oldl->Aout + ((seg_newr!=seg_oldl) ? seg_newr->Aout : 0. );
  GEN->Atotal += Adiff;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & AROU_DEBUG_SPLIT) 
    _unur_arou_debug_split_stop( gen,seg_oldl,seg_newr );
#endif
  return UNUR_SUCCESS;
} 
double
_unur_arou_compute_x( double v, double u )
{
  if (!_unur_iszero(u)) return v/u;
  else if (v<0.)        return -UNUR_INFINITY;
  else                  return UNUR_INFINITY;
} 
int
_unur_arou_run_dars( struct unur_gen *gen )
{
  struct unur_arou_segment *seg, *seg_next;
  double Atot, Asqueezetot;    
  double Alimit;               
  int n_splitted = 1;          
  int splitted;                
  double xl, xr;               
  double xsp, fxsp;            
  CHECK_NULL(gen,UNUR_ERR_NULL);     COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;
  Atot = 0.;            
  Asqueezetot = 0.;     
  for (seg = GEN->seg; seg != NULL; seg = seg->next ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
    Asqueezetot += seg->Ain;
    Atot += seg->Ain + seg->Aout;
  }
  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_segs < GEN->max_segs) ) {
    if (GEN->n_segs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_segs );
    else
      Alimit = 0.; 
    n_splitted = 0;
    for (seg = GEN->seg; seg->next != NULL; seg = seg->next ) {
      COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
      if (GEN->n_segs >= GEN->max_segs)
	break;
      if (seg->Aout <= Alimit) 
	continue;  
      seg_next = seg->next;
      xl = _unur_arou_compute_x(seg->ltp[0],seg->ltp[1]);
      xr = _unur_arou_compute_x(seg->rtp[0],seg->rtp[1]);
      if (xl>xr) xl = -UNUR_INFINITY;   
      if ( _unur_FP_is_minus_infinity(xl)
	   && _unur_FP_same(seg->dltp[0],-1.) && _unur_iszero(seg->dltp[2]) )
	xl = seg->dltp[1];
      if ( _unur_FP_is_infinity(xr)
	   && _unur_FP_same(seg->drtp[0],-1.) && _unur_iszero(seg->drtp[2]) )
	xr = seg->drtp[1];
      xsp = _unur_arcmean(xl,xr);
      fxsp = PDF(xsp);
      splitted = _unur_arou_segment_split(gen, seg, xsp, fxsp);
      if (splitted == UNUR_SUCCESS) {
      	++n_splitted;
	if (seg->next != seg_next)
	  seg = seg->next;
      }
      else if (splitted != UNUR_ERR_SILENT) {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	return UNUR_ERR_GEN_CONDITION;
      }
    }
    if (n_splitted == 0) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: no intervals could be splitted.");
      break;
    }
  }
  if ( GEN->max_ratio * GEN->Atotal > GEN->Asqueeze ) {
    if ( GEN->n_segs >= GEN->max_segs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    GEN->max_segs = GEN->n_segs;
  }
  return UNUR_SUCCESS;
} 
int
_unur_arou_make_guide_table( struct unur_gen *gen )
{
  struct unur_arou_segment *seg;
  double Acum, Aincum, Astep;
  int max_guide_size;
  int j;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  if (!GEN->guide) {
    max_guide_size = (GEN->guide_factor > 0.) ? ((int)(GEN->max_segs * GEN->guide_factor)) : 1;
    if (max_guide_size <= 0) max_guide_size = 1;   
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_arou_segment*) );
  }
  Acum = 0.;       
  Aincum = 0.;     
  for (seg = GEN->seg; seg != NULL; seg = seg->next ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
    Acum += seg->Ain + seg->Aout;
    Aincum += seg->Ain;
    seg->Acum = Acum;
  }
  GEN->Atotal = Acum;
  GEN->Asqueeze = Aincum;
  GEN->guide_size = (int)(GEN->n_segs * GEN->guide_factor);
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, seg=GEN->seg; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
    while( seg->Acum < Acum )
      if( seg->next != NULL )    
        seg = seg->next;
      else {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = seg;
    Acum += Astep;
  }
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = seg;
  return UNUR_SUCCESS;
} 
double
_unur_arou_segment_arcmean( struct unur_arou_segment *seg )
{
  double xl, xr;
  CHECK_NULL(seg,UNUR_INFINITY);  COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_INFINITY);
  xl = (seg->ltp[1] > 0.) ? (seg->ltp[0] / seg->ltp[1]) :
    ( _unur_iszero(seg->dltp[0]) ? -UNUR_INFINITY : (seg->dltp[1]) );
  xr = (seg->rtp[1] > 0.) ? (seg->rtp[0] / seg->rtp[1]) :
    ( _unur_iszero(seg->drtp[0]) ? UNUR_INFINITY : (seg->drtp[1]) );
  return _unur_arcmean(xl,xr);
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_arou_debug_init( const struct unur_par *par, const struct unur_gen *gen )
{
  FILE *LOG;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_AROU_PAR,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ratio-of-uniforms method with enveloping polygon\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_arou_sample",gen->genid);
  if (gen->variant & AROU_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: center = %g",gen->genid,GEN->center);
  _unur_print_if_default(gen,AROU_SET_CENTER);
  if (gen->variant & AROU_VARFLAG_USECENTER)
    fprintf(LOG,"\n%s: use center as construction point",gen->genid);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: maximum number of segments         = %d",gen->genid,GEN->max_segs);
  _unur_print_if_default(gen,AROU_SET_MAX_SEGS);
  fprintf(LOG,"\n%s: bound for ratio  Asqueeze / Atotal = %g%%",gen->genid,GEN->max_ratio*100.);
  _unur_print_if_default(gen,AROU_SET_MAX_SQHRATIO);
  fprintf(LOG,"\n%s:\n",gen->genid);
  if (gen->variant & AROU_VARFLAG_USEDARS) {
    fprintf(LOG,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(gen,AROU_SET_USE_DARS);
    fprintf(LOG,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(gen,AROU_SET_DARS_FACTOR);
  }
  else {
    fprintf(LOG,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(gen,AROU_SET_USE_DARS);
  }
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: sampling from list of segments: indexed search (guide table method)\n",gen->genid);
  fprintf(LOG,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN->guide_factor);
  _unur_print_if_default(gen,AROU_SET_GUIDEFACTOR);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fprintf(LOG,"%s: number of starting points = %d",gen->genid,PAR->n_starting_cpoints);
  _unur_print_if_default(gen,AROU_SET_N_STP);
  fprintf(LOG,"\n%s: starting points:",gen->genid);
  if (gen->set & AROU_SET_STP)
    for (i=0; i<PAR->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",PAR->starting_cpoints[i]);
    }
  else
    fprintf(LOG," use \"equidistribution\" rule [default]");
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_arou_debug_segments(gen);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void 
_unur_arou_debug_dars_start( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: DARS started **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(gen,AROU_SET_DARS_FACTOR);
  fprintf(LOG,"\n%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_arou_debug_dars( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS finished **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_arou_debug_segments(gen);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_arou_debug_free( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_arou_debug_segments(gen);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_arou_debug_segments( const struct unur_gen *gen )
{
  FILE *LOG;
  struct unur_arou_segment *seg;
  double sAin, sAout, Atotal;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:Segments: %d\n",gen->genid,GEN->n_segs);
  if ((gen->debug & AROU_DEBUG_SEGMENTS) && GEN->seg != NULL) {
    fprintf(LOG,"%s: Nr.\t    left touching point\t\t   intersection point\t\t tangent at left touching point\n",gen->genid);
    for (seg = GEN->seg, i=0; seg->next!=NULL; seg=seg->next, i++) {
      COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID); 
      fprintf(LOG,"%s:[%3d]: (%-12.6g,%-12.6g)   (%-12.6g,%-12.6g)   (%-12.6g,%-12.6g,%-12.6g)\n", gen->genid, i,
	      seg->ltp[0],seg->ltp[1],
	      seg->mid[0],seg->mid[1],
	      seg->dltp[0],seg->dltp[1],seg->dltp[2]);
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID); 
    fprintf(LOG,"%s:[...]: (%-12.6g,%-12.6g)\n", gen->genid,seg->ltp[0],seg->ltp[1]);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  if (GEN->Atotal <= 0.) {
    fprintf(LOG,"%s: Construction of enveloping polygon not successful\n",gen->genid);
    fprintf(LOG,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!!!!!!!\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    Atotal = -1.;   
  }
  else {
    Atotal = GEN->Atotal;
  }
  if ((gen->debug & AROU_DEBUG_SEGMENTS) && GEN->seg != NULL) {
    fprintf(LOG,"%s:Areas in segments:\n",gen->genid);
    fprintf(LOG,"%s: Nr.\t inside squeeze\t\t   outside squeeze\t     total segment\t\tcumulated\n",gen->genid);
    sAin = sAout = 0.;
    for (seg = GEN->seg, i=0; seg->next!=NULL; seg=seg->next, i++) {
      COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID); 
      sAin += seg->Ain;
      sAout += seg->Aout;
      fprintf(LOG,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      seg->Ain, seg->Ain * 100. / Atotal,
	      seg->Aout, seg->Aout * 100. / Atotal,
	      seg->Ain + seg->Aout, (seg->Ain + seg->Aout) * 100. / Atotal,
	      seg->Acum, seg->Acum * 100. / Atotal);
    }
    fprintf(LOG,"%s:\t----------  ---------  |  ----------  ---------  |  ----------  ---------  +\n",gen->genid);
    fprintf(LOG,"%s: Sum : %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-11.6g(%6.3f%%)\n",
	    gen->genid,
	    sAin, sAin * 100./Atotal,
	    sAout, sAout * 100./Atotal,
	    sAin + sAout, (sAin + sAout) * 100./Atotal);
    fprintf(LOG,"%s:\n",gen->genid);
  }
  fprintf(LOG,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(LOG,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(LOG,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_arou_debug_split_start( const struct unur_gen *gen,
			      const struct unur_arou_segment *seg, 
			      double x, double fx )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  CHECK_NULL(seg,RETURN_VOID);  COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: split segment at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(LOG,"%s: old segment:\n",gen->genid);
  fprintf(LOG,"%s:   left  construction point = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	  gen->genid, seg->ltp[0], seg->ltp[1], seg->ltp[0]/seg->ltp[1], sqrt(seg->ltp[1]) );
  fprintf(LOG,"%s:   intersection point       = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\n",
	  gen->genid, seg->mid[0], seg->mid[1], seg->mid[0]/seg->mid[1]);
  fprintf(LOG,"%s:   right construction point = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	  gen->genid, seg->rtp[0], seg->rtp[1], seg->rtp[0]/seg->rtp[1], sqrt(seg->rtp[1]) );
  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg->Ain, seg->Ain * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg->Aout, seg->Aout * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  (seg->Ain + seg->Aout), (seg->Ain +seg->Aout) * 100./GEN->Atotal);
  fflush(LOG);
} 
void
_unur_arou_debug_split_stop( const struct unur_gen *gen, 
			     const struct unur_arou_segment *seg_left,
			     const struct unur_arou_segment *seg_right )
{
  FILE *LOG;
  int chopped;
  CHECK_NULL(gen,RETURN_VOID);        COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  CHECK_NULL(seg_left,RETURN_VOID);   COOKIE_CHECK(seg_left,CK_AROU_SEG,RETURN_VOID);
  CHECK_NULL(seg_right,RETURN_VOID);  COOKIE_CHECK(seg_right,CK_AROU_SEG,RETURN_VOID);
  LOG = unur_get_stream();
  chopped = (seg_left==seg_right) ? 1 : 0;
  if (chopped)
    fprintf(LOG,"%s: new segment (chopped):\n",gen->genid);
  else
    fprintf(LOG,"%s: new segments:\n",gen->genid);
  fprintf(LOG,"%s:   left  construction point  = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	  gen->genid, seg_left->ltp[0], seg_left->ltp[1], seg_left->ltp[0]/seg_left->ltp[1], sqrt(seg_left->ltp[1]) );
  fprintf(LOG,"%s:   intersection point        = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\n",
	  gen->genid, seg_left->mid[0], seg_left->mid[1], seg_left->mid[0]/seg_left->mid[1] );
  if (chopped) {
    fprintf(LOG,"%s:   right construction point  = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	    gen->genid, seg_left->rtp[0], seg_left->rtp[1], seg_left->rtp[0]/seg_left->rtp[1], sqrt(seg_left->rtp[1]) );
  }
  else {
    fprintf(LOG,"%s:   middle construction point = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	    gen->genid, seg_left->rtp[0], seg_left->rtp[1], seg_left->rtp[0]/seg_left->rtp[1], sqrt(seg_left->rtp[1]) );
    fprintf(LOG,"%s:   intersection point        = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\n",
	    gen->genid, seg_right->mid[0], seg_right->mid[1], seg_right->mid[0]/seg_right->mid[1] );
    fprintf(LOG,"%s:   right construction point  = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	    gen->genid, seg_right->rtp[0], seg_right->rtp[1], seg_right->rtp[0]/seg_right->rtp[1], sqrt(seg_right->rtp[1]) );
  }
  if (!chopped) {
    fprintf(LOG,"%s: left segment:\n",gen->genid);
    fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_left->Ain, seg_left->Ain * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_left->Aout, seg_left->Aout * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    (seg_left->Ain + seg_left->Aout), (seg_left->Ain +seg_left->Aout) * 100./GEN->Atotal);
    fprintf(LOG,"%s: right segment:\n",gen->genid);
    fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_right->Ain, seg_right->Ain * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_right->Aout, seg_right->Aout * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    (seg_right->Ain + seg_right->Aout), (seg_right->Ain +seg_right->Aout) * 100./GEN->Atotal);
  }
  fprintf(LOG,"%s: total areas:\n",gen->genid);
  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (GEN->Atotal - GEN->Asqueeze) * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
#endif    
#ifdef UNUR_ENABLE_INFO
void
_unur_arou_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(distr->set & (UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_MODE )) ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
			  "You may provide a point near the mode as \"center\"."); 
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: AROU (Automatic Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g\n", GEN->Atotal);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", GEN->Atotal/(0.5*DISTR.area));
  else
    _unur_string_append(info,"<= %g\n", GEN->Atotal/GEN->Asqueeze);
  _unur_string_append(info,"   area ratio squeeze/hat = %g\n",
 		      GEN->Asqueeze/GEN->Atotal);
  _unur_string_append(info,"   # segments = %d\n", GEN->n_segs);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   max_sqhratio = %g  %s\n", GEN->max_ratio,
			(gen->set & AROU_SET_MAX_SQHRATIO) ? "" : "[default]");
    _unur_string_append(info,"   max_segments = %d  %s\n", GEN->max_segs_info,
			(gen->set & AROU_SET_MAX_SEGS) ? "" : "[default]");
    if (gen->variant & AROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    if (gen->variant & AROU_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & AROU_SET_MAX_SQHRATIO) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"max_sqhratio\" closer to 1 to decrease rejection constant." );
    if (GEN->Asqueeze/GEN->Atotal < GEN->max_ratio) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You should increase \"max_segments\" to obtain the desired rejection constant." );
    _unur_string_append(info,"\n");
  }
} 
#endif   
