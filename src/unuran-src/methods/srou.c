/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "srou.h"
#include "srou_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define SROU_VARFLAG_VERIFY   0x002u   
#define SROU_VARFLAG_SQUEEZE  0x004u   
#define SROU_VARFLAG_MIRROR   0x008u   
#define SROU_DEBUG_REINIT    0x00000010u   
#define SROU_SET_R            0x001u   
#define SROU_SET_CDFMODE      0x002u   
#define SROU_SET_PDFMODE      0x004u   
#define GENTYPE "SROU"         
static struct unur_gen *_unur_srou_init( struct unur_par *par );
static int _unur_srou_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_srou_create( struct unur_par *par );
static int _unur_srou_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_srou_clone( const struct unur_gen *gen );
static void _unur_srou_free( struct unur_gen *gen);
static double _unur_srou_sample( struct unur_gen *gen );
static double _unur_srou_sample_mirror( struct unur_gen *gen );
static double _unur_srou_sample_check( struct unur_gen *gen );
static double _unur_gsrou_sample( struct unur_gen *gen );
static double _unur_gsrou_sample_check( struct unur_gen *gen );
static int _unur_srou_rectangle( struct unur_gen *gen );
static int _unur_gsrou_envelope( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_srou_debug_init( const struct unur_gen *gen, int is_reinit );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_srou_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_srou_par*)par->datap) 
#define GEN       ((struct unur_srou_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    
static UNUR_SAMPLING_ROUTINE_CONT *
_unur_srou_getSAMPLE( struct unur_gen *gen )
{
  if (gen->variant & SROU_VARFLAG_VERIFY)
    return (gen->set & SROU_SET_R) ? _unur_gsrou_sample_check : _unur_srou_sample_check;
  else {
    if (gen->set & SROU_SET_R)
      return _unur_gsrou_sample;
    else
      return (gen->variant & SROU_VARFLAG_MIRROR) ? _unur_srou_sample_mirror : _unur_srou_sample;
  }
} 
#define SQRT2    (M_SQRT2)
struct unur_par *
unur_srou_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_srou_par) );
  COOKIE_SET(par,CK_SROU_PAR);
  par->distr    = distr;           
  PAR->r         = 1.;             
  PAR->Fmode     = -1.;            
  PAR->um        = -1.;            
  par->method   = UNUR_METH_SROU;  
  par->variant  = 0u;              
  par->set      = 0u;                  
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_srou_init;
  return par;
} 
int
unur_srou_set_r( struct unur_par *par, double r )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );
  if (r < 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r < 1");
    return UNUR_ERR_PAR_SET;
  }
  if (_unur_isone(r)) {
    PAR->r = r;
    par->set &= ~SROU_SET_R;
  }
  else {
    if (r<1.01) r = 1.01;
    PAR->r = r;
    par->set |= SROU_SET_R;
  }
  par->set &= ~SROU_SET_PDFMODE;
  return UNUR_SUCCESS;
} 
int 
unur_srou_set_cdfatmode( struct unur_par *par, double Fmode )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  PAR->Fmode = Fmode;
  par->set |= SROU_SET_CDFMODE;
  return UNUR_SUCCESS;
} 
int 
unur_srou_set_pdfatmode( UNUR_PAR *par, double fmode )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_isfinite(fmode)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return UNUR_ERR_PAR_SET;
  }
  PAR->um = (par->set & SROU_SET_R) ? pow(fmode,1./(PAR->r+1.)) : sqrt(fmode);
  par->set |= SROU_SET_PDFMODE;
  return UNUR_SUCCESS;
} 
int
unur_srou_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );
  par->variant = (verify) ? (par->variant | SROU_VARFLAG_VERIFY) : (par->variant & (~SROU_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_srou_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= SROU_VARFLAG_VERIFY;
  else
    gen->variant &= ~SROU_VARFLAG_VERIFY;
  SAMPLE = _unur_srou_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
int 
unur_srou_set_usesqueeze( struct unur_par *par, int usesqueeze )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );
  par->variant = (usesqueeze) ? (par->variant | SROU_VARFLAG_SQUEEZE) : (par->variant & (~SROU_VARFLAG_SQUEEZE));
  return UNUR_SUCCESS;
} 
int 
unur_srou_set_usemirror( struct unur_par *par, int usemirror )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, SROU );
  par->variant = (usemirror) ? (par->variant | SROU_VARFLAG_MIRROR) : (par->variant & (~SROU_VARFLAG_MIRROR));
  return UNUR_SUCCESS;
} 
int
unur_srou_chg_cdfatmode( struct unur_gen *gen, double Fmode )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  GEN->Fmode = Fmode;
  gen->set |= SROU_SET_CDFMODE;
  return UNUR_SUCCESS;
} 
int
unur_srou_chg_pdfatmode( struct unur_gen *gen, double fmode )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, SROU, UNUR_ERR_GEN_INVALID );
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_isfinite(fmode)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return UNUR_ERR_PAR_SET;
  }
  GEN->um = (gen->set & SROU_SET_R) ? pow(fmode,1./(GEN->r+1.)) : sqrt(fmode);
  gen->set |= SROU_SET_PDFMODE;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_srou_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  int rcode;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_SROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_SROU_PAR,NULL);
  if (par->set & SROU_SET_R) {
    par->variant &= ~SROU_VARFLAG_MIRROR;
    par->variant &= ~SROU_VARFLAG_SQUEEZE;
  }
  if (par->set & SROU_SET_CDFMODE)
    par->variant &= ~SROU_VARFLAG_MIRROR;
  else
    par->variant &= ~SROU_VARFLAG_SQUEEZE;
  gen = _unur_srou_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_srou_check_par(gen) != UNUR_SUCCESS) {
    _unur_srou_free(gen); return NULL;
  }
  if (gen->set & SROU_SET_R)
    rcode = _unur_gsrou_envelope( gen );
  else
    rcode = _unur_srou_rectangle( gen );
  if (rcode!=UNUR_SUCCESS) {
    _unur_srou_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_srou_debug_init(gen, FALSE);
#endif
  return gen;
} 
int
_unur_srou_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_srou_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  if (gen->set & SROU_SET_R)
    rcode = _unur_gsrou_envelope( gen );
  else
    rcode = _unur_srou_rectangle( gen );
  SAMPLE = _unur_srou_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & SROU_DEBUG_REINIT) _unur_srou_debug_init(gen,TRUE);
#endif
  return rcode;
} 
struct unur_gen *
_unur_srou_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_SROU_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_srou_gen) );
  COOKIE_SET(gen,CK_SROU_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_srou_getSAMPLE(gen);
  gen->destroy = _unur_srou_free;
  gen->clone = _unur_srou_clone;
  gen->reinit = _unur_srou_reinit;
  GEN->r     = PAR->r;                
  GEN->Fmode = PAR->Fmode;            
  GEN->um    = PAR->um;               
  GEN->vl = GEN->vr = 0.;
  GEN->xl = GEN->xr = 0.;
  GEN->p = 0.;
  GEN->a = GEN->b = 0.;
  GEN->log_ab = 0.;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_srou_info;
#endif
  return gen;
} 
int
_unur_srou_check_par( struct unur_gen *gen )
{
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)");
    if (unur_distr_cont_upd_mode(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA)) {
    if (unur_distr_cont_upd_pdfarea(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF");
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"area and/or CDF at mode");
    DISTR.mode = _unur_max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = _unur_min(DISTR.mode,DISTR.BD_RIGHT);
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_srou_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_srou_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_SROU_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_srou_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_SROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_SROU_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_srou_sample( struct unur_gen *gen )
{ 
  double U,V,X,x,xx;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,UNUR_INFINITY);
  while (1) {
    while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
    U *= GEN->um;
    V = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);
    X = V/U;
    x = X + DISTR.mode;
    if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
      continue;
    if ( (gen->variant & SROU_VARFLAG_SQUEEZE) &&
	 (X >= GEN->xl) && 
	 (X <= GEN->xr ) && 
	 (U < GEN->um) ) {
      xx = V / (GEN->um - U);
      if ( (xx >= GEN->xl) && (xx <= GEN->xr ) )
	return x;
    }
    if (U*U <= PDF(x))
      return x;
  }
} 
double
_unur_srou_sample_mirror( struct unur_gen *gen )
{ 
  double U,V,X,x,fx,fnx,uu;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,UNUR_INFINITY);
  while (1) {
    while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
    U *= GEN->um * SQRT2;
    V = 2. * (_unur_call_urng(gen->urng) - 0.5) * GEN->vr;
    X = V/U;
    x = X + DISTR.mode;
    fx  = (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) ? 0. : PDF(x);
    uu = U * U;
    if (uu <= fx)
      return x;
    x = -X + DISTR.mode;
    fnx  = (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) ? 0. : PDF(x);
    if (uu <= fx + fnx)
      return x;
  }
} 
double
_unur_srou_sample_check( struct unur_gen *gen )
{ 
  double U,uu,V,X,x,nx,fx,sfx,fnx,xfx,xfnx,xx;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,UNUR_INFINITY);
  if (gen->variant & SROU_VARFLAG_MIRROR) {
    while (1) {
      while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
      U *= GEN->um * SQRT2;
      V = 2. * (_unur_call_urng(gen->urng) - 0.5) * GEN->vr;
      X = V/U;
      x = X + DISTR.mode;
      nx = -X + DISTR.mode;
      fx  = (x  < DISTR.BD_LEFT || x  > DISTR.BD_RIGHT) ? 0. : PDF(x);
      fnx = (nx < DISTR.BD_LEFT || nx > DISTR.BD_RIGHT) ? 0. : PDF(nx);
      uu = U * U;
      xfx  = (x  - DISTR.mode) * sqrt(fx);
      xfnx = (nx - DISTR.mode) * sqrt(fnx);
      if ( ((2.+4.*DBL_EPSILON) * GEN->um*GEN->um < fx + fnx)    
	   || (xfx < (1.+UNUR_EPSILON) * GEN->vl) 
	   || (xfx > (1.+UNUR_EPSILON) * GEN->vr)
	   || (xfnx < (1.+UNUR_EPSILON) * GEN->vl) 
	   || (xfnx > (1.+UNUR_EPSILON) * GEN->vr) )
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
      if (uu <= fx)
	return x;
      if (uu <= fx + fnx)
	return nx;
    }
  }
  else { 
    while (1) {
      while ( _unur_iszero(U = _unur_call_urng(gen->urng)) );
      U *= GEN->um;
      V = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);
      X = V/U;
      x = X + DISTR.mode;
      if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	continue;
      fx = PDF(x);
      sfx = sqrt(fx);
      xfx = X * sfx;
      if ( ( sfx > (1.+DBL_EPSILON) * GEN->um )   
	   || (xfx < (1.+UNUR_EPSILON) * GEN->vl) 
	   || (xfx > (1.+UNUR_EPSILON) * GEN->vr) )
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
      if ( (gen->variant & SROU_VARFLAG_SQUEEZE) &&
	   (X >= GEN->xl) && 
	   (X <= GEN->xr ) && 
	   (U < GEN->um) ) {
	xx = xfx / (GEN->um - sfx);
	if ( (xx > (1.-UNUR_EPSILON) * GEN->xl) &&
	     (xx < (1.-UNUR_EPSILON) * GEN->xr) )
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");
	xx = V / (GEN->um - U);
	if ( (xx >= GEN->xl) && (xx <= GEN->xr ) )
	  return x;
      }
      if (U*U <= PDF(x))
	return x;
    }
  }
} 
double
_unur_gsrou_sample( struct unur_gen *gen )
{ 
  double U,Ur,V,W,X,Z;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,UNUR_INFINITY);
  while (1) {
    W = GEN->log_ab *_unur_call_urng(gen->urng);
    Z = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);
    U = (exp(-W)-1.) * GEN->a/GEN->b;
    V = -Z/(GEN->a + GEN->b*U);
    U *= GEN->um;
    Ur = pow(U,GEN->r);
    X = V/Ur + DISTR.mode;
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    if (Ur*U <= PDF(X))
      return X;
  }
} 
double
_unur_gsrou_sample_check( struct unur_gen *gen )
{ 
  double U,Ur,V,W,X,x,Z;
  double fx,uf,vf,vhl,vhr;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_SROU_GEN,UNUR_INFINITY);
  while (1) {
    W = GEN->log_ab *_unur_call_urng(gen->urng);
    Z = GEN->vl + _unur_call_urng(gen->urng) * (GEN->vr - GEN->vl);
    U = (exp(-W)-1.) * GEN->a/GEN->b;
    V = -Z/(GEN->a + GEN->b*U);
    U *= GEN->um;
    Ur = pow(U,GEN->r);
    X = V/Ur;
    x = X + DISTR.mode;
    if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
      continue;
    fx = PDF(x);
    uf = pow(fx,1./(GEN->r+1));
    vf = X * pow(fx,GEN->r/(GEN->r+1.));
    vhl = - GEN->vl /(GEN->a + GEN->b*(uf/GEN->um));
    vhr = - GEN->vr /(GEN->a + GEN->b*(uf/GEN->um));
    if ( ( uf > (1.+DBL_EPSILON) * GEN->um )   
	 || (vf < (1.+UNUR_EPSILON) * vhl )
	 || (vf > (1.+UNUR_EPSILON) * vhr ) )
      {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
      }
    if (Ur*U <= fx)
      return x;
  }
} 
int
_unur_srou_rectangle( struct unur_gen *gen )
{ 
  double vm, fm;             
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_SROU_GEN, UNUR_ERR_COOKIE );
  if (!(gen->set & SROU_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    if (fm <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return UNUR_ERR_GEN_DATA;
    }
    if (!_unur_isfinite(fm)) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
      return UNUR_ERR_PAR_SET;
    }
    GEN->um = sqrt(fm);    
  }
  vm = DISTR.area / GEN->um;
  if (gen->set & SROU_SET_CDFMODE) {
    GEN->vl = -GEN->Fmode * vm;
    GEN->vr = vm + GEN->vl;
    GEN->xl = GEN->vl/GEN->um;
    GEN->xr = GEN->vr/GEN->um;
  }
  else {
    GEN->vl = -vm;
    GEN->vr = vm;
    GEN->xl = GEN->vl/GEN->um;
    GEN->xr = GEN->vr/GEN->um;
    gen->variant &= ~SROU_VARFLAG_SQUEEZE;
  }
  return UNUR_SUCCESS;
} 
int
_unur_gsrou_envelope( struct unur_gen *gen )
{ 
  double fm;             
  double vm;             
  double pr;             
  double p;              
  double r = GEN->r;
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen, CK_SROU_GEN, UNUR_ERR_COOKIE );
  if (!(gen->set & SROU_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    if (fm <= 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return UNUR_ERR_GEN_DATA;
    }
    if (!_unur_isfinite(fm)) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
      return UNUR_ERR_PAR_SET;
    }
    GEN->um = pow(fm,1./(r+1.));    
  }
  vm = DISTR.area / (GEN->r*GEN->um);
  if (gen->set & SROU_SET_CDFMODE) {
    GEN->vl = -GEN->Fmode * vm;
    GEN->vr = vm + GEN->vl;
  }
  else {
    GEN->vl = -vm;
    GEN->vr = vm;
  }
  GEN->p = p = 1. - 2.187/pow(r + 5 - 1.28/r, 0.9460 );
  pr = pow(p,r);
  GEN->b = (1. - r * pr/p + (r-1.)*pr) / ((pr-1.)*(pr-1));
  GEN->a = -(p-1.)/(pr-1.) - p * GEN->b;
  GEN->log_ab = log(GEN->a/(GEN->a+GEN->b));
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_srou_debug_init( const struct unur_gen *gen, int is_reinit )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_SROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
    fprintf(LOG,"%s: method  = srou (simple universal ratio-of-uniforms)\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: reinit!\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  if (gen->set & SROU_SET_R) {
    fprintf(LOG,"%s: Generalized version: r = %g\n",gen->genid,GEN->r);
    fprintf(LOG,"%s:\n",gen->genid);
    fprintf(LOG,"%s: sampling routine = _unur_gsrou_sample",gen->genid);
    if (gen->variant & SROU_VARFLAG_VERIFY)
      fprintf(LOG,"_check");
    fprintf(LOG,"()\n%s:\n",gen->genid);
  }
  else {
    fprintf(LOG,"%s: Simple version (r = 1)  [default]\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    fprintf(LOG,"%s: sampling routine = _unur_srou_sample",gen->genid);
    if (gen->variant & SROU_VARFLAG_VERIFY)
      fprintf(LOG,"_check");
    else if (gen->variant & SROU_VARFLAG_MIRROR)
      fprintf(LOG,"_mirror");
    fprintf(LOG,"()\n%s:\n",gen->genid);
  }
  if (gen->set & SROU_SET_CDFMODE)
    fprintf(LOG,"%s: F(mode) = %g\n",gen->genid,GEN->Fmode);
  else
    fprintf(LOG,"%s: F(mode) unknown\n",gen->genid);
  if (gen->variant & SROU_VARFLAG_SQUEEZE)
    fprintf(LOG,"%s: use universal squeeze\n",gen->genid);
  else
    fprintf(LOG,"%s: no (universal) squeeze\n",gen->genid);
  if (gen->variant & SROU_VARFLAG_MIRROR)
    fprintf(LOG,"%s: use mirror principle\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->set & SROU_SET_R) {
    fprintf(LOG,"%s: Enveloping region:\n",gen->genid);
    fprintf(LOG,"%s:    um = %g\n",gen->genid,GEN->um);
    fprintf(LOG,"%s:    vl = %g\n",gen->genid,GEN->vl);
    fprintf(LOG,"%s:    vr = %g\n",gen->genid,GEN->vr);
    fprintf(LOG,"%s:    p  = %g\n",gen->genid,GEN->p);
    fprintf(LOG,"%s:    a  = %g\n",gen->genid,GEN->a);
    fprintf(LOG,"%s:    b  = %g\n",gen->genid,GEN->b);
  }
  else {
    fprintf(LOG,"%s: Rectangle:\n",gen->genid);
    fprintf(LOG,"%s:    left upper point  = (%g,%g)\n",gen->genid,GEN->vl,GEN->um);
    fprintf(LOG,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN->vr,GEN->um);
  }
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_srou_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  double h_area, rc;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %g   %s\n", DISTR.mode,
		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   area(PDF) = %g\n", DISTR.area);
  if (gen->set & SROU_SET_CDFMODE)
    _unur_string_append(info,"   F(mode)   = %g\n", GEN->Fmode); 
  else
    _unur_string_append(info,"   F(mode)   = [unknown]\n"); 
  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
			  "You may provide the \"mode\"");
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: SROU (Simple Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"   r = %g  %s\n", GEN->r,
		      (gen->set & SROU_SET_R) ? "[generalized version]" : "");
  if (gen->set & SROU_SET_CDFMODE)
    _unur_string_append(info,"   use CDF at mode\n");
  if (gen->variant & SROU_VARFLAG_SQUEEZE)
    _unur_string_append(info,"   use squeeze\n");
  if (gen->variant & SROU_VARFLAG_MIRROR)
    _unur_string_append(info,"   use mirror principle\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  if (gen->set & SROU_SET_R) {
    rc = unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize);
    _unur_string_append(info,"   enveloping rectangle = (%g,%g) x (%g,%g)\n",
			GEN->vl,GEN->vr, 0.,GEN->um);
    _unur_string_append(info,"   rejection constant = %.2f  [approx.]\n", rc);
  }
  else {
    _unur_string_append(info,"   bounding rectangle = (%g,%g) x (%g,%g)\n",
			GEN->vl,GEN->vr, 0.,GEN->um);
    h_area = (GEN->vr - GEN->vl) * GEN->um;
    _unur_string_append(info,"   area(hat) = %g\n", h_area);
    if (gen->set & SROU_SET_CDFMODE) 
      rc = 2.;
    else 
      rc = (gen->variant & SROU_VARFLAG_MIRROR) ? 2.829 : 4.;
    _unur_string_append(info,"   rejection constant = %g\n", rc);
  }
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"     r = %g  %s\n", GEN->r,
			(gen->set & SROU_SET_R) ? "" : "[default]");
    if (gen->set & SROU_SET_CDFMODE)
      _unur_string_append(info,"   cdfatmode = %g\n", GEN->Fmode); 
    else
      _unur_string_append(info,"   cdfatmode = [not set]\n"); 
    if (gen->variant & SROU_VARFLAG_SQUEEZE)
      _unur_string_append(info,"   usesqueeze\n");
    if (gen->variant & SROU_VARFLAG_MIRROR)
      _unur_string_append(info,"   usemirror\n");
    if (gen->variant & SROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & SROU_SET_CDFMODE))
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"cdfatmode\" to reduce the rejection constant.");
    _unur_string_append(info,"\n");
  }
} 
#endif   
