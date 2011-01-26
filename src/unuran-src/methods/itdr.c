/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <utils/unur_fp_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "itdr.h"
#include "itdr_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define C_MAX  (-0.1)
#define DX (1.e-6)      
#define NEAR_POLE  (1.e-8)
#define RESOLUTION_XI (1.e-5)
#define ITDR_VARFLAG_VERIFY   0x001u   
#define ITDR_DEBUG_REINIT    0x00000010u   
#define ITDR_SET_XI      0x001u     
#define ITDR_SET_CP      0x002u     
#define ITDR_SET_CT      0x004u     
#define GENTYPE "ITDR"         
static struct unur_gen *_unur_itdr_init( struct unur_par *par );
static int _unur_itdr_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_itdr_create( struct unur_par *par );
static int _unur_itdr_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_itdr_clone( const struct unur_gen *gen );
static void _unur_itdr_free( struct unur_gen *gen);
static double _unur_itdr_sample( struct unur_gen *gen );
static double _unur_itdr_sample_check( struct unur_gen *gen );
static int _unur_itdr_get_hat( struct unur_gen *gen );
static int _unur_itdr_get_hat_pole( struct unur_gen *gen );
static int _unur_itdr_get_hat_tail( struct unur_gen *gen );
static double _unur_itdr_lc( struct unur_gen *gen, double x );
static double _unur_itdr_ilc( struct unur_gen *gen, double x );
static double _unur_itdr_find_xt( struct unur_gen *gen, double b );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_itdr_debug_init( const struct unur_gen *gen, int error );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_itdr_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_itdr_par*)par->datap) 
#define GEN       ((struct unur_itdr_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define I2O(x)    ( GEN->sign*(x) + GEN->pole )
#define O2I(x)    ( GEN->sign * ((x)-GEN->pole) )
#define PDF(x)    ( _unur_cont_PDF(I2O(x), gen->distr) )
#define dPDF(x)   ( GEN->sign * _unur_cont_dPDF(I2O(x), gen->distr) )
#define PDFo(x)   ( _unur_cont_PDF((x), gen->distr) )
#define dPDFo(x)  ( _unur_cont_dPDF((x), gen->distr) )
#define logPDF(x)   ( _unur_cont_logPDF(I2O(x), gen->distr) )   
#define dlogPDF(x)  ( GEN->sign * _unur_cont_dlogPDF(I2O(x), gen->distr) )  
#define logPDFo(x)  ( _unur_cont_logPDF((x), gen->distr) )  
#define dlogPDFo(x) ( _unur_cont_dlogPDF((x), gen->distr) )  
#define T(c,x)   ( -pow((x), (c)) )
#define DT(c,x)  ( -(c)*pow((x), ((c)-1.)) )
#define TI(c,x)  ( pow(-(x), 1./(c)) )
#define FT(c,x)  ( -pow(-(x), ((c)+1.)/(c))*((c)/((c)+1.)) )
#define FTI(c,x) ( -pow(-(x)*((c)+1.)/(c), (c)/((c)+1.)) )
#define logTI(c,x)  ( -log(-(x)) / (c) )
#define TsI(c,x)  ( 1./((x)*(x)) )
#define FTs(c,x)  ( -1./(x) )
#define FTsI(c,x) ( -1./(x) )
#define _unur_itdr_getSAMPLE(gen) \
   ( ((gen)->variant & ITDR_VARFLAG_VERIFY) \
     ? _unur_itdr_sample_check : _unur_itdr_sample )
struct unur_par *
unur_itdr_new( const struct unur_distr *distr )
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
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dPDF");
    return NULL;
  }
  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode (pole)");
    return NULL; 
  }
  if ( ! (_unur_isfinite(DISTR_IN.mode) &&
	  (_unur_FP_equal(DISTR_IN.mode,DISTR_IN.domain[0]) ||
	   _unur_FP_equal(DISTR_IN.mode,DISTR_IN.domain[1]))) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_PROP,"pole not on boundary of domain");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_itdr_par) );
  COOKIE_SET(par,CK_ITDR_PAR);
  par->distr    = distr;          
  PAR->xi = INFINITY;       
  PAR->cp = INFINITY;       
  PAR->ct = INFINITY;       
  par->method   = UNUR_METH_ITDR;     
  par->variant  = 0u;                 
  par->set      = 0u;                 
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_itdr_init;
  return par;
} 
int
unur_itdr_set_xi( struct unur_par *par, double xi )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );
  if (xi <= par->distr->data.cont.BD_LEFT || 
      xi >= par->distr->data.cont.BD_RIGHT) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"xi out of domain");
    return UNUR_ERR_PAR_SET;
  }
  PAR->xi = xi;
  par->set |= ITDR_SET_XI;
  return UNUR_SUCCESS;
} 
int
unur_itdr_set_cp( struct unur_par *par, double cp )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );
  if ( cp > C_MAX || cp <= -1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp > -0.1 or <= -1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->cp = cp;
  par->set |= ITDR_SET_CP;
  return UNUR_SUCCESS;
} 
int
unur_itdr_set_ct( struct unur_par *par, double ct )
{
  double range;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );
  range = ( par->distr->data.cont.BD_RIGHT
	    - par->distr->data.cont.BD_LEFT );
  if ( ct > C_MAX || (ct <= -1. && !_unur_isfinite(range)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ct > -0.1 or <= -1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->ct = ct;
  par->set |= ITDR_SET_CT;
  return UNUR_SUCCESS;
} 
double unur_itdr_get_xi( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );
  return GEN->xi;
} 
double unur_itdr_get_cp( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );
  return GEN->cp;
} 
double unur_itdr_get_ct( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );
  return GEN->ct;
} 
double unur_itdr_get_area( struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );
  return GEN->Atot;
} 
int
unur_itdr_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );
  par->variant = (verify) 
    ? (par->variant | ITDR_VARFLAG_VERIFY) 
    : (par->variant & (~ITDR_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_itdr_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ITDR, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= ITDR_VARFLAG_VERIFY;
  else
    gen->variant &= ~ITDR_VARFLAG_VERIFY;
  SAMPLE = _unur_itdr_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_itdr_init( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_ITDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_ITDR_PAR,NULL);
  gen = _unur_itdr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_itdr_check_par(gen) != UNUR_SUCCESS) {
    _unur_itdr_free(gen); return NULL;
  }
  if (_unur_itdr_get_hat(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_itdr_debug_init(gen,UNUR_FAILURE);
#endif
    _unur_itdr_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_itdr_debug_init(gen,UNUR_SUCCESS);
#endif
  return gen;
} 
int
_unur_itdr_reinit( struct unur_gen *gen )
{
  int rcode;
  gen->set &= ~(ITDR_SET_XI | ITDR_SET_CP | ITDR_SET_CT);
  if ( (rcode = _unur_itdr_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  rcode = _unur_itdr_get_hat(gen);
  SAMPLE = _unur_itdr_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & ITDR_DEBUG_REINIT)
      _unur_itdr_debug_init(gen,rcode);
#endif
  return rcode;
} 
struct unur_gen *
_unur_itdr_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_ITDR_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_itdr_gen) );
  COOKIE_SET(gen,CK_ITDR_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_itdr_getSAMPLE(gen);
  gen->destroy = _unur_itdr_free;
  gen->clone = _unur_itdr_clone;
  gen->reinit = _unur_itdr_reinit;
  GEN->pole = DISTR.mode;  
  GEN->xi = PAR->xi;       
  GEN->cp = PAR->cp;       
  GEN->ct = PAR->ct;       
  GEN->bx = INFINITY;      
  GEN->xp = INFINITY;      
  GEN->xt = INFINITY;      
  GEN->alphap = INFINITY;  
  GEN->betap = INFINITY;
  GEN->Tfxt = INFINITY;    
  GEN->dTfxt = INFINITY;   
  GEN->by = INFINITY;      
  GEN->Ap = INFINITY;           
  GEN->Ac = INFINITY;           
  GEN->At = INFINITY;           
  GEN->Atot = INFINITY;    
  GEN->sy = 0.;            
  GEN->sign = 1.;          
  GEN->bd_right = INFINITY; 
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_itdr_info;
#endif
  return gen;
} 
int
_unur_itdr_check_par( struct unur_gen *gen )
{
  GEN->pole = DISTR.mode;  
  do {
    if (_unur_isfinite(DISTR.BD_LEFT) && !_unur_isfinite(DISTR.BD_RIGHT)) {
      GEN->sign = 1.; 
      if (dPDFo(DISTR.BD_LEFT) <= 0.) break;
    }
    if (!_unur_isfinite(DISTR.BD_LEFT) && _unur_isfinite(DISTR.BD_RIGHT)) {
      GEN->sign = -1.;
      if (dPDFo(DISTR.BD_RIGHT) >= 0.) break;
    }
    if (_unur_isfinite(DISTR.BD_LEFT) && _unur_isfinite(DISTR.BD_RIGHT)) {
      GEN->sign = (PDFo(DISTR.BD_LEFT)>=PDFo(DISTR.BD_RIGHT)) ? 1. : -1.;
      if ( GEN->sign*dPDFo(DISTR.BD_LEFT) <= 0. &&
	   GEN->sign*dPDFo(DISTR.BD_RIGHT) <= 0. )
	break;
    }
    _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute sign of region");
    return UNUR_ERR_DISTR_PROP;
  } while (1);
  GEN->bd_right = ( (GEN->sign > 0) 
		    ? DISTR.BD_RIGHT - GEN->pole
		    : GEN->pole - DISTR.BD_LEFT );
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_itdr_clone( const struct unur_gen *gen )
{
#define CLONE  ((struct unur_itdr_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ITDR_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_itdr_free( struct unur_gen *gen )
{
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_ITDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_ITDR_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_itdr_sample( struct unur_gen *gen )
{
  double U, V, X, Y;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ITDR_GEN,INFINITY);
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN->Atot;
    if (U < GEN->Ap) {
      V = _unur_call_urng(gen->urng) * GEN->Ap;
      if (_unur_isfsame(GEN->cp, -0.5)) {
	Y = ( FTsI(GEN->cp, GEN->betap*V + FTs(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TsI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
      else {
	Y = ( FTI(GEN->cp, GEN->betap*V + FT(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
    }
    else if ((U -= GEN->Ap) < GEN->Ac) {
      X = U * GEN->bx / GEN->Ac;
      Y = _unur_call_urng(gen->urng) * GEN->by;
      if (Y <= GEN->sy)
	return (I2O(X));
    }
    else {
      U -= GEN->Ac;
      if (_unur_isfsame(GEN->ct, -0.5)) {
	X = GEN->xt + (FTsI(GEN->ct,
			   GEN->dTfxt*U
			   + FTs(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TsI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
      else {
	X = GEN->xt + (FTI(GEN->ct,
			   GEN->dTfxt*U
			   + FT(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
    }
    X = I2O(X);
    if (Y <= PDFo(X))
      return X;
  }
} 
double
_unur_itdr_sample_check( struct unur_gen *gen )
{
#define ht(x)  ( TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*((x)-GEN->xt)) )
#define hp(x)  ( (T(GEN->cp,(x)) - GEN->alphap) / GEN->betap )
  double U, V, X, Y;
  double fx, hx, sqx;  
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ITDR_GEN,INFINITY);
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN->Atot;
    if (U < GEN->Ap) {
      V = _unur_call_urng(gen->urng) * GEN->Ap;
      if (_unur_isfsame(GEN->cp, -0.5)) {
	Y = ( FTsI(GEN->cp, GEN->betap*V + FTs(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TsI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
      else {
	Y = ( FTI(GEN->cp, GEN->betap*V + FT(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
      hx = hp(X);
      sqx = 0.;
    }
    else if ((U -= GEN->Ap) < GEN->Ac) {
      X = U * GEN->bx / GEN->Ac;
      Y = _unur_call_urng(gen->urng) * GEN->by;
      hx = hp(X);
      sqx = GEN->sy;
    }
    else {
      U -= GEN->Ac;
      if (_unur_isfsame(GEN->ct, -0.5)) {
	X = GEN->xt + (FTsI(GEN->ct,
			   GEN->dTfxt*U
			   + FTs(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TsI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
      else {
	X = GEN->xt + (FTI(GEN->ct,
			   GEN->dTfxt*U
			   + FT(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
      hx = ht(X);
      sqx = 0.;
    }
    X = I2O(X);
    fx = PDFo(X);
    if ( (1.+UNUR_EPSILON) * hx < fx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    if ( (1.-UNUR_EPSILON) * sqx > fx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");
    if (Y <= PDFo(X))
      return X;
  }
#undef ht
#undef hp
} 
int
_unur_itdr_get_hat( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ITDR_GEN,UNUR_ERR_COOKIE);
  if (gen->set & ITDR_SET_XI) {
    GEN->bx = O2I(GEN->xi);
  }
  else {
    GEN->bx = _unur_itdr_find_xt( gen, 0. );
    GEN->xi = I2O(GEN->bx);
    if (!_unur_isfinite(GEN->bx)) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute bx");
      return UNUR_ERR_DISTR_PROP;
    }
  }
  if (_unur_itdr_get_hat_pole(gen) != UNUR_SUCCESS)
    return UNUR_ERR_DISTR_PROP;
  if (_unur_FP_equal(GEN->bx, GEN->bd_right)) {
    GEN->At = 0.;
  }
  else {
    if (_unur_itdr_get_hat_tail(gen) != UNUR_SUCCESS)
      return UNUR_ERR_DISTR_PROP;
  }
  GEN->Atot = GEN->Ap + GEN->Ac + GEN->At;
  return UNUR_SUCCESS;
} 
int
_unur_itdr_get_hat_pole( struct unur_gen *gen )
{
#define hp(x)  ( (T(cp,(x)) - GEN->alphap) / GEN->betap )
  double cp, xp;
  double pdf_bx;
  double near_pole, ilc_near_pole, pdf_near_pole, logpdf_near_pole;
  double ilc_bx = -INFINITY;
  if (gen->set & ITDR_SET_CP) {
    cp = GEN->cp;
  }
  else {
    ilc_bx = _unur_itdr_ilc(gen, GEN->bx);
    near_pole = GEN->bx*NEAR_POLE + fabs(GEN->pole)*DBL_EPSILON;
    ilc_near_pole = (DISTR.logpdf) 
      ? logPDF(near_pole) / log(near_pole)
      : log(PDF(near_pole)) / log(near_pole);
    cp = ilc_near_pole;
    if (cp > C_MAX) cp = C_MAX;
    if (cp <= -1.) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->cp = cp;
  }
  if (cp < -0.5)
    GEN->bx = _unur_min(2.*GEN->bx, GEN->bd_right);
  pdf_bx = PDF(GEN->bx);
  near_pole = fabs(GEN->pole)*DBL_EPSILON;
  if (near_pole < 1.e-100) near_pole = 1.e-100;
  pdf_near_pole = logpdf_near_pole = INFINITY;
  while (1) {
    if (DISTR.logpdf) {
      logpdf_near_pole = logPDF(near_pole);
      if (_unur_isfinite(logpdf_near_pole)) 
	break;
    }
    else {
      pdf_near_pole = PDF(near_pole);
      if (_unur_isfinite(pdf_near_pole)) 
	break;
    }
    near_pole *= 1000.;
    if (!_unur_isfinite(near_pole)) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
      return UNUR_ERR_DISTR_PROP;
    }
  }
  while (1) {
    xp = GEN->bx * pow(1.+cp, -1./cp);
    if ( !(xp > 0. && xp < GEN->bx) ) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: xp");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->betap = DT(cp,xp) / dPDF(xp);
    GEN->alphap = T(cp,xp) - GEN->betap * PDF(xp);
    if ( hp(GEN->bx) < pdf_bx ||
	 (DISTR.logpdf && _unur_FP_less(log(hp(near_pole)), logpdf_near_pole)) ||
	 (DISTR.logpdf==NULL && _unur_FP_less(hp(near_pole), pdf_near_pole)) ) {
      if (gen->set & ITDR_SET_CP) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"inverse pdf not T_cp concave");
	return UNUR_ERR_DISTR_PROP;
      }
      GEN->cp = cp = 0.9*cp-0.1;
      if (cp < ilc_bx) {
	GEN->cp = cp = ilc_bx; 
	ilc_bx = -INFINITY;
      }
      if (cp < -0.999) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
	return UNUR_ERR_DISTR_PROP;
      }	
    }
    else break;
  }
  GEN->xp = xp;
  GEN->by = hp(GEN->bx);
  GEN->Ap = -FT(cp, GEN->alphap + GEN->betap * GEN->by) / GEN->betap;
  GEN->Ac = GEN->by * GEN->bx;
  GEN->sy = PDF(GEN->bx);
  return UNUR_SUCCESS;
#undef hp
} 
int
_unur_itdr_get_hat_tail( struct unur_gen *gen )
{
#define ht(x)  ( TI(ct, GEN->Tfxt + GEN->dTfxt*((x)-xt)) )
#define loght(x)  ( TI(ct, GEN->Tfxt + GEN->dTfxt*((x)-xt)) )
  double ct, xt;
  double lc_bx, lc_inf;
  double br;
  double bx = GEN->bx;
  GEN->xt = xt = _unur_itdr_find_xt( gen, bx );
  if (gen->set & ITDR_SET_CT) {
    ct = GEN->ct;
  }
  else {
    ct = _unur_itdr_lc(gen, 0.5*(bx + xt));
    if ( _unur_isfinite(GEN->bd_right)) 
      lc_inf = _unur_itdr_lc(gen, GEN->bd_right);
    else { 
      if (DISTR.logpdf) {
	lc_inf = log(1.e100) / logPDF(1.e100);
	lc_inf += -0.01;
      }
      else {
	lc_inf = log(1.e10*bx) / log(PDF(1.e10*bx));
	lc_inf += -0.05;
      }
    }
    if (lc_inf < ct) ct = lc_inf;
    if (ct > C_MAX) ct = C_MAX;
    if (ct <= -1.) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for tail: ct");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->ct = ct;
  }
  lc_bx = _unur_itdr_lc(gen, bx);
  while (1) {
    GEN->Tfxt = T(ct, PDF(xt));
    GEN->dTfxt = DT(ct, PDF(xt)) * dPDF(xt);
    br = 1000.*bx;  
    if (br > GEN->bd_right) br = GEN->bd_right;
    if ( ((GEN->Tfxt + GEN->dTfxt*(bx-xt)) >= 0.) ||
	 (DISTR.logpdf && (_unur_FP_less(loght(br),logPDF(br)) ||
			   _unur_FP_less(loght(bx), logPDF(bx)) )) ||
	 (DISTR.logpdf==NULL && (_unur_FP_less(ht(br),PDF(br)) ||
				 _unur_FP_less(ht(bx), PDF(bx)) )) ) {
      if (gen->set & ITDR_SET_CT) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"pdf not T_ct concave");
	return UNUR_ERR_DISTR_PROP;
      }
      ct = 0.5*(ct + lc_bx);
      if (ct > GEN->ct || ct < -0.999 || _unur_FP_approx(ct,lc_bx)) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for tail: ct");
	return UNUR_ERR_DISTR_PROP;
      }
      GEN->ct = ct;
    }
    else
      break;
  }
  GEN->At = (!_unur_isfinite(GEN->bd_right)) ? 0. 
    : FT(ct, GEN->Tfxt + GEN->dTfxt * (GEN->bd_right - xt)) / GEN->dTfxt;
  GEN->At += -FT(ct, GEN->Tfxt + GEN->dTfxt * (GEN->bx - xt)) / GEN->dTfxt;
  return UNUR_SUCCESS;
#undef ht
#undef loght
} 
double
_unur_itdr_lc( struct unur_gen *gen, double x )
{
  double dx, f, df, ddf;
  if (DISTR.dlogpdf == NULL) {
    f = PDF(x);
    df = dPDF(x);
    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;
    if (x+dx > GEN->bd_right)
      ddf = (dPDF(x)-dPDF(x-dx))/dx;
    else
      ddf = (dPDF(x+dx)-dPDF(x-dx))/(2.*dx);
    return 1. - ddf*f/(df*df); 
  }
  else {
    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;
    if (x+dx > GEN->bd_right)
      return (1./dlogPDF(x) - 1./dlogPDF(x-dx))/dx;
    else
      return (1./dlogPDF(x+dx) - 1./dlogPDF(x-dx))/(2.*dx);
  }
} 
double
_unur_itdr_ilc( struct unur_gen *gen, double x )
{
  if (DISTR.dlogpdf == NULL) {
    double dx, df, ddf;
    df = dPDF(x);
    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;
    if (x+dx > GEN->bd_right)
      ddf = (dPDF(x)-dPDF(x-dx))/dx;
    else
      ddf = (dPDF(x+dx)-dPDF(x-dx))/(2.*dx);
    return 1.+x*ddf/(df); 
  }
  else {
    double dx, dlf, ddlf;
    dlf = dlogPDF(x);
    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;
    if (x+dx > GEN->bd_right)
      ddlf = (dlogPDF(x)-dlogPDF(x-dx))/dx;
    else
      ddlf = (dlogPDF(x+dx)-dlogPDF(x-dx))/(2.*dx);
    return 1.+x*(dlf + ddlf/dlf); 
  }
} 
double
_unur_itdr_find_xt( struct unur_gen *gen, double b )
{
#define FKT(x) ( DISTR.dlogpdf \
                 ? (1./((x)-b) + dlogPDF(x)) \
                 : (((x)-b)*dPDF(x) + PDF(x)) )
  double xl, xu;  
  double xn;      
  if (b < 0.) return INFINITY;
  xl = b + _unur_max(1., (fabs(GEN->pole)+b)*UNUR_SQRT_DBL_EPSILON);
  if (xl > GEN->bd_right) xl = GEN->bd_right;
  while (!_unur_isfinite(FKT(xl)) || _unur_iszero(PDF(xl)) ) {
    xl = 0.5*(xl + b);
    if (!_unur_isfinite(xl) || _unur_FP_same(xl,b)) return INFINITY;
  }
  xu = xl;
  if (_unur_FP_greater(xu,GEN->bd_right)) return GEN->bd_right;
  if (FKT(xl)>0.) {
    do {
      xl = xu;
      xu += xu - b;
      if (!_unur_isfinite(xu) || xu < (1.+2.*DBL_EPSILON)*xl)
	return INFINITY;
      if (xu >= GEN->bd_right) 
	return GEN->bd_right;
    } while(FKT(xu) > 0.);
  }
  else { 
    do {
      xu = xl;
      xl = 0.5*(xl + b);
      if (!_unur_isfinite(xl)) return INFINITY;
    } while(FKT(xl) < 0.);
  }
  while(xu > (1.+RESOLUTION_XI)*xl) {
    xn = 0.5*(xl+xu);
    if(FKT(xn)>0.) 
      xl = xn;
    else 
      xu = xn;
  }
  return 0.5*(xl+xu);
#undef FKT
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_itdr_debug_init( const struct unur_gen *gen, int error )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ITDR_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = itdr (inverse transformed density rejection)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_itdr_sample",gen->genid);
  if (gen->variant & ITDR_VARFLAG_VERIFY) fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);
  fprintf(LOG,"%s: sign = %g\n",gen->genid, GEN->sign);
  fprintf(LOG,"%s: pole = %g\n",gen->genid, GEN->pole);
  fprintf(LOG,"%s: bd_right = %g\n",gen->genid, GEN->bd_right);
  fprintf(LOG,"%s: xi = %g",gen->genid, GEN->xi);
  fprintf(LOG,"%s\n", (gen->set & ITDR_SET_XI) ? "" : " [computed]");
  fprintf(LOG,"%s: bx = %g\n",gen->genid, GEN->bx);
  fprintf(LOG,"%s: pole region:\n",gen->genid);
  fprintf(LOG,"%s:\tcp = %g",gen->genid, GEN->cp);
  fprintf(LOG,"%s\n", (gen->set & ITDR_SET_CP) ? "" : " [computed]");
  fprintf(LOG,"%s:\txp = %g\n",gen->genid, GEN->xp);
  fprintf(LOG,"%s:\talphap = %g, betap = %g\n",gen->genid, GEN->alphap, GEN->betap);
  fprintf(LOG,"%s:\tby = %g\n",gen->genid, GEN->by);
  fprintf(LOG,"%s:\tsy = %g\n",gen->genid, GEN->sy);
  fprintf(LOG,"%s: tail region:\n",gen->genid);
  fprintf(LOG,"%s:\tct = %g",gen->genid, GEN->ct);
  fprintf(LOG,"%s\n", (gen->set & ITDR_SET_CT) ? "" : " [computed]");
  fprintf(LOG,"%s:\txt = %g\n",gen->genid, GEN->xt);
  fprintf(LOG,"%s:\tTfxt = %g, dTfxt = %g\n",gen->genid, GEN->Tfxt, GEN->dTfxt);
  fprintf(LOG,"%s: Area = %g + %g + %g = %g\n",gen->genid,
	  GEN->Ap, GEN->Ac, GEN->At, GEN->Atot);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: **** INIT %s ***\n",gen->genid,
	  (error==UNUR_SUCCESS) ? "successful" : "failed" );   
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_itdr_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   pole/mode = %g\n", DISTR.mode);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: ITDR (Inverse Transformed Density Rejection -- 2 point method)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g  [ = %g + %g + %g ]\n",
		      GEN->Atot, GEN->Ap, GEN->Ac, GEN->At);
  _unur_string_append(info,"   rejection constant = ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"%g\n", GEN->Atot/DISTR.area);
  else
    _unur_string_append(info,"%.2f  [approx. ]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   cp = %g  %s\n", GEN->cp,
			(gen->set & ITDR_SET_CP) ? "" : " [computed]");
    _unur_string_append(info,"   ct = %g  %s\n", GEN->cp,
			(gen->set & ITDR_SET_CT) ? "" : " [computed]");
    _unur_string_append(info,"   xi = %g  %s\n", GEN->xi,
			(gen->set & ITDR_SET_XI) ? "" : " [computed]");
    if (gen->variant & ITDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
