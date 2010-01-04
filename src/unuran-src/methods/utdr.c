/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "utdr.h"
#include "utdr_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define UTDR_VARFLAG_VERIFY     0x01u   
#define UTDR_DEBUG_REINIT    0x00000010u   
#define UTDR_SET_CPFACTOR       0x001u
#define UTDR_SET_DELTA          0x002u
#define UTDR_SET_PDFMODE        0x004u   
#define GENTYPE "UTDR"         
static struct unur_gen *_unur_utdr_init( struct unur_par *par );
static int _unur_utdr_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_utdr_create( struct unur_par *par );
static int _unur_utdr_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_utdr_clone( const struct unur_gen *gen );
static void _unur_utdr_free( struct unur_gen *gen);
static double _unur_utdr_sample( struct unur_gen *generator );
static double _unur_utdr_sample_check( struct unur_gen *generator );
static int _unur_utdr_hat( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_utdr_debug_init( const struct unur_gen *gen, 
				   double ttly, double ttlys, double ttry, double ttrys,
				   double cfac, int setupok, double c );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_utdr_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_utdr_par*)par->datap) 
#define GEN       ((struct unur_utdr_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    
#define _unur_utdr_getSAMPLE(gen) \
   ( ((gen)->variant & UTDR_VARFLAG_VERIFY) \
     ? _unur_utdr_sample_check : _unur_utdr_sample )
struct unur_par *
unur_utdr_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_utdr_par) );
  COOKIE_SET(par,CK_UTDR_PAR);
  par->distr       = distr;   
  PAR->c_factor     = 0.664; 
  PAR->delta_factor = 0.00001;
  PAR->fm        = -1.;                
  PAR->hm        = -1.;                
  par->method   = UNUR_METH_UTDR;     
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_utdr_init;
  return par;
} 
int 
unur_utdr_set_pdfatmode( UNUR_PAR *par, double fmode )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  PAR->fm = fmode;             
  PAR->hm = -1./sqrt(fmode);   
  par->set |= UTDR_SET_PDFMODE;
  return UNUR_SUCCESS;
} 
int
unur_utdr_set_cpfactor( struct unur_par *par, double cp_factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );
  if (cp_factor <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (cp_factor > 2.1)
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor > 2 not recommended. skip");
  PAR->c_factor = cp_factor;
  par->set |= UTDR_SET_CPFACTOR;
  return UNUR_SUCCESS;
} 
int
unur_utdr_set_deltafactor( struct unur_par *par, double delta )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );
  if (delta <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"delta <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (delta > 0.1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"delta must be small");
    return UNUR_ERR_PAR_SET;
  }
  PAR->delta_factor = delta;
  par->set |= UTDR_SET_DELTA;
  return UNUR_SUCCESS;
} 
int
unur_utdr_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, UTDR );
  par->variant = (verify) ? (par->variant | UTDR_VARFLAG_VERIFY) : (par->variant & (~UTDR_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_utdr_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= UTDR_VARFLAG_VERIFY;
  else
    gen->variant &= ~UTDR_VARFLAG_VERIFY;
  SAMPLE = _unur_utdr_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
int
unur_utdr_chg_pdfatmode( struct unur_gen *gen, double fmode )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, UTDR, UNUR_ERR_GEN_INVALID );
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"PDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  GEN->fm = fmode;             
  GEN->hm = -1./sqrt(fmode);   
  gen->set |= UTDR_SET_PDFMODE;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_utdr_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_UTDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);
  gen = _unur_utdr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_utdr_check_par(gen) != UNUR_SUCCESS) {
    _unur_utdr_free(gen); return NULL;
  }
  if ( _unur_utdr_hat(gen)!=UNUR_SUCCESS ) {
    _unur_utdr_free(gen); return NULL;
  }
  return gen;
} 
int
_unur_utdr_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_utdr_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  GEN->il = DISTR.BD_LEFT;
  GEN->ir = DISTR.BD_RIGHT;
  SAMPLE = _unur_utdr_getSAMPLE(gen);
  return _unur_utdr_hat( gen );
} 
struct unur_gen *
_unur_utdr_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_utdr_gen) );
  COOKIE_SET(gen,CK_UTDR_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_utdr_getSAMPLE(gen);
  gen->destroy = _unur_utdr_free;
  gen->clone = _unur_utdr_clone;
  gen->reinit = _unur_utdr_reinit;
  GEN->il = DISTR.BD_LEFT;           
  GEN->ir = DISTR.BD_RIGHT;          
  GEN->fm = PAR->fm;                  
  GEN->hm = PAR->hm;                  
  GEN->c_factor = PAR->c_factor;
  GEN->delta_factor = PAR->delta_factor;
  GEN->vollc = 0.; 
  GEN->volcompl = 0.; 
  GEN->voll = 0.; 
  GEN->al = 0.; 
  GEN->ar = 0.; 
  GEN->col = 0.; 
  GEN->cor = 0.; 
  GEN->sal = 0.; 
  GEN->sar = 0.; 
  GEN->bl = 0.; 
  GEN->br = 0.; 
  GEN->ttlx = 0.; 
  GEN->ttrx = 0.; 
  GEN->brblvolc = 0.; 
  GEN->drar = 0.; 
  GEN->dlal = 0.; 
  GEN->ooar2 = 0.; 
  GEN->ooal2 = 0.;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_utdr_info;
#endif
  return gen;
} 
int
_unur_utdr_check_par( struct unur_gen *gen )
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
_unur_utdr_clone( const struct unur_gen *gen )
{
#define CLONE  ((struct unur_utdr_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_UTDR_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_utdr_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_UTDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_UTDR_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
double
_unur_utdr_sample( struct unur_gen *gen )
{ 
  double u,v,x,help,linx;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_UTDR_GEN,INFINITY);
  while (1) {
    u = _unur_call_urng(gen->urng) * GEN->volcompl;
    if (u <= GEN->voll) {
      u = GEN->voll-u; 
      x = -GEN->dlal+GEN->ooal2/(u-GEN->col);
      help = GEN->al*(u-GEN->col);
      linx = help*help;
    }
    else {
      if (u <= GEN->vollc) {
	x = (u-GEN->voll) * GEN->brblvolc + GEN->bl;
	linx = GEN->fm;
      }
      else {
	x = - GEN->drar - GEN->ooar2 / (u-GEN->vollc - GEN->cor);
	help = GEN->ar * (u-GEN->vollc - GEN->cor);
	linx = help*help;
      }
    }
    v = _unur_call_urng(gen->urng) * linx;
    if (x<DISTR.mode) {
      if (x >= GEN->ttlx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sal;
	if (v * help * help <= 1.) return x;
      } 
    }
    else {
      if (x <= GEN->ttrx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sar;
	if (v * help * help <= 1.) return x; 
      }
    }
    if (v <= PDF(x)) return x; 
  }
} 
double
_unur_utdr_sample_check( struct unur_gen *gen )
{ 
  double u,v,x,help,linx,pdfx,squeezex;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_UTDR_GEN,INFINITY);
  while (1) {
    u = _unur_call_urng(gen->urng) * GEN->volcompl;
    if (u <= GEN->voll) {
      u = GEN->voll-u; 
      x = -GEN->dlal+GEN->ooal2/(u-GEN->col);
      help = GEN->al*(u-GEN->col);
      linx = help*help;
    }
    else {
      if (u <= GEN->vollc) {
	x = (u-GEN->voll) * GEN->brblvolc + GEN->bl;
	linx = GEN->fm;
      }
      else {
	x = - GEN->drar - GEN->ooar2 / (u-GEN->vollc - GEN->cor);
	help = GEN->ar * (u-GEN->vollc - GEN->cor);
	linx = help*help;
      }
    }
    v = _unur_call_urng(gen->urng) * linx;
    squeezex=0.;
    if (x<DISTR.mode) {
      if (x >= GEN->ttlx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sal;
        squeezex=1./(help*help);
      } 
    }
    else {
      if (x <= GEN->ttrx) {
	help = GEN->hm - (DISTR.mode - x) * GEN->sar;
        squeezex=1./(help*help);
      }
    }
    pdfx=PDF(x);
    if(_unur_FP_less(linx,pdfx))
      { _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
      _unur_log_printf(gen->genid,__FILE__,__LINE__,"x %e PDF(x) %e hat(x) %e squeeze(x) %e", \
		       x,pdfx,linx,squeezex ); 
      }
    if(_unur_FP_less(pdfx,squeezex))
      { _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");
      _unur_log_printf(gen->genid,__FILE__,__LINE__,"x %e PDF(x) %e hat(x) %e squeeze(x) %e", \
		       x,pdfx,linx,squeezex ); 
      }
    if (v <= PDF(x)) return x;
  }
} 
#define SMALL_VAL 1.e-50
int
_unur_utdr_hat( struct unur_gen *gen )
{ 
  double fm;
  int setupok=1;
  double c,cfac,volc,volr,ttly,ttlys,ttry,ttrys,dl,dr,delta,delta1,delta2,pdfx;
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen, CK_UTDR_GEN, UNUR_ERR_COOKIE );
  if (!(gen->set & UTDR_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    if (fm <= 0.) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return UNUR_ERR_GEN_DATA;
    }
    GEN->fm = fm;           
    GEN->hm = -1/sqrt(fm);  
  }
  ttry = 0.;
  ttrys = 0.;
  ttly = 0.;
  ttlys = 0.;
  dl = 0.;
  dr = 0.;
  volr = 0.;
  do {
    cfac = (setupok) ? GEN->c_factor : 2.;     
    c = cfac * DISTR.area/GEN->fm;
    setupok=1;         
    GEN->ttlx = DISTR.mode - c;
    GEN->ttrx = DISTR.mode + c;
    if ( GEN->ttlx < GEN->il) { 
      GEN->bl = GEN->il;
      GEN->al = 0.;
      GEN->voll = 0.;
      if (GEN->il < DISTR.mode) {
        GEN->ttlx = DISTR.mode + (GEN->il - DISTR.mode) * 0.6;
        pdfx=PDF(GEN->ttlx);
        if (pdfx > SMALL_VAL)
          GEN->sal = (GEN->hm + 1./sqrt(pdfx)) / (DISTR.mode - GEN->ttlx);
        else 
	  GEN->ttlx = DISTR.mode;
      }  
    }
    else {
     ttlys = PDF(GEN->ttlx);
     if (ttlys < SMALL_VAL) { 
       GEN->il = GEN->ttlx;
       GEN->bl = GEN->il;
       GEN->al = 0.;
       GEN->voll = 0.;
       GEN->ttlx=DISTR.mode;
     }
     else {
       ttlys = -1./sqrt(ttlys);
       GEN->sal =  (GEN->hm - ttlys) / (DISTR.mode - GEN->ttlx);
       delta2 = ( GEN->sal > 0. ) ? -ttlys/GEN->sal : -ttlys;
       delta1 = fabs(GEN->ttlx);
       delta = GEN->delta_factor * ((delta1<=delta2) ? delta2 : delta1);
       if (delta > c * 0.01) {
	 delta = UNUR_SQRT_DBL_EPSILON * ((delta1<=delta2) ? delta2 : delta1);
	 if (delta > c * 0.01) {
	   delta = c * 0.01;
	   _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,
			 "Delta larger than c/100!!, perhaps you can use a mode closer to 0 to remove this problem?");
         }
       }
       ttly = -1./sqrt(PDF(GEN->ttlx+delta));
       GEN->al = (ttly-ttlys)/delta;
       if (GEN->al <= 0.) 
	 setupok = 0; 
       else {
	 GEN->bl = GEN->ttlx + (GEN->hm - ttly)/GEN->al;
	 dl = ttly - GEN->al * GEN->ttlx;
	 GEN->voll = -1./(GEN->al * GEN->hm);
	 GEN->col = GEN->voll;
	 if (GEN->il > -INFINITY)
	   GEN->voll += 1./(GEN->al * (GEN->al * GEN->il + dl));
       }
     }
    }
    if(setupok) {
      if ( GEN->ttrx > GEN->ir) {
        GEN->br = GEN->ir;
        GEN->ar = 0.;
        volr = 0.;
        if (GEN->ir > DISTR.mode) {
          GEN->ttrx = DISTR.mode + (GEN->ir - DISTR.mode) * 0.6;
          pdfx = PDF(GEN->ttrx);
          if (pdfx > SMALL_VAL)
            GEN->sar = (GEN->hm + 1./sqrt(PDF(GEN->ttrx))) / (DISTR.mode - GEN->ttrx);
          else 
	    GEN->ttrx = DISTR.mode;
        } 
      }
      else {
        ttrys = PDF(GEN->ttrx);
        if (ttrys < SMALL_VAL){
          GEN->ir = GEN->ttrx;
          GEN->br = GEN->ir;
          GEN->ar = 0.;
          volr = 0.;
          GEN->ttrx = DISTR.mode;
	}
	else {
	  ttrys= -1./sqrt(ttrys);
	  GEN->sar = (GEN->hm - ttrys) / (DISTR.mode - GEN->ttrx);
	  delta2 = (GEN->sar<0.) ? ttrys/GEN->sar : -ttrys;
	  delta1 = fabs(GEN->ttrx);
	  delta = GEN->delta_factor * ((delta1<=delta2) ? delta2 : delta1);
	  if (delta > c*0.01) { 
	    delta = UNUR_SQRT_DBL_EPSILON * ((delta1<=delta2) ? delta2 : delta1);
	    if (delta > c*0.01) {
	      delta=c*0.01;
	      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
			    "Delta larger than c/100!!, perhaps you can use a mode closer to 0 to remove this problem?");
	    }
	  }
	  ttry = -1./sqrt(PDF(GEN->ttrx-delta));
	  GEN->ar = (ttrys - ttry)/delta;
	  if (GEN->ar >= 0.) 
	    setupok = 0;
	  else { 
	    GEN->br = GEN->ttrx + (GEN->hm - ttry) / GEN->ar;
	    dr = ttry - GEN->ar * GEN->ttrx;
	    volr = 1./(GEN->ar * GEN->hm);
	    GEN->cor = volr;
	    if (GEN->ir<INFINITY)
	      volr -= 1./(GEN->ar * (GEN->ar * GEN->ir + dr));
	  }
	}
      }
    }
    if(setupok) {
      volc = (GEN->br - GEN->bl) * GEN->fm;
      GEN->vollc = GEN->voll + volc;
      GEN->volcompl = GEN->vollc + volr;
      if (volc>0.) 
        GEN->brblvolc = (GEN->br - GEN->bl)/volc;
      if (!_unur_iszero(GEN->ar)) {
        GEN->drar = dr/GEN->ar;
        GEN->ooar2 = 1./(GEN->ar*GEN->ar);
      }
      if (!_unur_iszero(GEN->al)) {
        GEN->dlal = dl/GEN->al;
        GEN->ooal2 = 1./(GEN->al*GEN->al);
      }
    }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_utdr_debug_init(gen,ttly,ttlys,ttry,ttrys,cfac,setupok,c);
#endif
    if (!_unur_isfsame(cfac,2.)) {
      if(setupok)
        if (GEN->volcompl > 4. * DISTR.area || GEN->volcompl < 0.5 * DISTR.area)
        setupok=0;
    }
    else { 
      if (setupok==0 || GEN->volcompl > 8. * DISTR.area || GEN->volcompl < 0.5 * DISTR.area) {
        _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"; Area below hat too large or zero!! possible reasons: PDF, mode or area below PDF wrong;  density not T-concave\n");
        return 0;
      }
    }
  } while (!setupok);
  return UNUR_SUCCESS;
} 
#undef SMALL_VAL
#ifdef UNUR_ENABLE_LOGGING
static void
_unur_utdr_debug_init( const struct unur_gen *gen,
		       double ttly, double ttlys, double ttry, double ttrys,
		       double cfac, int setupok, double c )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_UTDR_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = utdr(transformed density rejection with 3 points of contact)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_utdr_sample",gen->genid);
  if (gen->variant & UTDR_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: Data for hat and squeeze:\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:\tc_factor=%e delta_factor=%e real c=%e\n",gen->genid,GEN->c_factor,GEN->delta_factor,c);
  fprintf(LOG,"%s:\ttlx=%e bl=%e mode=%e\n",gen->genid,GEN->ttlx,GEN->bl,DISTR.mode);
  fprintf(LOG,"%s:\tbr=%e trx=%e\n",gen->genid,GEN->br,GEN->ttrx);
  fprintf(LOG,"%s:\ttly=%e tlys=%e al=%e \n",gen->genid,ttly,ttlys,GEN->al);
  fprintf(LOG,"%s:\ttry=%e trys=%e ar=%e \n",gen->genid,ttry,ttrys,GEN->ar);
  fprintf(LOG,"%s:\tcfac=%e setupok=%d volcompl=%e pdf_area=%e\n",gen->genid,cfac,setupok,GEN->volcompl,DISTR.area);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_utdr_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %g   %s\n", unur_distr_cont_get_mode(distr),
		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   area(PDF) = %g\n", DISTR.area);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: UTDR (Universal Transformed Density Rejection -- 3 point method)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   rejection constant = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   deltafactor = %g  %s\n", GEN->delta_factor,
			(gen->set & UTDR_SET_DELTA) ? "" : "[default]");
    if (gen->set & UTDR_SET_PDFMODE)
      _unur_string_append(info,"   pdfatmode = %g\n", GEN->fm);
    if (gen->set & UTDR_SET_CPFACTOR)
      _unur_string_append(info,"   cpfactor = %g\n", GEN->c_factor);
    if (gen->variant & UTDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
