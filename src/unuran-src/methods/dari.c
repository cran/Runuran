/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dari.h"
#include "dari_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define DARI_VARFLAG_VERIFY     0x01u   
#define DARI_DEBUG_REINIT    0x00000010u   
#define DARI_SET_CFACTOR        0x001u
#define DARI_SET_TABLESIZE      0x002u
#define GENTYPE "DARI"         
static struct unur_gen *_unur_dari_init( struct unur_par *par );
static int _unur_dari_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dari_create( struct unur_par *par );
static int _unur_dari_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_dari_clone( const struct unur_gen *gen );
static void _unur_dari_free( struct unur_gen *gen);
static int _unur_dari_sample( struct unur_gen *gen );
static int _unur_dari_sample_check( struct unur_gen *gen );
static int _unur_dari_hat( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dari_debug_init( struct unur_gen *gen, const char *status );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dari_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr     
#define PAR       ((struct unur_dari_par*)par->datap) 
#define GEN       ((struct unur_dari_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.discr          
#define PMF(x)    _unur_discr_PMF((x),(gen->distr))    
#define T(x) (-1./sqrt(x))    
#define F(x) (-1./(x))        
#define FM(x) (-1./(x))       
#define N0 (GEN->n[0])        
#define _unur_dari_getSAMPLE(gen) \
   ( ((gen)->variant & DARI_VARFLAG_VERIFY) \
     ? _unur_dari_sample_check : _unur_dari_sample )
struct unur_par *
unur_dari_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  if (DISTR_IN.pmf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF"); 
    return NULL;
  }
  if (DISTR_IN.domain[0] < 0) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_PROP,"domain contains negative numbers");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_dari_par) );
  COOKIE_SET(par,CK_DARI_PAR);
  par->distr       = distr;   
  PAR->c_factor  = 0.664;
  PAR->squeeze   = 0; 
  PAR->size      = 100; 
  par->method   = UNUR_METH_DARI;     
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_dari_init;
  return par;
} 
int
unur_dari_set_cpfactor( struct unur_par *par, double cpfactor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DARI );
  if (cpfactor <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (cpfactor > 2.1)
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor > 2 not recommended. skip");
  PAR->c_factor = cpfactor;
  par->set |= DARI_SET_CFACTOR;
  return UNUR_SUCCESS;
} 
int
unur_dari_set_squeeze( struct unur_par *par, int squeeze )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DARI );
  PAR->squeeze = squeeze;
  return UNUR_SUCCESS;
} 
int
unur_dari_set_tablesize( struct unur_par *par, int size )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DARI );
  if (size < 0) {  
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid table size");
    return UNUR_ERR_PAR_SET;
  }
  PAR->size = size;
  par->set |= DARI_SET_TABLESIZE;
  return UNUR_SUCCESS;
} 
int
unur_dari_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DARI );
  par->variant = (verify) ? (par->variant | DARI_VARFLAG_VERIFY) : (par->variant & (~DARI_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_dari_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DARI, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_discr_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= DARI_VARFLAG_VERIFY;
  else
    gen->variant &= ~DARI_VARFLAG_VERIFY;
  SAMPLE = _unur_dari_getSAMPLE(gen); 
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dari_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE, par, NULL );
  if ( par->method != UNUR_METH_DARI ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DARI_PAR,NULL);
  gen = _unur_dari_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_dari_check_par(gen) != UNUR_SUCCESS) {
    _unur_dari_free(gen); return NULL;
  }
  if ( _unur_dari_hat(gen)!=UNUR_SUCCESS ) {
    _unur_dari_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_dari_debug_init(gen,"INIT completed");
#endif
  return gen;
} 
int
_unur_dari_reinit( struct unur_gen *gen )
{
  int result;
  if ( (result = _unur_dari_check_par(gen)) != UNUR_SUCCESS)
    return result;
  if ( (result = _unur_dari_hat( gen )) != UNUR_SUCCESS)
    return result;
  SAMPLE = _unur_dari_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & DARI_DEBUG_REINIT)
    _unur_dari_debug_init(gen,"REINIT completed");
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dari_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DARI_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_dari_gen) );
  COOKIE_SET(gen,CK_DARI_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_dari_getSAMPLE(gen);
  gen->destroy = _unur_dari_free;
  gen->clone = _unur_dari_clone;
  gen->reinit = _unur_dari_reinit;
  GEN->squeeze = PAR->squeeze;        
  GEN->c_factor = PAR->c_factor;      
  if ((unsigned)DISTR.BD_RIGHT - (unsigned)DISTR.BD_LEFT < INT_MAX)
    GEN->size = _unur_min(PAR->size,DISTR.BD_RIGHT-DISTR.BD_LEFT+1);
  else 
    GEN->size = PAR->size;
  GEN->hp = (GEN->size > 0) ? _unur_xmalloc( GEN->size * sizeof(double) ) : NULL;
  GEN->hb = (GEN->size > 0) ? _unur_xmalloc( GEN->size * sizeof(char) )   : NULL;
  GEN->vt=0.;            
  GEN->vc=0.;            
  GEN->vcr=0.;           
  GEN->xsq[0]=0.;        
  GEN->xsq[1]=0.;        
  GEN->y[0]=0.;          
  GEN->y[1]=0.;          
  GEN->ys[0]=0.;         
  GEN->ys[1]=0.;         
  GEN->ac[0]=0.;         
  GEN->ac[1]=0.;         
  GEN->pm=0.;            
  GEN->Hat[0]=0.;        
  GEN->Hat[1]=0.;        
  GEN->m=0;              
  GEN->x[0]=0;           
  GEN->x[1]=0;           
  GEN->s[0]=0;           
  GEN->s[1]=0;           
  GEN->n[0]=0;           
  GEN->n[1]=0;           
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dari_info;
#endif
  return gen;
} 
int
_unur_dari_check_par( struct unur_gen *gen )
{
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (unur_distr_discr_upd_mode( gen->distr )!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }
  if (DISTR.BD_LEFT > DISTR.mode)
    DISTR.mode = DISTR.BD_LEFT;
  else if (DISTR.BD_RIGHT < DISTR.mode)
    DISTR.mode = DISTR.BD_RIGHT;
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (unur_distr_discr_upd_pmfsum(gen->distr)!=UNUR_SUCCESS)
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF; use default");
  if (DISTR.sum <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"sum <= 0");
    return UNUR_ERR_GEN_DATA;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dari_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_dari_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DARI_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  if (GEN->size > 0) {
    CLONE->hp = _unur_xmalloc( GEN->size * sizeof(double) );
    memcpy( CLONE->hp, GEN->hp, GEN->size * sizeof(double) );
    CLONE->hb = _unur_xmalloc( GEN->size * sizeof(char) );
    memcpy( CLONE->hb, GEN->hb, GEN->size * sizeof(char) );
  }
  return clone;
#undef CLONE
} 
void
_unur_dari_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_DARI ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DARI_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->hp)   free(GEN->hp);
  if (GEN->hb)   free(GEN->hb);
  _unur_generic_free(gen);
} 
int
_unur_dari_sample( struct unur_gen *gen )
{
  double U, h;
  double X = 0.;
  int k,i;
  static const int sign[2] = {-1,1};
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DARI_GEN,INT_MAX);
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN->vt;
    if (U<=GEN->vc) {
      X = U * (GEN->ac[1]-GEN->ac[0]) / GEN->vc + GEN->ac[0]; 
      k = (int)(X+0.5);
      i = (k<GEN->m) ? 0 : 1;
      if (GEN->squeeze && sign[i]*(GEN->ac[i]-GEN->s[i]) > sign[i]*(X-k))
	return k;
      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if (!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = 0.5 - PMF(k)/GEN->pm;
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
      }
      else {
	h = 0.5-PMF(k)/GEN->pm;
      }
      if (h <= sign[i]*(k-X))
	return k;
    }
    else {
      if (U<= GEN->vcr) {
	i = 1;
	U -= GEN->vc;
      } 
      else {
	i = 0;
	U -= GEN->vcr;
      }
      U = GEN->Hat[i] + sign[i]*U; 
      X = GEN->x[i] + (FM(U*GEN->ys[i])-GEN->y[i]) / GEN->ys[i];
      k = (int)(X+0.5);
      if (GEN->squeeze && (sign[i]*k <= sign[i]*GEN->x[i]+1) && (GEN->xsq[i] <= sign[i]*(X-k))) 
	return k;
      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if (!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i] - PMF(k);
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
      }
      else {
	h = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i]-PMF(k);
      }
      if (sign[i]*U >= h)
	return k;
    }
  }
} 
int
_unur_dari_sample_check( struct unur_gen *gen )
{
  double U, h;
  double X = 0.;
  double hkm05;
  int k,i;
  static const int sign[2] = {-1,1};
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DARI_GEN,INT_MAX);
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN->vt;
    if (U <= GEN->vc) {
      X = U * (GEN->ac[1]-GEN->ac[0]) / GEN->vc + GEN->ac[0]; 
      k = (int)(X+0.5);
      i = (k<GEN->m) ? 0 : 1;
      if (GEN->squeeze && sign[i]*(GEN->ac[i]-GEN->s[i]) > sign[i]*(X-k))
	return k;
      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if (!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = 0.5 - PMF(k)/GEN->pm;
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
	if (h+UNUR_EPSILON*100.<-0.5) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PMF(i) > hat(i) for centerpart");
	  _unur_log_printf(gen->genid,__FILE__,__LINE__,
			   "i %d PMF(x) %.20e hat(x) %.20e", k,PMF(k),GEN->pm ); 
        }
      }
      else {
	h = 0.5 - PMF(k)/GEN->pm;
	if (h+UNUR_EPSILON*100.<-0.5) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PMF(i) > hat(i) for centerpart");
	  _unur_log_printf(gen->genid,__FILE__,__LINE__,
			   "i %d PMF(x) %.20e hat(x) %.20e", k,PMF(k),GEN->pm ); 
        }
      }
      if (h <= sign[i]*(k-X))
	return k;
    }
    else {
      if (U<= GEN->vcr) {
	i = 1;
	U -= GEN->vc;
      } 
      else {
	i = 0;
	U -= GEN->vcr;
      }
      U = GEN->Hat[i] + sign[i]*U; 
      X = GEN->x[i] + (FM(U*GEN->ys[i])-GEN->y[i]) / GEN->ys[i];
      k = (int)(X+0.5);
      if(k==GEN->s[i]) 
	k += sign[i];
      if (GEN->squeeze && (sign[i]*k <= sign[i]*GEN->x[i]+1) && (GEN->xsq[i] <= sign[i]*(X-k))) 
	return k;
      if (sign[i]*k <= sign[i]*GEN->n[i]) {
	if(!GEN->hb[k-N0]) {
	  GEN->hp[k-N0] = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i] - PMF(k); 
          if(k != GEN->s[i]+sign[i]) {
            hkm05 = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k-sign[i]*0.5-GEN->x[i])) / GEN->ys[i];
	    if (GEN->hp[k-N0]+UNUR_EPSILON < hkm05) {
	      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
			  "for tailpart hat too low, ie hp[k] < H(k-0.5)");
	      _unur_log_printf(gen->genid,__FILE__,__LINE__,
			       "k %d hp  %.20e H(k-0.5) %.20e ", k,GEN->hp[k-N0],hkm05 ); 
            }
	  }
	  GEN->hb[k-N0] = 1;
	}
	h = GEN->hp[k-N0];
      }
      else {
	h = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k+sign[i]*0.5-GEN->x[i])) / GEN->ys[i] - PMF(k);
        hkm05 = sign[i] * F(GEN->y[i]+GEN->ys[i]*(k-sign[i]*0.5-GEN->x[i])) / GEN->ys[i];
        if(k != GEN->s[i]+sign[i]) {
	  if (h+UNUR_EPSILON < hkm05) {
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
			"PMF(i) > hat(i) for tailpart");
	    _unur_log_printf(gen->genid,__FILE__,__LINE__,
			     "k %d h  %.20e H(k-0.5) %.20e ", k,h,hkm05 ); 
	  }
        }
      }
      if (sign[i]*U >= h)
	return k;
    }
  }
} 
int
_unur_dari_hat( struct unur_gen *gen )
{
  int sign[2] = {-1,1};
  int b[2], d, i, j;
  double v[2], at[2];
  double t0 = 1.;
  int setup = 1;
  int rep = 1;
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen, CK_DARI_GEN, UNUR_ERR_COOKIE );
  GEN->m = DISTR.mode;
  b[0] = DISTR.BD_LEFT;
  b[1] = DISTR.BD_RIGHT;
  GEN->pm = PMF(GEN->m);
  d = _unur_max(2, (int)( GEN->c_factor/(GEN->pm/DISTR.sum)));
  if (_unur_iszero(GEN->pm)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PMF(mode)=0");
    return UNUR_ERR_GEN_DATA;
  }
  do {
    for(i=0; i<=1; i++) {
      GEN->x[i] = GEN->m + sign[i] * d;
      if (sign[i]*GEN->x[i]+1 > sign[i]*b[i]) {
	v[i] = 0; 
	GEN->s[i] = b[i];
      }
      else {
	GEN->y[i] = T( PMF(GEN->x[i]) );
	GEN->ys[i] = sign[i] * (T( PMF(GEN->x[i]+sign[i])) - GEN->y[i]);
	if (GEN->ys[i]*sign[i] > -DBL_EPSILON) {
	  setup = -setup; 
	  i = 1; 
	}
        else {
	  GEN->s[i] = (int)(0.5+GEN->x[i]+(T(GEN->pm)-GEN->y[i])/GEN->ys[i]);
	  GEN->Hat[i] = ( F(GEN->y[i]+GEN->ys[i]*(GEN->s[i]+sign[i]*1.5-GEN->x[i])) /
			  GEN->ys[i]-sign[i]*PMF(GEN->s[i]+sign[i]) ); 
	  at[i] = GEN->x[i] + (FM(GEN->ys[i]*GEN->Hat[i])-GEN->y[i]) / GEN->ys[i]; 
          if(GEN->squeeze)
	    GEN->xsq[i] = sign[i]*(at[i]-(GEN->s[i]+sign[i]));
	  v[i] = sign[i]*(F(GEN->y[i]+GEN->ys[i]*(b[i]+sign[i]*0.5-GEN->x[i]))/
			  GEN->ys[i]-F(GEN->y[i]+GEN->ys[i]*(at[i]-GEN->x[i]))/GEN->ys[i]);
	}
      }
      if (setup>0)
	GEN->ac[i] = GEN->s[i] + sign[i]*(PMF(GEN->s[i])/GEN->pm-0.5);
    }
    if(setup>0) {
      GEN->vc = GEN->pm*(GEN->ac[1]-GEN->ac[0]); 
      GEN->vt = GEN->vc+v[0]+v[1];
      GEN->vcr = GEN->vc+v[1];
      GEN->n[0] = _unur_max(b[0],GEN->m - GEN->size/2);
      GEN->n[1] = GEN->n[0] + GEN->size - 1;
      if (GEN->n[1] > b[1]) {
	GEN->n[1] = b[1];
	GEN->n[0] = GEN->n[1]- GEN->size + 1;
      }
      for (j=0; j<GEN->size; j++)
	GEN->hb[j] = 0;
    }
    if (setup == 1 || setup == -1) {
      t0= 2. * DISTR.sum;
      if (setup==1 && GEN->vt<=t0)
	rep=0;
      else { 
	setup = 2;
	d = (int) (t0 / GEN->pm);
      }
    }
    else 
      rep=0; 
  } while(rep);
  if (setup == -2 || GEN->vt > 100.*t0 || !(GEN->vt > 0.)) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug)
      _unur_dari_debug_init(gen,"RE/INIT failed try again");
#endif
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"Area below hat too large or zero!! possible reasons: PDF, mode or area below PMF wrong;  or PMF not T-concave");
    return UNUR_ERR_GEN_DATA;
  }
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_dari_debug_init( struct unur_gen *gen, const char *status )
{
  FILE *LOG;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DARI_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = dari (discrete automatic rejection inversion)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );
  fprintf(LOG,"%s: sampling routine = _unur_dari_sample",gen->genid);
  if (gen->variant & DARI_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: Data for hat and squeeze:\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:area below hat: total %f center: %f, left tail %f, right tail %f\n", gen->genid,
	  GEN->vt, GEN->vc, GEN->vt-GEN->vcr, GEN->vcr-GEN->vc);
  fprintf(LOG,"%s: mode %d and mode probability %f\n",gen->genid, GEN->m, GEN->pm); 
  for(i=0;i<=1;i++) {
    fprintf(LOG,"%s:i=%d: x=%d; Hat=%f; ac=%f; s=%d;\n", gen->genid,
	    i, GEN->x[i], GEN->Hat[i], GEN->ac[i], GEN->s[i]);
    fprintf(LOG,"%s:i=%d: xsq=%f; y=%f; ys=%f; n:=%d (for aux.table)\n", gen->genid,
	    i, GEN->xsq[i], GEN->y[i], GEN->ys[i], GEN->n[i]);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: %s ************\n",gen->genid, status );
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dari_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PMF\n");
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %d   %s\n", DISTR.mode,
                      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   sum(PMF)  = %g   %s\n", DISTR.sum,
                      (distr->set & UNUR_DISTR_SET_PMFSUM) ? "" : "[unknown]");
  _unur_string_append(info,"\n");
  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may provide the \"mode\".");
    if (!(distr->set & UNUR_DISTR_SET_PMFSUM))
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may provide the \"pmfsum\".");
    _unur_string_append(info,"\n");
  }
  _unur_string_append(info,"method: DARI (Discrete Automatic Rejection Inversion)\n");
  if (GEN->size == 0) 
    _unur_string_append(info,"   no table\n");
  else
    _unur_string_append(info,"   use table of size %d\n", GEN->size);
  if (GEN->squeeze)
    _unur_string_append(info,"   use squeeze\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   sum(hat) = %g\n",GEN->vt);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PMFSUM)
    _unur_string_append(info,"= %g\n", GEN->vt/DISTR.sum);
  else
    _unur_string_append(info,"= %.2f  [approx.]\n", 
			unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   tablesize = %d  %s\n", GEN->size,
			(gen->set & DARI_SET_TABLESIZE) ? "" : "[default]");
    if (GEN->squeeze)
      _unur_string_append(info,"   squeeze = on\n");
    if (gen->set & DARI_SET_CFACTOR)
      _unur_string_append(info,"   cpfactor = %g\n",   GEN->c_factor);
    if (gen->variant & DARI_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
