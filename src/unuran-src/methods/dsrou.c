/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dsrou.h"
#include "dsrou_struct.h"
#define DSROU_VARFLAG_VERIFY   0x002u  
#define DSROU_DEBUG_REINIT    0x00000010u  
#define DSROU_SET_CDFMODE     0x001u   
#define GENTYPE "DSROU"        
static struct unur_gen *_unur_dsrou_init( struct unur_par *par );
static int _unur_dsrou_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dsrou_create( struct unur_par *par );
static int _unur_dsrou_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_dsrou_clone( const struct unur_gen *gen );
static void _unur_dsrou_free( struct unur_gen *gen);
static int _unur_dsrou_sample( struct unur_gen *gen );
static int _unur_dsrou_sample_check( struct unur_gen *gen );
static int _unur_dsrou_rectangle( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dsrou_debug_init( const struct unur_gen *gen, int is_reinit );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dsrou_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr     
#define PAR       ((struct unur_dsrou_par*)par->datap) 
#define GEN       ((struct unur_dsrou_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.discr          
#define PMF(x)    _unur_discr_PMF((x),(gen->distr))   
#define _unur_dsrou_getSAMPLE(gen) \
   ( ((gen)->variant & DSROU_VARFLAG_VERIFY) \
     ? _unur_dsrou_sample_check : _unur_dsrou_sample )
#define SQRT2     (M_SQRT2)
struct unur_par *
unur_dsrou_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_dsrou_par) );
  COOKIE_SET(par,CK_DSROU_PAR);
  par->distr    = distr;      
  PAR->Fmode     = -1.;                
  par->method   = UNUR_METH_DSROU;    
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_dsrou_init;
  return par;
} 
int 
unur_dsrou_set_cdfatmode( struct unur_par *par, double Fmode )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DSROU );
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
   }
  PAR->Fmode = Fmode;
  par->set |= DSROU_SET_CDFMODE;
  return UNUR_SUCCESS;
} 
int
unur_dsrou_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DSROU );
  par->variant = (verify) ? (par->variant | DSROU_VARFLAG_VERIFY) : (par->variant & (~DSROU_VARFLAG_VERIFY));
  return UNUR_ERR_NULL;
} 
int
unur_dsrou_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_discr_error) 
    return UNUR_FAILURE;
  if (verify)
    gen->variant |= DSROU_VARFLAG_VERIFY;
  else
    gen->variant &= ~DSROU_VARFLAG_VERIFY;
  SAMPLE = _unur_dsrou_getSAMPLE(gen); 
  return UNUR_SUCCESS;
} 
int
unur_dsrou_chg_cdfatmode( struct unur_gen *gen, double Fmode )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  GEN->Fmode = Fmode;
  gen->set |= DSROU_SET_CDFMODE;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dsrou_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_DSROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DSROU_PAR,NULL);
  gen = _unur_dsrou_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_dsrou_check_par(gen) != UNUR_SUCCESS) {
    _unur_dsrou_free(gen); return NULL;
  }
  if ( _unur_dsrou_rectangle(gen)!=UNUR_SUCCESS ) {
    _unur_dsrou_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_dsrou_debug_init(gen, FALSE);
#endif
  return gen;
} 
int
_unur_dsrou_reinit( struct unur_gen *gen )
{
  int result;
  if ( (result = _unur_dsrou_check_par(gen)) != UNUR_SUCCESS)
    return result;
  if ( (result = _unur_dsrou_rectangle(gen)) != UNUR_SUCCESS)
    return result;
  SAMPLE = _unur_dsrou_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & DSROU_DEBUG_REINIT) _unur_dsrou_debug_init(gen,TRUE);
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dsrou_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSROU_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_dsrou_gen) );
  COOKIE_SET(gen,CK_DSROU_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_dsrou_getSAMPLE(gen);
  gen->destroy = _unur_dsrou_free;
  gen->clone = _unur_dsrou_clone;
  gen->reinit = _unur_dsrou_reinit;
  GEN->Fmode = PAR->Fmode;            
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dsrou_info;
#endif
  return gen;
} 
int
_unur_dsrou_check_par( struct unur_gen *gen )
{
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (unur_distr_discr_upd_mode(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (unur_distr_discr_upd_pmfsum(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF");
      return UNUR_ERR_DISTR_REQUIRED;
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
_unur_dsrou_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_dsrou_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DSROU_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  return clone;
#undef CLONE
} 
void
_unur_dsrou_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_DSROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DSROU_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_generic_free(gen);
} 
int
_unur_dsrou_sample( struct unur_gen *gen )
{ 
  double U,V,X;
  int I;
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DSROU_GEN,INT_MAX);
  while (1) {
    V = GEN->al + _unur_call_urng(gen->urng) * (GEN->ar - GEN->al);
    V /= (V<0.) ? GEN->ul : GEN->ur;    
    while ( _unur_iszero(U = _unur_call_urng(gen->urng)));
    U *= (V<0.) ? GEN->ul : GEN->ur;
    X = floor(V/U) + DISTR.mode;
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    I = (int) X;
    if (U*U <= PMF(I))
      return I;
  }
} 
int
_unur_dsrou_sample_check( struct unur_gen *gen )
{ 
  double U,V,pI,VI,X;
  double um2, vl, vr;
  int I;
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DSROU_GEN,INT_MAX);
  while (1) {
    V = GEN->al + _unur_call_urng(gen->urng) * (GEN->ar - GEN->al);
    V /= (V<0.) ? GEN->ul : GEN->ur;
    while ( _unur_iszero(U = _unur_call_urng(gen->urng)));
    U *= (V<0.) ? GEN->ul : GEN->ur;
    X = floor(V/U) + DISTR.mode;
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    I = (int) X;
    pI = PMF(I);
    VI = V/U * sqrt(pI);
    um2 = (2.+4.*DBL_EPSILON) * ((V<0) ? GEN->ul*GEN->ul : GEN->ur*GEN->ur);
    vl = (GEN->ul>0.) ? (1.+UNUR_EPSILON) * GEN->al/GEN->ul : 0.;
    vr = (1.+UNUR_EPSILON) * GEN->ar/GEN->ur;
    if ( pI > um2 || VI < vl || VI > vr ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PMF(x) > hat(x)");
    }
    if (U*U <= pI)
      return I;
  }
} 
int
_unur_dsrou_rectangle( struct unur_gen *gen )
{ 
  double pm, pbm;               
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_DSROU_GEN, UNUR_ERR_COOKIE );
  pm = PMF(DISTR.mode);
  pbm = (DISTR.mode-1 < DISTR.BD_LEFT) ? 0. : PMF(DISTR.mode-1);
  if (pm <= 0. || pbm < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PMF(mode) <= 0.");
    return UNUR_ERR_GEN_DATA;
  }
  GEN->ul = sqrt(pbm);
  GEN->ur = sqrt(pm);
  if (_unur_iszero(GEN->ul)) {
    GEN->al = 0.;
    GEN->ar = DISTR.sum;
  }
  else if (gen->set & DSROU_SET_CDFMODE) {
    GEN->al = -(GEN->Fmode * DISTR.sum)+pm;
    GEN->ar = DISTR.sum + GEN->al;
  }
  else {
    GEN->al = -(DISTR.sum - pm);
    GEN->ar = DISTR.sum;
  }    
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
static void
_unur_dsrou_debug_init( const struct unur_gen *gen, int is_reinit )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSROU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
    fprintf(LOG,"%s: method  = dsrou (discrete simple universal ratio-of-uniforms)\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: reinit!\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );
  fprintf(LOG,"%s: sampling routine = _unur_dsrou_sample",gen->genid);
  if (gen->variant & DSROU_VARFLAG_VERIFY)
    fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);
  if (gen->set & DSROU_SET_CDFMODE)
    fprintf(LOG,"%s: CDF(mode) = %g\n",gen->genid,GEN->Fmode);
  else
    fprintf(LOG,"%s: CDF(mode) unknown\n",gen->genid);
  fprintf(LOG,"%s: no (universal) squeeze\n",gen->genid);
  fprintf(LOG,"%s: no mirror principle\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: Rectangles:\n",gen->genid);
  if (GEN->ul > 0.)
    fprintf(LOG,"%s:    left upper point  = (%g,%g) \tarea = %g   (%5.2f%%)\n",
	    gen->genid,GEN->al/GEN->ul,GEN->ul,fabs(GEN->al),100.*fabs(GEN->al)/(-GEN->al+GEN->ar));
  else
    fprintf(LOG,"%s:    left upper point  = (0,0) \tarea = 0   (0.00%%)\n",gen->genid);
  fprintf(LOG,"%s:    right upper point = (%g,%g) \tarea = %g   (%5.2f%%)\n",
	  gen->genid,GEN->ar/GEN->ur,GEN->ur,GEN->ar,100.*GEN->ar/(-GEN->al+GEN->ar));
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dsrou_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PMF\n");
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %d   %s\n", DISTR.mode,
                      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   sum(PMF)  = %g\n", DISTR.sum);
  if (gen->set & DSROU_SET_CDFMODE)
    _unur_string_append(info,"   F(mode)   = %g\n", GEN->Fmode);
  else
    _unur_string_append(info,"   F(mode)   = [unknown]\n");
  _unur_string_append(info,"\n");
  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX ) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may provide the \"mode\"");
    _unur_string_append(info,"\n");
  }
  _unur_string_append(info,"method: DSROU (Discrete Simple Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   enveloping rectangle = (%g,%g) x (%g,%g)  [left]\n",
		      (GEN->ul > 0.)?GEN->al/GEN->ul:0., 0.,
		      0., (GEN->ul > 0.)?GEN->ul:0.);
  _unur_string_append(info,"                          (%g,%g) x (%g,%g)  [right]\n",
		      0.,GEN->ar/GEN->ur,  0., GEN->ur);
  _unur_string_append(info,"   area(hat) = %g + %g = %g\n",
		      fabs(GEN->al), GEN->ar, -GEN->al+GEN->ar);
  _unur_string_append(info,"   rejection constant = %g\n",
		      2. * (-GEN->al+GEN->ar) / DISTR.sum);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    if (gen->set & DSROU_SET_CDFMODE)
      _unur_string_append(info,"   cdfatmode = %g\n", GEN->Fmode);
    else
      _unur_string_append(info,"   cdfatmode = [not set]\n");
    if (gen->variant & DSROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & DSROU_SET_CDFMODE))
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"cdfatmode\" to reduce the rejection constant.");
    _unur_string_append(info,"\n");
  }
} 
#endif   
