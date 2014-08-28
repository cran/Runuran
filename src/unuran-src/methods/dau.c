/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dau.h"
#include "dau_struct.h"
#define DAU_DEBUG_REINIT       0x00000010u  
#define DAU_DEBUG_PRINTVECTOR  0x00000100u
#define DAU_DEBUG_TABLE        0x00000200u
#define DAU_SET_URNFACTOR       0x01u
#define GENTYPE "DAU"         
static struct unur_gen *_unur_dau_init( struct unur_par *par );
static int _unur_dau_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dau_create( struct unur_par *par );
static int _unur_dau_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_dau_clone( const struct unur_gen *gen );
static void _unur_dau_free( struct unur_gen *gen);
static int _unur_dau_sample( struct unur_gen *gen );
static int _unur_dau_create_tables( struct unur_gen *gen );
static int _unur_dau_make_urntable( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dau_debug_init( struct unur_gen *gen );
static void _unur_dau_debug_table( struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dau_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr      
#define PAR       ((struct unur_dau_par*)par->datap) 
#define GEN       ((struct unur_dau_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define SAMPLE    gen->sample.discr     
#define _unur_dau_getSAMPLE(gen)   (_unur_dau_sample)
struct unur_par *
unur_dau_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL(GENTYPE,distr,NULL);
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  if (DISTR_IN.pv == NULL) {
    if ( DISTR_IN.pmf
	 && ( (((unsigned)DISTR_IN.domain[1] - (unsigned)DISTR_IN.domain[0]) < UNUR_MAX_AUTO_PV)
	      || ( (distr->set & UNUR_DISTR_SET_PMFSUM) && DISTR_IN.domain[0] > INT_MIN ) ) ) {
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV. Try to compute it.");
    }
    else {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); return NULL;
    }
  }
  par = _unur_par_new( sizeof(struct unur_dau_par) );
  COOKIE_SET(par,CK_DAU_PAR);
  par->distr     = distr;            
  PAR->urn_factor = 1.;              
  par->method    = UNUR_METH_DAU;    
  par->variant   = 0u;               
  par->set       = 0u;                   
  par->urng      = unur_get_default_urng(); 
  par->urng_aux  = NULL;                    
  par->debug     = _unur_default_debugflag; 
  par->init = _unur_dau_init;
  return par;
} 
int
unur_dau_set_urnfactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DAU );
  if (factor < 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"relative urn size < 1.");
    return UNUR_ERR_PAR_SET;
  }
  PAR->urn_factor = factor;
  par->set |= DAU_SET_URNFACTOR;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dau_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_DAU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DAU_PAR,NULL);
  gen = _unur_dau_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if ( _unur_dau_check_par(gen) != UNUR_SUCCESS ) {
    _unur_dau_free(gen); return NULL;
  }
  if ( (_unur_dau_create_tables(gen) != UNUR_SUCCESS) ||
       (_unur_dau_make_urntable(gen) != UNUR_SUCCESS) ) {
    _unur_dau_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_dau_debug_init(gen);
#endif
  return gen;
} 
int
_unur_dau_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_dau_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  if ( ((rcode = _unur_dau_create_tables(gen)) != UNUR_SUCCESS) ||
       ((rcode = _unur_dau_make_urntable(gen)) != UNUR_SUCCESS) ) {
    return rcode;
  }
  SAMPLE = _unur_dau_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & DAU_DEBUG_REINIT) _unur_dau_debug_init(gen);
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dau_create( struct unur_par *par)
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DAU_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_dau_gen) );
  COOKIE_SET(gen,CK_DAU_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_dau_getSAMPLE(gen);
  gen->destroy = _unur_dau_free;
  gen->clone = _unur_dau_clone;
  gen->reinit = _unur_dau_reinit;
  GEN->urn_factor = PAR->urn_factor; 
  GEN->len = 0;             
  GEN->urn_size = 0;
  GEN->jx = NULL;
  GEN->qx = NULL;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dau_info;
#endif
  return gen;
} 
int
_unur_dau_check_par( struct unur_gen *gen )
{
  if (DISTR.pv == NULL) {
    if (unur_distr_discr_make_pv( gen->distr ) <= 0) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dau_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_dau_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DAU_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->jx = _unur_xmalloc( GEN->urn_size * sizeof(int) );
  memcpy( CLONE->jx, GEN->jx, GEN->urn_size * sizeof(int) );
  CLONE->qx = _unur_xmalloc( GEN->urn_size * sizeof(double) );
  memcpy( CLONE->qx, GEN->qx, GEN->urn_size * sizeof(double) );
  return clone;
#undef CLONE
} 
void
_unur_dau_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_DAU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DAU_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->jx) free(GEN->jx);
  if (GEN->qx) free(GEN->qx);
  _unur_generic_free(gen);
} 
int
_unur_dau_sample( struct unur_gen *gen )
{ 
  int iu;
  double u;
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DAU_GEN,INT_MAX);
  u = _unur_call_urng(gen->urng);
  u *= GEN->urn_size;
  iu = (int) u;
  if (iu >= GEN->len) return (GEN->jx[iu] + DISTR.domain[0]);
  u -= iu;   
  return (((u <= GEN->qx[iu]) ? iu : GEN->jx[iu] ) + DISTR.domain[0]);
} 
int
_unur_dau_create_tables( struct unur_gen *gen )
{ 
  GEN->len = DISTR.n_pv;
  GEN->urn_size = (int)(GEN->len * GEN->urn_factor);
  if (GEN->urn_size < GEN->len)
    GEN->urn_size = GEN->len;
  GEN->jx = _unur_xrealloc( GEN->jx, GEN->urn_size * sizeof(int) );
  GEN->qx = _unur_xrealloc( GEN->qx, GEN->urn_size * sizeof(double) );
  return UNUR_SUCCESS;
} 
int
_unur_dau_make_urntable( struct unur_gen *gen )
{ 
  int *begin, *poor, *rich;     
  int *npoor;                   
  double *pv;                   
  int n_pv;                     
  double sum, ratio;
  int i;                        
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;
  for( sum=0, i=0; i<n_pv; i++ ) {
    sum += pv[i];
    if (pv[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      return UNUR_ERR_GEN_DATA;
    }
  }
  begin = _unur_xmalloc( (GEN->urn_size+2) * sizeof(int) );
  poor = begin;                    
  rich = begin + GEN->urn_size + 1; 
  ratio = GEN->urn_size / sum;
  for( i=0; i<n_pv; i++ ) {
    GEN->qx[i] = pv[i] * ratio;  
    if (GEN->qx[i] >= 1.) {      
      *rich = i;                
      --rich;                   
      GEN->jx[i] = i;            
    }
    else {                      
      *poor = i;                
      ++poor;                   
    }
  }
  for( ; i<GEN->urn_size; i++ ) {
    GEN->qx[i] = 0.;
    *poor = i; 
    ++poor;
  }
  if (rich == begin + GEN->urn_size + 1 ) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    free (begin);
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  ++rich;
  while (poor != begin) {
    if (rich > begin + GEN->urn_size + 1) {
      break;
    }
    npoor = poor - 1;                       
    GEN->jx[*npoor] = *rich;                 
    GEN->qx[*rich] -= 1. - GEN->qx[*npoor];   
    if (GEN->qx[*rich] < 1.) {
      *npoor = *rich;      
      ++rich;              
    }
    else
      --poor;              
  }
  if (poor != begin) {
    sum = 0.;                   
    while (poor != begin) {
      npoor = poor - 1;         
      sum += 1. - GEN->qx[*npoor];
      GEN->jx[*npoor] = *npoor;  
      GEN->qx[*npoor] = 1.;      
      --poor;                   
    }
    if (fabs(sum) > UNUR_SQRT_DBL_EPSILON)
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"squared histogram");
  }
  free(begin);
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_dau_debug_init( struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DAU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variate\n",gen->genid);
  fprintf(LOG,"%s: method  = alias and alias-urn method\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_discr_debug( gen->distr,gen->genid,(gen->debug & DAU_DEBUG_PRINTVECTOR));
  fprintf(LOG,"%s: sampling routine = _unur_dau_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: length of probability vector = %d\n",gen->genid,GEN->len);
  fprintf(LOG,"%s: size of urn table = %d   (rel. = %g%%",
	  gen->genid,GEN->urn_size,100.*GEN->urn_factor);
  _unur_print_if_default(gen,DAU_SET_URNFACTOR);
  if (GEN->urn_size == GEN->len)
    fprintf(LOG,")   (--> alias method)\n");
  else
    fprintf(LOG,")   (--> alias-urn method)\n");
  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->debug & DAU_DEBUG_TABLE) {
    _unur_dau_debug_table(gen);
    fprintf(LOG,"%s:\n",gen->genid);
  }
} 
#define HIST_WIDTH   40  
void
_unur_dau_debug_table( struct unur_gen *gen )
{
  FILE *LOG;
  int i, j, m;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DAU_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: alias table:\n", gen->genid); 
  fprintf(LOG,"%s:\n", gen->genid);
  fprintf(LOG,"%s:         ratio donor/acceptor",gen->genid);
  for (i=0; i<HIST_WIDTH-17; i++)
    fprintf(LOG," ");
  fprintf(LOG,"jx:     qx:\n");
  for (i=0; i<GEN->urn_size; i++){
    m = HIST_WIDTH * GEN->qx[i] + 0.5;
    fprintf(LOG,"%s:[%4d]: ", gen->genid,i); 
    for (j=0; j<HIST_WIDTH; j++)
      if (j<m)
	fprintf(LOG, "*"); 
      else                
	fprintf(LOG,"-");
    fprintf(LOG," %5d  ", GEN->jx[i]);           
    fprintf(LOG,"  %6.3f%%\n", GEN->qx[i]*100);  
  }
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dau_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PV  [length=%d%s]\n",
		      DISTR.domain[1]-DISTR.domain[0]+1,
		      (DISTR.pmf==NULL) ? "" : ", created from PMF");
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: DAU (Alias-Urn)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#look-ups] = %g\n", 1+1./GEN->urn_factor);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   urnfactor = %g  %s\n", GEN->urn_factor,
			(gen->set & DAU_SET_URNFACTOR) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }
} 
#endif   
