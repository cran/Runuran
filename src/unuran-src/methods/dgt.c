/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dgt.h"
#include "dgt_struct.h"
#define DGT_VARFLAG_DIV     0x01u     
#define DGT_VARFLAG_ADD     0x02u     
#define DGT_VAR_THRESHOLD   1000      
#define DGT_DEBUG_REINIT       0x00000010u  
#define DGT_DEBUG_PRINTVECTOR  0x00000100u
#define DGT_DEBUG_TABLE        0x00000200u
#define DGT_SET_GUIDEFACTOR    0x010u
#define DGT_SET_VARIANT        0x020u
#define GENTYPE "DGT"         
static struct unur_gen *_unur_dgt_init( struct unur_par *par );
static int _unur_dgt_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_dgt_create( struct unur_par *par );
static int _unur_dgt_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_dgt_clone( const struct unur_gen *gen );
static void _unur_dgt_free( struct unur_gen *gen);
static int _unur_dgt_sample( struct unur_gen *gen );
static int _unur_dgt_create_tables( struct unur_gen *gen );
static int _unur_dgt_make_guidetable( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_dgt_debug_init( struct unur_gen *gen );
static void _unur_dgt_debug_table( struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_dgt_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.discr      
#define PAR       ((struct unur_dgt_par*)par->datap) 
#define GEN       ((struct unur_dgt_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define SAMPLE    gen->sample.discr     
#define _unur_dgt_getSAMPLE(gen)  (_unur_dgt_sample)
struct unur_par *
unur_dgt_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
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
  par = _unur_par_new( sizeof(struct unur_dgt_par) );
  COOKIE_SET(par,CK_DGT_PAR);
  par->distr       = distr;          
  PAR->guide_factor = 1.;            
  par->method      = UNUR_METH_DGT;  
  par->variant     = 0u;             
  par->set         = 0u;                 
  par->urng        = unur_get_default_urng(); 
  par->urng_aux    = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_dgt_init;
  return par;
} 
int
unur_dgt_set_variant( struct unur_par *par, unsigned variant )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DGT );
  if (variant != DGT_VARFLAG_ADD && variant != DGT_VARFLAG_DIV) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
    return UNUR_ERR_PAR_VARIANT;
  }
  par->set |= DGT_SET_VARIANT;
  par->variant = variant;
  return UNUR_SUCCESS;
} 
int
unur_dgt_set_guidefactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DGT );
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"relative table size < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->guide_factor = factor;
  par->set |= DGT_SET_GUIDEFACTOR;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dgt_init( struct unur_par *par )
{ 
  struct unur_gen *gen;         
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_DGT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DGT_PAR,NULL);
  gen = _unur_dgt_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if ( _unur_dgt_check_par(gen) != UNUR_SUCCESS ) {
    _unur_dgt_free(gen); return NULL;
  }
  if ( (_unur_dgt_create_tables(gen) != UNUR_SUCCESS) ||
       (_unur_dgt_make_guidetable(gen) != UNUR_SUCCESS) ) {
    _unur_dgt_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_dgt_debug_init(gen);
#endif
  return gen;
} 
int
_unur_dgt_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_dgt_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  if ( ((rcode = _unur_dgt_create_tables(gen)) != UNUR_SUCCESS) ||
       ((rcode = _unur_dgt_make_guidetable(gen)) != UNUR_SUCCESS) ) {
    return rcode;
  }
  SAMPLE = _unur_dgt_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & DGT_DEBUG_REINIT) _unur_dgt_debug_init(gen);
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dgt_create( struct unur_par *par )
{
  struct unur_gen *gen;       
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DGT_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_dgt_gen) );
  COOKIE_SET(gen,CK_DGT_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_dgt_getSAMPLE(gen);
  gen->destroy = _unur_dgt_free;
  gen->clone = _unur_dgt_clone;
  gen->reinit = _unur_dgt_reinit;
  GEN->guide_factor = PAR->guide_factor;
  GEN->cumpv = NULL;
  GEN->guide_table = NULL;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_dgt_info;
#endif
  return gen;
} 
int
_unur_dgt_check_par( struct unur_gen *gen )
{
  if (DISTR.pv == NULL) {
    if (unur_distr_discr_make_pv( gen->distr ) <= 0) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }
  if (gen->variant == 0)   
    gen->variant = (DISTR.n_pv > DGT_VAR_THRESHOLD) 
      ? DGT_VARFLAG_DIV : DGT_VARFLAG_ADD;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_dgt_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_dgt_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DGT_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->cumpv = _unur_xmalloc( DISTR.n_pv * sizeof(double) );
  memcpy( CLONE->cumpv, GEN->cumpv, DISTR.n_pv * sizeof(double) );
  CLONE->guide_table = _unur_xmalloc( GEN->guide_size * sizeof(int) );
  memcpy( CLONE->guide_table, GEN->guide_table, GEN->guide_size * sizeof(int) );
  return clone;
#undef CLONE
} 
void
_unur_dgt_free( struct unur_gen *gen )
{ 
  if (!gen) 
    return;
  if ( gen->method != UNUR_METH_DGT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DGT_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->guide_table) free(GEN->guide_table);
  if (GEN->cumpv)       free(GEN->cumpv);
  _unur_generic_free(gen);
} 
int
_unur_dgt_sample( struct unur_gen *gen )
{ 
  int j;
  double u;
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DGT_GEN,INT_MAX);
  u = _unur_call_urng(gen->urng);
  j = GEN->guide_table[(int)(u * GEN->guide_size)];
  u *= GEN->sum;
  while (GEN->cumpv[j] < u) j++;
  return (j + DISTR.domain[0]);
} 
int
unur_dgt_eval_invcdf_recycle( const struct unur_gen *gen, double u, double *recycle )
{
  int j;
  if (recycle) *recycle = 0.;
  _unur_check_NULL( GENTYPE, gen, INT_MAX );
  if ( gen->method != UNUR_METH_DGT ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INT_MAX;
  }
  COOKIE_CHECK(gen,CK_DGT_GEN,INT_MAX);
  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return INT_MAX;  
  }
  j = GEN->guide_table[(int)(u * GEN->guide_size)];
  u *= GEN->sum;
  while (GEN->cumpv[j] < u) j++;
  if (recycle) {
    *recycle = 1. - (GEN->cumpv[j] - u) / DISTR.pv[j];
  }
  j+=DISTR.domain[0];
  if (j<DISTR.domain[0]) j = DISTR.domain[0];
  if (j>DISTR.domain[1]) j = DISTR.domain[1];
  return j;
} 
int
unur_dgt_eval_invcdf( const struct unur_gen *gen, double u )
{
  return unur_dgt_eval_invcdf_recycle(gen,u,NULL);
} 
int
_unur_dgt_create_tables( struct unur_gen *gen )
{ 
  GEN->guide_size = (int)( DISTR.n_pv * GEN->guide_factor);
  if (GEN->guide_size <= 0)
    GEN->guide_size = 1;
  GEN->cumpv = _unur_xrealloc( GEN->cumpv, DISTR.n_pv * sizeof(double) );
  GEN->guide_table = _unur_xrealloc( GEN->guide_table, GEN->guide_size * sizeof(int) );
  return UNUR_SUCCESS;
} 
int
_unur_dgt_make_guidetable( struct unur_gen *gen )
{ 
  double *pv;                   
  int n_pv;                     
  double pvh;                   
  double gstep;                 
  int i,j;
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;
  for( i=0, pvh=0.; i<n_pv; i++ ) {
    GEN->cumpv[i] = ( pvh += pv[i] );
    if (pv[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      return UNUR_ERR_GEN_DATA;
    }
  }
  GEN->sum = GEN->cumpv[n_pv-1];
  if (gen->variant == DGT_VARFLAG_DIV) {
    GEN->guide_table[0] = 0;
    for( j=1, i=0; j<GEN->guide_size ;j++ ) {
      while( GEN->cumpv[i]/GEN->sum < ((double)j)/GEN->guide_size ) 
	i++;
      if (i >= n_pv) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN->guide_table[j]=i;
    }
  }
  else { 
    gstep = GEN->sum / GEN->guide_size;
    pvh = 0.;
    for( j=0, i=0; j<GEN->guide_size ;j++ ) {
      while (GEN->cumpv[i] < pvh) 
	i++;
      if (i >= n_pv) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN->guide_table[j] = i;
      pvh += gstep;
    }
  }
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide_table[j] = n_pv - 1;
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_dgt_debug_init( struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DGT_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = indexed search (guide table)\n",gen->genid);
  fprintf(LOG,"%s: variant = %u ",gen->genid,gen->variant);
  _unur_print_if_default(gen,DGT_SET_VARIANT);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_distr_discr_debug( gen->distr,gen->genid,(gen->debug & DGT_DEBUG_PRINTVECTOR));
  fprintf(LOG,"%s: sampling routine = _unur_dgt_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: length of probability vector = %d\n",gen->genid,DISTR.n_pv);
  fprintf(LOG,"%s: length of guide table = %d   (rel. = %g%%",
	  gen->genid,GEN->guide_size,100.*GEN->guide_factor);
  _unur_print_if_default(gen,DGT_SET_GUIDEFACTOR);
  if (GEN->guide_size == 1) 
    fprintf(LOG,") \t (-->sequential search");
  fprintf(LOG,")\n%s:\n",gen->genid);
  fprintf(LOG,"%s: sum over PMF (as computed) = %#-20.16g\n",gen->genid,GEN->sum);
  if (gen->debug & DGT_DEBUG_TABLE)
    _unur_dgt_debug_table(gen);
} 
void
_unur_dgt_debug_table( struct unur_gen *gen )
{   
  FILE *LOG;
  int i,j,m;
  int n_asts;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DGT_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: guide table:\n", gen->genid); 
  fprintf(LOG,"%s:\n", gen->genid);
  n_asts = 0;
  for (i=0; i<GEN->guide_size; i++){
    fprintf(LOG,"%s: [%5d] -> %5d ", gen->genid, i, GEN->guide_table[i]);
    if (i == GEN->guide_size-1)
      j = GEN->guide_size - GEN->guide_table[i];
    else
      j = GEN->guide_table[i+1] - GEN->guide_table[i] + 1;
    for (m=0; m<j && m<10; m++ ) {
      fprintf(LOG," *");
      ++n_asts;
    }
    if (m<j){
      n_asts += j-m;
      fprintf(LOG," ... %d", j);
    }
    fprintf(LOG,"\n");
  }
  fprintf(LOG,"%s:\n", gen->genid);
  fprintf(LOG,"%s: expected number of comparisons = %g\n",gen->genid,
          ((double)n_asts)/GEN->guide_size);
  fprintf(LOG,"%s:\n", gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_dgt_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"method: DGT (Guide Table)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#look-ups] = %g\n", 1+1./GEN->guide_factor);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   guidefactor = %g  %s\n", GEN->guide_factor,
			(gen->set & DGT_SET_GUIDEFACTOR) ? "" : "[default]");
    if (gen->set & DGT_SET_VARIANT)
      _unur_string_append(info,"   variant = %d\n", gen->variant);
    _unur_string_append(info,"\n");
  }
} 
#endif   
