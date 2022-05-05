/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "hist.h"
#include "hist_struct.h"
#define HIST_DEBUG_PRINTHIST   0x00000100u
#define GENTYPE "HIST"         
static struct unur_gen *_unur_hist_init( struct unur_par *par );
static struct unur_gen *_unur_hist_create( struct unur_par *par );
static struct unur_gen *_unur_hist_clone( const struct unur_gen *gen );
static void _unur_hist_free( struct unur_gen *gen);
static double _unur_hist_sample( struct unur_gen *gen );
static int _unur_hist_create_tables( struct unur_gen *gen );
static int _unur_hist_make_guidetable( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_hist_debug_init( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_hist_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cemp      
#define PAR       ((struct unur_hist_par*)par->datap) 
#define GEN       ((struct unur_hist_gen*)gen->datap) 
#define DISTR     gen->distr->data.cemp 
#define SAMPLE    gen->sample.cont           
#define _unur_hist_getSAMPLE(gen)   (_unur_hist_sample)
struct unur_par *
unur_hist_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CEMP,NULL);
  if (DISTR_IN.hist_prob == NULL || !(distr->set & UNUR_DISTR_SET_DOMAIN)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"histogram"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_hist_par) );
  COOKIE_SET(par,CK_HIST_PAR);
  par->distr    = distr;          
  par->method   = UNUR_METH_HIST; 
  par->variant  = 0u;             
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init     = _unur_hist_init;
  return par;
} 
struct unur_gen *
_unur_hist_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_HIST ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HIST_PAR,NULL);
  gen = _unur_hist_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if ( (_unur_hist_create_tables(gen) != UNUR_SUCCESS) ||
       (_unur_hist_make_guidetable(gen) != UNUR_SUCCESS) ) {
    _unur_hist_free(gen); return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_hist_debug_init(gen);
#endif
  return gen;
} 
struct unur_gen *
_unur_hist_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HIST_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_hist_gen) );
  COOKIE_SET(gen,CK_HIST_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_hist_getSAMPLE(gen);
  gen->destroy = _unur_hist_free;
  gen->clone = _unur_hist_clone;
  if (DISTR.hist_bins) {
    DISTR.hmin = DISTR.hist_bins[0];
    DISTR.hmax = DISTR.hist_bins[DISTR.n_hist];
  }
  GEN->n_hist = DISTR.n_hist;      
  GEN->prob   = DISTR.hist_prob;   
  GEN->hmin   = DISTR.hmin;        
  GEN->hmax   = DISTR.hmax;        
  GEN->hwidth = (DISTR.hmax - DISTR.hmin) / DISTR.n_hist;
  GEN->bins   = (DISTR.hist_bins) ? DISTR.hist_bins : NULL;
  GEN->sum = 0.;
  GEN->cumpv = NULL;
  GEN->guide_table = NULL;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_hist_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_hist_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_hist_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HIST_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->prob = clone->distr->data.cemp.hist_prob;   
  CLONE->bins = clone->distr->data.cemp.hist_bins;   
  CLONE->cumpv = _unur_xmalloc( GEN->n_hist * sizeof(double) );
  memcpy( CLONE->cumpv, GEN->cumpv, GEN->n_hist * sizeof(double) );
  CLONE->guide_table = _unur_xmalloc( GEN->n_hist * sizeof(int) );
  memcpy( CLONE->guide_table, GEN->guide_table, GEN->n_hist * sizeof(int) );
  return clone;
#undef CLONE
} 
void
_unur_hist_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_HIST ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HIST_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->guide_table) free(GEN->guide_table);
  if (GEN->cumpv)       free(GEN->cumpv);
  _unur_generic_free(gen);
} 
double
_unur_hist_sample( struct unur_gen *gen )
{ 
  double U;
  int J;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_HIST_GEN,UNUR_INFINITY);
  U = _unur_call_urng(gen->urng);
  J = GEN->guide_table[(int)(U * GEN->n_hist)];
  U *= GEN->sum;
  while (GEN->cumpv[J] < U) J++;
  U = (U - (J ? GEN->cumpv[J-1] : 0.)) / GEN->prob[J];
  if (GEN->bins) 
    return (U * GEN->bins[J+1] + (1.-U) * GEN->bins[J]);
  else
    return (GEN->hmin + (U+J)*GEN->hwidth);
} 
int
_unur_hist_create_tables( struct unur_gen *gen )
{ 
  GEN->cumpv = _unur_xrealloc( GEN->cumpv, GEN->n_hist * sizeof(double) );
  GEN->guide_table = _unur_xrealloc( GEN->guide_table, GEN->n_hist * sizeof(int) );
  return UNUR_SUCCESS;
} 
int
_unur_hist_make_guidetable( struct unur_gen *gen )
{ 
  double *pv;                   
  int n_pv;                     
  double pvh;                   
  double gstep;                 
  int i,j;
  pv = GEN->prob;
  n_pv = GEN->n_hist;
  for( i=0, pvh=0.; i<n_pv; i++ ) {
    GEN->cumpv[i] = ( pvh += pv[i] );
    if (pv[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      return UNUR_ERR_GEN_DATA;
    }
  }
  GEN->sum = GEN->cumpv[n_pv-1];
  gstep = GEN->sum / GEN->n_hist;
  pvh = 0.;
  for( j=0, i=0; j<GEN->n_hist;j++ ) {
    while (GEN->cumpv[i] < pvh) 
      i++;
    if (i >= n_pv) {
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
      break;
    }
    GEN->guide_table[j] = i;
    pvh += gstep;
  }
  for( ; j<GEN->n_hist; j++ )
    GEN->guide_table[j] = n_pv - 1;
  return UNUR_SUCCESS;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_hist_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HIST_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = HIST (HISTogram of empirical distribution)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cemp_debug( gen->distr, gen->genid, (gen->debug & HIST_DEBUG_PRINTHIST));
  fprintf(LOG,"%s: sampling routine = _unur_hist_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_hist_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = DATA  [histogram of size=%d]\n", DISTR.n_hist);
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: HIST (HISTogram of empirical distribution)\n");
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }
} 
#endif   
