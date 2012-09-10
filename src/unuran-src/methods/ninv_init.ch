/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_gen *
_unur_ninv_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);
  if (par->variant == NINV_VARFLAG_NEWTON && ! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   
  }
  gen = _unur_ninv_create(par);
  _unur_par_free(par);
  if (!gen) { return NULL; }
  if (_unur_ninv_check_par(gen) != UNUR_SUCCESS) {
    _unur_ninv_free(gen); return NULL;
  }
  if (GEN->table_on) {
    if (_unur_ninv_create_table(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }
  else {
    if (_unur_ninv_compute_start(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_ninv_debug_init(gen);
#endif
  return gen;
} 
int
_unur_ninv_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_ninv_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  if (DISTR.upd_area != NULL)
    if ((DISTR.upd_area)(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cannot compute normalization constant");
      return UNUR_ERR_GEN_DATA;
    }
  if (GEN->table != NULL)
    rcode = _unur_ninv_create_table(gen);
  else 
    rcode = unur_ninv_chg_start( gen, 0., 0. );
  SAMPLE = _unur_ninv_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_REINIT) {
    _unur_distr_cont_debug( gen->distr, gen->genid );
    if (rcode==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
  }
#endif
  return UNUR_SUCCESS;
} 
static struct unur_gen *
_unur_ninv_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NINV_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_ninv_gen) );
  COOKIE_SET(gen,CK_NINV_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_ninv_getSAMPLE(gen);
  gen->destroy = _unur_ninv_free;
  gen->clone = _unur_ninv_clone;
  gen->reinit = _unur_ninv_reinit;
  GEN->max_iter = PAR->max_iter;      
  GEN->x_resolution = PAR->x_resolution; 
  GEN->u_resolution = PAR->u_resolution; 
  GEN->table_on = PAR->table_on;      
  GEN->table_size = PAR->table_size;  
  GEN->s[0] = PAR->s[0];              
  GEN->s[1] = PAR->s[1];
  GEN->table = NULL;
  GEN->f_table = NULL;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_ninv_info;
#endif
  return gen;
} 
int
_unur_ninv_check_par( struct unur_gen *gen )
{
  if ( GEN->x_resolution < 0. && GEN->u_resolution < 0. ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"both x-resolution and u-resolution negativ. using defaults.");
    GEN->x_resolution = 1.e-8;
  }
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];
  GEN->CDFmin = GEN->Umin = (DISTR.trunc[0] > -UNUR_INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN->CDFmax = GEN->Umax = (DISTR.trunc[1] < UNUR_INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;
  if (_unur_FP_greater(GEN->CDFmin, GEN->CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    return UNUR_ERR_GEN_DATA;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_ninv_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_ninv_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NINV_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  if (GEN->table) {
    CLONE->table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->table, GEN->table, GEN->table_size * sizeof(double) );
    CLONE->f_table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->f_table, GEN->f_table, GEN->table_size * sizeof(double) );
  }
  return clone;
#undef CLONE
} 
void
_unur_ninv_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->table)   free(GEN->table);
  if (GEN->f_table) free(GEN->f_table);
  _unur_generic_free(gen);
} 
int
_unur_ninv_create_table( struct unur_gen *gen )
{
  int i;
  double x;
  int table_size = GEN->table_size;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);
  GEN->table    = _unur_xrealloc( GEN->table,   table_size * sizeof(double));
  GEN->f_table  = _unur_xrealloc( GEN->f_table, table_size * sizeof(double));
  GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
  GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
  GEN->CDFs[0]  = CDF(GEN->s[0]);
  GEN->CDFs[1]  = CDF(GEN->s[1]);
  GEN->table_on = FALSE;
  GEN->table[0]              = DISTR.domain[0];
  GEN->f_table[0]            = GEN->CDFmin;    
  GEN->table[table_size-1]   = DISTR.domain[1];
  GEN->f_table[table_size-1] = GEN->CDFmax;    
  for (i=1; i<table_size/2; i++){
    x = GEN->CDFmin + i * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[i]   = _unur_ninv_regula(gen,x);
    GEN->f_table[i] = CDF(GEN->table[i]);
    x = GEN->CDFmin + (table_size-i-1) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size-1-i] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size-1-i] = CDF(GEN->table[table_size-1-i]);
    if (GEN->table[i] > -UNUR_INFINITY) {
      GEN->s[0] = GEN->table[i];
      GEN->CDFs[0] = GEN->f_table[i];
    }
    if (GEN->table[table_size-1-i] < UNUR_INFINITY) {
      GEN->s[1] = GEN->table[table_size-1-i];
      GEN->CDFs[1] = GEN->f_table[table_size-1-i];
    }
  }  
  if (table_size & 1) { 
    x = GEN->CDFmin + (table_size/2) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size/2] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size/2] = CDF(GEN->table[table_size/2]);
  }  
  GEN->table_on = TRUE;
  return UNUR_SUCCESS;
}  
int
_unur_ninv_compute_start( struct unur_gen *gen )
{
  double u;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);
  if( GEN->table_on )
    return UNUR_SUCCESS;
  if( !_unur_FP_same(GEN->s[0],GEN->s[1]) ) {
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    return UNUR_SUCCESS;
  }
  switch (gen->variant) {
  case NINV_VARFLAG_BISECT:
  case NINV_VARFLAG_REGULA:
    GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);    
    u = GEN->CDFmin + 0.5*(1.-INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[0] = _unur_ninv_regula(gen,u);
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    u = GEN->CDFmin + 0.5*(1.+INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[1] = _unur_ninv_regula(gen,u);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    break;    
  case NINV_VARFLAG_NEWTON:
    GEN->s[0] = _unur_max( DISTR.domain[0], -9.987655 );
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]); 
    GEN->CDFs[1] = CDF(GEN->s[1]);
    u = 0.5 * (GEN->CDFmin + GEN->CDFmax);
    GEN->s[0] = _unur_ninv_regula(gen, u); 
    GEN->CDFs[0] = CDF(GEN->s[0]);
    break;    
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }  
  return UNUR_SUCCESS;
}  
