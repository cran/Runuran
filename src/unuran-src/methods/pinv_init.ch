/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_gen *
_unur_pinv_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  double lfc;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_PINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_PINV_PAR,NULL);
  gen = _unur_pinv_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  if (_unur_pinv_check_par(gen) != UNUR_SUCCESS) {
    _unur_pinv_free(gen); return NULL;
  }
  if (DISTR.logpdf != NULL && (gen->variant & PINV_VARIANT_PDF) ) {
    lfc = UNUR_INFINITY;
    if ( (gen->distr->set & UNUR_DISTR_SET_MODE) &&
	 !_unur_FP_less(DISTR.mode,DISTR.domain[0]) &&
	 !_unur_FP_greater(DISTR.mode,DISTR.domain[1]) ) {
      lfc = (DISTR.logpdf)(DISTR.mode,gen->distr);
    }
    if (!_unur_isfinite(lfc))
      lfc = (DISTR.logpdf)(DISTR.center,gen->distr);
    if (lfc < -3.)
      GEN->logPDFconstant = lfc;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_pinv_debug_init_start(gen);
#endif
  if (_unur_pinv_preprocessing(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_pinv_debug_init(gen,FALSE);
#endif
    _unur_pinv_free(gen); return NULL;
  }
  if (_unur_pinv_create_table(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_pinv_debug_init(gen,FALSE);
#endif
    _unur_pinv_free(gen); return NULL;
  }
  _unur_lobatto_free(&(GEN->aCDF));
  _unur_pinv_make_guide_table(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_pinv_debug_init(gen,TRUE);
#endif
  return gen;
} 
struct unur_gen *
_unur_pinv_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_PINV_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_pinv_gen) );
  COOKIE_SET(gen,CK_PINV_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_pinv_getSAMPLE(gen);
  gen->destroy = _unur_pinv_free;
  gen->clone = _unur_pinv_clone;
  GEN->order = PAR->order;            
  GEN->smooth = PAR->smooth;          
  GEN->u_resolution = PAR->u_resolution; 
  GEN->bleft_par  = PAR->bleft;          
  GEN->bright_par = PAR->bright;
  GEN->sleft  = PAR->sleft;              
  GEN->sright = PAR->sright;
  GEN->max_ivs = PAR->max_ivs;           
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  GEN->dleft = -INFINITY;
  GEN->dright = INFINITY;
  GEN->Umax = 1.;
  GEN->iv = NULL;
  GEN->n_ivs = -1;        
  GEN->guide_size = 0; 
  GEN->guide = NULL;
  GEN->area = DISTR.area; 
  GEN->logPDFconstant = 0.;   
  GEN->aCDF = NULL;           
  GEN->iv = _unur_xmalloc(GEN->max_ivs * sizeof(struct unur_pinv_interval) );
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_pinv_info;
#endif
  return gen;
} 
int
_unur_pinv_check_par( struct unur_gen *gen )
{
  switch (GEN->smooth) {
  case 2:
    if (GEN->order < 5) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"order must be >= 5 when smoothness equals 2");
      GEN->order = 5;
      gen->set |= PINV_SET_ORDER_COR;
    }
    if (GEN->order % 3 != 2) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"order must be 2 mod 3 when smoothness equals 2");
      GEN->order = 2 + 3 * (GEN->order / 3);
      gen->set |= PINV_SET_ORDER_COR;
    }
    if (DISTR.pdf == NULL || DISTR.dpdf == NULL) {
      _unur_warning(gen->genid,UNUR_ERR_DISTR_REQUIRED,"PDF or dPDF --> try smoothness=1 instead");
      GEN->smooth = 1;
      gen->set |= PINV_SET_SMOOTH_COR;
    }
    else {
      break;
    }
  case 1:
    if (GEN->order % 2 != 1) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"order must be odd when smoothness equals 1");
      GEN->order += 1;
      gen->set |= PINV_SET_ORDER_COR;
    }
    if (DISTR.pdf == NULL) {
      _unur_warning(gen->genid,UNUR_ERR_DISTR_REQUIRED,"PDF --> use smoothness=0 instead");
      GEN->smooth = 0;
      gen->set |= PINV_SET_SMOOTH_COR;
    }
    else { 
      break;
    }
  case 0:
    break;
  default:
    _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"smoothness must be 0, 1, or 2");
    GEN->smooth = 0; 
  }
  GEN->bleft = _unur_max(GEN->bleft_par,DISTR.domain[0]);
  GEN->bright = _unur_min(GEN->bright_par,DISTR.domain[1]);
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];
  GEN->dleft =  DISTR.domain[0];
  GEN->dright =  DISTR.domain[1];
  DISTR.center = unur_distr_cont_get_center(gen->distr);
  if (DISTR.center < GEN->dleft || DISTR.center > GEN->dright) {
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,
		"center moved into domain of distribution");
    DISTR.center = _unur_max(DISTR.center,GEN->dleft);
    DISTR.center = _unur_min(DISTR.center,GEN->dright);
  }
  if (gen->variant & PINV_VARIANT_PDF) {
    if (PDF(DISTR.center)<=0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "PDF(center) <= 0.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_pinv_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_pinv_gen*)clone->datap)
  struct unur_gen *clone;
  int i;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->aCDF = NULL;
  CLONE->iv =  _unur_xmalloc((GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  memcpy( CLONE->iv, GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  for(i=0; i<=GEN->n_ivs; i++) {
    CLONE->iv[i].ui = _unur_xmalloc( GEN->order * sizeof(double) );
    CLONE->iv[i].zi = _unur_xmalloc( GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].ui, GEN->iv[i].ui, GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].zi, GEN->iv[i].zi, GEN->order * sizeof(double) );
  }
  CLONE->guide = _unur_xmalloc( GEN->guide_size * sizeof(int) );
  memcpy( CLONE->guide, GEN->guide, GEN->guide_size * sizeof(int) );
  return clone;
#undef CLONE
} 
void
_unur_pinv_free( struct unur_gen *gen )
{ 
  int i;
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->guide) free (GEN->guide);
  _unur_lobatto_free(&(GEN->aCDF));
  if (GEN->iv) {
    for(i=0; i<=GEN->n_ivs; i++){
      free(GEN->iv[i].ui);
      free(GEN->iv[i].zi);
    }
    free (GEN->iv);
  }
  _unur_generic_free(gen);
} 
int
_unur_pinv_make_guide_table (struct unur_gen *gen)
{
  int i,j, imax;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);
  GEN->guide_size = (int) (GEN->n_ivs * PINV_GUIDE_FACTOR);
  if (GEN->guide_size <= 0) GEN->guide_size = 1;
  GEN->guide = _unur_xrealloc( GEN->guide, GEN->guide_size * sizeof(int) );
  imax = GEN->n_ivs;
  i = 0;
  GEN->guide[0] = 0;
  for( j=1; j<GEN->guide_size ;j++ ) {
    while(GEN->iv[i+1].cdfi/GEN->Umax < j/(double)GEN->guide_size && i < imax)
      i++;
    if (i >= imax) break;
    GEN->guide[j]=i;
  }
  i = _unur_min(i,imax);
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = i;
  return UNUR_SUCCESS;
} 
double
_unur_pinv_eval_PDF (double x, struct unur_gen *gen)
{
  struct unur_distr *distr = gen->distr;
  double fx, dx;
  int i;
  for (i=1; i<=2; i++) {
    if (DISTR.logpdf != NULL)
      fx = exp((DISTR.logpdf)(x,distr) - GEN->logPDFconstant);
    else
      fx = (DISTR.pdf)(x,distr);
    if (fx >= INFINITY) {
      dx = 2.*fabs(x)*DBL_EPSILON;
      dx = _unur_max(dx,2.*DBL_MIN);
      x += ((x - GEN->bleft) < (GEN->bright - x)) ? dx : -dx; 
    }
    else 
      break;
  }
  return fx;
} 
