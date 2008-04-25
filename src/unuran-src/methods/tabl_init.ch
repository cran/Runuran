/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_gen *
_unur_tabl_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_TABL ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  gen = _unur_tabl_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_tabl_debug_init_start(par,gen);
#endif
  do {
    if (PAR->n_slopes > 0) {
      if (_unur_tabl_get_intervals_from_slopes(par,gen)!=UNUR_SUCCESS) {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make hat function");
	_unur_par_free(par); _unur_tabl_free(gen); return NULL;
      }
      break;
    }
    if (PAR->cpoints == NULL) {
      if (! ( (par->distr->set & UNUR_DISTR_SET_MODE) &&
	      unur_tabl_set_cpoints(par,1,&(DISTR.mode)) == UNUR_SUCCESS) ) {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot compute slopes");
	_unur_par_free(par); _unur_tabl_free(gen); return NULL;
      }
    }
    if (PAR->cpoints != NULL) {
      if (_unur_tabl_get_intervals_from_cpoints(par,gen)!=UNUR_SUCCESS) {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make hat function");
	_unur_par_free(par); _unur_tabl_free(gen); return NULL;
      }
      break;
    }
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  } while(0);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TABL_DEBUG_IV_START) _unur_tabl_debug_intervals(gen,"starting intervals:",FALSE);
#endif
  PAR->n_slopes = GEN->n_ivs;
  if (_unur_tabl_compute_intervals(par,gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot split intervals");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  }
  if (_unur_tabl_make_guide_table(gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  }
  gen->status = UNUR_SUCCESS;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_tabl_debug_init_finished(gen);
#endif
  _unur_par_free(par);
  return gen;
} 
struct unur_gen *
_unur_tabl_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_tabl_gen) );
  COOKIE_SET(gen,CK_TABL_GEN);
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA))
    if (unur_distr_cont_upd_pdfarea(gen->distr)!=UNUR_SUCCESS)
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF, use default instead");
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_tabl_getSAMPLE(gen);
  gen->destroy = _unur_tabl_free;
  gen->clone = _unur_tabl_clone;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  if (par->distr->set & UNUR_DISTR_SET_DOMAIN) {
    PAR->bleft  = _unur_max(PAR->bleft, DISTR.BD_LEFT);
    PAR->bright = _unur_min(PAR->bright,DISTR.BD_RIGHT);
  }
  GEN->bleft       = PAR->bleft;         
  GEN->bright      = PAR->bright;        
  GEN->Umin = 0.;
  GEN->Umax = 1.;
  GEN->guide_factor = PAR->guide_factor; 
  GEN->max_ivs   = PAR->max_ivs;         
#ifdef UNUR_ENABLE_INFO
  GEN->max_ivs_info = PAR->max_ivs;      
#endif
  GEN->max_ratio = PAR->max_ratio;       
  GEN->darsfactor = PAR->darsfactor;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_tabl_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_tabl_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_tabl_gen*)clone->datap)
  struct unur_gen *clone;
  struct unur_tabl_interval *iv,*next, *clone_iv, *clone_prev;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = next) {
    clone_iv = _unur_xmalloc( sizeof(struct unur_tabl_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tabl_interval) );
    if (clone_prev == NULL) {
      CLONE->iv = clone_iv;
    }
    else {
      clone_prev->next = clone_iv;
    }
    next = iv->next;
    clone_prev = clone_iv;
  }
  if (clone_iv) clone_iv->next = NULL;
  CLONE->guide = NULL;
  if (_unur_tabl_make_guide_table(clone) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
  }
  return clone;
#undef CLONE
} 
void
_unur_tabl_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_TABL ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  SAMPLE = NULL;   
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_tabl_debug_free(gen);
#endif
  {
    struct unur_tabl_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }
  if (GEN->guide)  free(GEN->guide);
  _unur_generic_free(gen);
} 
int
_unur_tabl_get_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen )
{
  struct unur_tabl_interval *iv;
  double xmax, xmin;
  double sl, sr;
  int i;
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  GEN->n_ivs = 0;
  iv = NULL;
  GEN->bleft = INFINITY;
  GEN->bright = -INFINITY;
  for ( i=0; i < 2*PAR->n_slopes; i+=2 ) {
    xmax = PAR->slopes[i];      
    xmin = PAR->slopes[i+1];    
    if (xmax > xmin) { 
      sl = xmin; sr = xmax; }
    else {  
      sl = xmax; sr = xmin; }
    if (_unur_FP_greater(DISTR.BD_LEFT,sr)) continue;   
    if (_unur_FP_less(DISTR.BD_RIGHT,sl)) continue;     
    if (_unur_FP_greater(DISTR.BD_LEFT,sl)) {
      if (xmax > xmin) xmin = DISTR.BD_LEFT; else xmax = DISTR.BD_LEFT; }
    if (_unur_FP_less(DISTR.BD_RIGHT,sr)) {
      if (xmax < xmin) xmin = DISTR.BD_RIGHT; else xmax = DISTR.BD_RIGHT; }
    if (GEN->iv==NULL)  
      iv = GEN->iv = _unur_xmalloc(sizeof(struct unur_tabl_interval));   
    else       
      iv = iv->next = _unur_xmalloc(sizeof(struct unur_tabl_interval));
    ++(GEN->n_ivs);
    COOKIE_SET(iv,CK_TABL_IV);
    iv->xmax = xmax;
    iv->fmax = PDF(iv->xmax);
    iv->xmin = xmin;
    iv->fmin = PDF(iv->xmin);
    if (! (_unur_isfinite(iv->fmax) && _unur_isfinite(iv->fmin)
	   && iv->fmax >= 0. && iv->fmin >= 0.) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      iv->next = NULL;  
      return UNUR_ERR_GEN_DATA;
    }
    if (_unur_FP_less(iv->fmax,iv->fmin)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"PDF discontinuous or slope not monotone");
    }
    iv->fmin = 0.;
    iv->Ahat = fabs(xmax - xmin) * iv->fmax;
    iv->Asqueeze = fabs(xmax - xmin) * iv->fmin;
    iv->Acum = 0.;
    if (xmax > xmin) {
      GEN->bleft = _unur_min(GEN->bleft,xmin);
      GEN->bright = _unur_max(GEN->bright,xmax);
    }
    else {
      GEN->bleft = _unur_min(GEN->bleft,xmax);
      GEN->bright = _unur_max(GEN->bright,xmin);
    }
  }
  if (GEN->iv==NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"invalid slopes");
    return UNUR_ERR_GEN_DATA;
  }
  iv->next = NULL;
  gen->distr->set &= ~UNUR_DISTR_SET_PDFAREA;
  unur_distr_cont_upd_pdfarea( gen->distr );
  return UNUR_SUCCESS;
} 
int
_unur_tabl_get_intervals_from_cpoints( struct unur_par *par, struct unur_gen *gen )
{
  struct unur_tabl_interval *iv;
  double sl, sr; 
  double fl, fr; 
  double cp;     
  int i;
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  GEN->n_ivs = 0;
  iv = NULL;
  sr = GEN->bleft;
  fr = PDF(sr);
  for (i=0; i<=PAR->n_cpoints; i++) {
    if (i < PAR->n_cpoints) {
      cp = PAR->cpoints[i];
      if (! _unur_FP_less(GEN->bleft,cp))
	continue;   
      if (! _unur_FP_greater(GEN->bright,cp)) {
	i = (PAR->n_cpoints)-1;  
	continue; 
      }
    }
    else { 
      cp = GEN->bright;
    }
    sl = sr; fl = fr;
    sr = cp; fr = PDF(sr);
    if (GEN->iv==NULL)  
      iv = GEN->iv = _unur_xmalloc(sizeof(struct unur_tabl_interval));
    else       
      iv = iv->next = _unur_xmalloc(sizeof(struct unur_tabl_interval));
    ++(GEN->n_ivs);
    COOKIE_SET(iv,CK_TABL_IV);
    if (! (_unur_isfinite(fr) && _unur_isfinite(fl)
	   && fr >= 0. && fl >= 0.) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      iv->next = NULL;  
      return UNUR_ERR_GEN_DATA;
    }
    if (fr > fl) {
      iv->xmax = sr; iv->fmax = fr;
      iv->xmin = sl; iv->fmin = fl;
    }
    else {
      iv->xmax = sl; iv->fmax = fl;
      iv->xmin = sr; iv->fmin = fr;
    }
    iv->Ahat = fabs(sr - sl) * iv->fmax;
    iv->Asqueeze = fabs(sr - sl) * iv->fmin;
    iv->Acum = 0.;
  }
  if (GEN->iv==NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"invalid slopes");
    return UNUR_ERR_GEN_DATA;
  }
  iv->next = NULL;
  DISTR.trunc[0] = DISTR.BD_LEFT = GEN->bleft;
  DISTR.trunc[1] = DISTR.BD_RIGHT = GEN->bright;
  gen->distr->set &= ~UNUR_DISTR_SET_PDFAREA;
  unur_distr_cont_upd_pdfarea( gen->distr );
  return UNUR_SUCCESS;
} 
int
_unur_tabl_compute_intervals( struct unur_par *par, struct unur_gen *gen )
{
  struct unur_tabl_interval *iv;
  int i,k;
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  if (par->variant & TABL_VARFLAG_USEEAR) {
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
      iv = _unur_tabl_run_equalarearule( par, gen, iv );
      if (iv == NULL) return UNUR_ERR_GEN_DATA;
    }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TABL_DEBUG_IV_START) _unur_tabl_debug_intervals(gen,"equal area rule applied:",FALSE);
#endif
  }
  if (par->variant & TABL_VARFLAG_USEDARS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TABL_DEBUG_DARS) _unur_tabl_debug_dars_start(par,gen);
#endif
    for (i=0; i<TABL_N_RETRY_DARS; i++) {
      if (_unur_tabl_run_dars(gen)!=UNUR_SUCCESS)
	return UNUR_ERR_GEN_DATA;
      if (GEN->n_ivs >= GEN->max_ivs)
	break;
      if (_unur_tabl_make_guide_table(gen) != UNUR_SUCCESS)
	return UNUR_ERR_GEN_CONDITION;
      for (k=0; k<TABL_N_RUN_ARS; k++)
	_unur_sample_cont(gen);
    }
  }
  return UNUR_SUCCESS;
} 
struct unur_tabl_interval *
_unur_tabl_run_equalarearule( struct unur_par *par, struct unur_gen *gen, 
			      struct unur_tabl_interval *iv_slope )
{
  struct unur_tabl_interval *iv, *iv_last;
  double bar_area, x, fx;
  double slope;
  CHECK_NULL(par,NULL);       COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  CHECK_NULL(gen,NULL);       COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  CHECK_NULL(iv_slope,NULL);  COOKIE_CHECK(iv_slope,CK_TABL_IV,NULL);
  if (GEN->n_ivs >= GEN->max_ivs) 
    return NULL;
  iv = iv_slope;        
  iv_last = iv_slope;   
  bar_area = DISTR.area * PAR->area_fract;
  while (_unur_FP_greater(iv->Ahat, bar_area)) {
    slope = (iv->xmax > iv->xmin) ? 1. : -1.;
    x = iv->xmax - slope * bar_area / iv->fmax;
    fx = PDF(x);
    if (! (_unur_isfinite(fx) && fx >= 0.)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      return NULL;
    }
    switch (_unur_tabl_split_interval( gen, iv, x, fx, TABL_VARFLAG_SPLIT_POINT )) {
    case UNUR_SUCCESS:  
      if (slope > 0.) {
	if (iv_last == iv_slope)
	  iv_last = iv->next;
      }
      else { 
	iv = iv->next; break;
      }
      break;
    case UNUR_ERR_SILENT: 
      break; 
    default: 
      return NULL;
    }
    if (GEN->n_ivs >= GEN->max_ivs) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"split A stopped, maximal number of intervals reached.");
      break;
    }
  }
  return ((iv->xmax > iv->xmin) ? iv_last : iv);
} 
int
_unur_tabl_run_dars( struct unur_gen *gen )
{
  struct unur_tabl_interval *iv;
  double Atot, Asqueezetot;    
  double Alimit;               
  int n_splitted = 1;          
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;
  Atot = 0.;            
  Asqueezetot = 0.;     
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }
  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_ivs < GEN->max_ivs) ) {
    if (GEN->n_ivs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_ivs );
    else
      Alimit = 0.; 
    n_splitted = 0;
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
      if (GEN->n_ivs >= GEN->max_ivs)
	break;
      if ((iv->Ahat - iv->Asqueeze) <= Alimit) 
	continue;  
      switch (_unur_tabl_split_interval( gen, iv, 0., 0., TABL_VARFLAG_SPLIT_ARC )) {
      case UNUR_SUCCESS:  
      case UNUR_ERR_SILENT: 
	++n_splitted;
	break; 
      default: 
	return UNUR_ERR_GEN_DATA;
      }
    }
    if (n_splitted == 0) {
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: no intervals could be splitted.");
      break;
    }
  }
  if ( GEN->max_ratio * GEN->Atotal > GEN->Asqueeze ) {
    if ( GEN->n_ivs >= GEN->max_ivs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    GEN->max_ivs = GEN->n_ivs;
  }
  return UNUR_SUCCESS;
} 
int
_unur_tabl_split_interval( struct unur_gen *gen,
			   struct unur_tabl_interval *iv_old, 
			   double x, double fx, 
			   unsigned split_mode )
{
  struct unur_tabl_interval *iv_new;
  double A_hat_old, A_squ_old;
  CHECK_NULL(gen,UNUR_ERR_NULL);     COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_old,UNUR_ERR_NULL);  COOKIE_CHECK(iv_old,CK_TABL_IV,UNUR_ERR_COOKIE);
  switch( split_mode ) {
  case TABL_VARFLAG_SPLIT_POINT:    
    break;
  case TABL_VARFLAG_SPLIT_MEAN:     
    x = 0.5 * (iv_old->xmin + iv_old->xmax); 
    fx = PDF(x);
    break;
  case TABL_VARFLAG_SPLIT_ARC:      
    x = _unur_arcmean(iv_old->xmin, iv_old->xmax); 
    fx = PDF(x);
    break;
  default: 
    _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    break;
  }
  if (! (_unur_isfinite(fx) && fx >= 0.)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return UNUR_ERR_GEN_DATA;
  }
  if (_unur_FP_greater(fx,iv_old->fmax) || _unur_FP_less(fx,iv_old->fmin)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF not monotone in slope");
    return UNUR_ERR_GEN_DATA;
  }
  A_hat_old = iv_old->Ahat;
  A_squ_old = iv_old->Asqueeze;
  if (fx <= 0.) {
    if (iv_old->fmin > 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not monotone in slope");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_old->xmin = x;
    iv_old->Ahat = fabs(iv_old->xmax - iv_old->xmin) * iv_old->fmax;
    GEN->Atotal += iv_old->Ahat - A_hat_old;
    if (!_unur_isfinite(GEN->Atotal)) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_INF;
    }
    return UNUR_ERR_SILENT;
  }
  iv_new = _unur_xmalloc(sizeof(struct unur_tabl_interval));
  ++(GEN->n_ivs);
  COOKIE_SET(iv_new,CK_TABL_IV);
  if (iv_old->xmax > iv_old->xmin) {
    iv_new->xmax  = iv_old->xmax;  
    iv_new->fmax = iv_old->fmax;
    iv_old->xmax  = iv_new->xmin = x; 
    iv_old->fmax = iv_new->fmin = fx; 
  }
  else {
    iv_new->xmin  = iv_old->xmin;  
    iv_new->fmin = iv_old->fmin;
    iv_old->xmin  = iv_new->xmax = x; 
    iv_old->fmin = iv_new->fmax = fx; 
  }
  iv_new->Ahat     = fabs(iv_new->xmax - iv_new->xmin) * iv_new->fmax;
  iv_new->Asqueeze = fabs(iv_new->xmax - iv_new->xmin) * iv_new->fmin;
  iv_old->Ahat     = fabs(iv_old->xmax - iv_old->xmin) * iv_old->fmax;
  iv_old->Asqueeze = fabs(iv_old->xmax - iv_old->xmin) * iv_old->fmin;
  GEN->Atotal += iv_old->Ahat + iv_new->Ahat - A_hat_old;
  GEN->Asqueeze += iv_old->Asqueeze + iv_new->Asqueeze - A_squ_old;
  iv_new->next = iv_old->next;
  iv_old->next = iv_new;
  if (! (_unur_isfinite(GEN->Atotal) && _unur_isfinite(GEN->Asqueeze)) ) {
    _unur_error(gen->genid,UNUR_ERR_INF,"hat unbounded");
    return UNUR_ERR_INF;
  }
  return UNUR_SUCCESS;
} 
int
_unur_tabl_make_guide_table( struct unur_gen *gen )
{
  struct unur_tabl_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  if (!GEN->guide) {
    int max_guide_size = (GEN->guide_factor > 0.) ? (GEN->max_ivs * GEN->guide_factor) : 1;
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tabl_interval*) );
  }
  Acum = 0.;            
  Asqueezecum = 0.;     
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;
  GEN->guide_size = GEN->n_ivs;
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      if( iv->next != NULL )    
        iv = iv->next;
      else {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = iv;
    Acum += Astep;
  }
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = iv;
  if (! (_unur_isfinite(GEN->Atotal) && _unur_isfinite(GEN->Asqueeze)
	 && GEN->Atotal > 0.
	 && (!_unur_FP_less(GEN->Atotal,DISTR.area) || !(gen->distr->set & UNUR_DISTR_SET_PDFAREA)) ) 
      ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"sum of areas not valid");
    return UNUR_ERR_GEN_DATA;
  }
  return UNUR_SUCCESS;
} 
