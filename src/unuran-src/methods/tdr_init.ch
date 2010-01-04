/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include "tdr_gw_init.ch"
#include "tdr_ps_init.ch"
struct unur_gen *
_unur_tdr_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_TDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TDR_PAR,NULL);
  gen = _unur_tdr_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }
  _unur_par_free(par);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_tdr_debug_init_start(gen);
#endif
  if (_unur_tdr_make_gen( gen ) != UNUR_SUCCESS) {
    _unur_tdr_free(gen); return NULL;
  }
  if (GEN->Atotal <= 0. || !_unur_isfinite(GEN->Atotal)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_tdr_free(gen);
    return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_tdr_debug_init_finished(gen);
#endif
  gen->status = UNUR_SUCCESS;
  return gen;
} 
int
_unur_tdr_make_gen( struct unur_gen *gen )
{ 
  int i,k;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  if (_unur_tdr_starting_cpoints(gen)!=UNUR_SUCCESS) return UNUR_FAILURE;
  if (_unur_tdr_starting_intervals(gen)!=UNUR_SUCCESS) return UNUR_FAILURE;
  if (GEN->n_ivs > GEN->max_ivs) GEN->max_ivs = GEN->n_ivs;
  if (gen->variant & TDR_VARFLAG_USEDARS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TDR_DEBUG_DARS) {
      _unur_tdr_make_guide_table(gen);
      _unur_tdr_debug_dars_start(gen);
    }
#endif
    for (i=0; i<3; i++) {
      if (_unur_tdr_run_dars(gen)!=UNUR_SUCCESS) return UNUR_FAILURE;
      _unur_tdr_make_guide_table(gen);
      if (GEN->n_ivs < GEN->max_ivs) {
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_tdr_debug_dars_finished(gen);
#endif
  }
  else { 
    _unur_tdr_make_guide_table(gen);
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_tdr_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TDR_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_tdr_gen) );
  COOKIE_SET(gen,CK_TDR_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  if (_unur_iszero(PAR->c_T))
    gen->variant = (gen->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_LOG;
  else if (_unur_FP_same(PAR->c_T, -0.5))
    gen->variant = (gen->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_SQRT;
  else
    gen->variant = (gen->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_POW;
  if ((gen->variant & TDR_VARMASK_T) == TDR_VAR_T_POW) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"c != 0. and c != -0.5 not implemented!");
    _unur_generic_free(gen);
    return NULL;
  }
  SAMPLE = _unur_tdr_getSAMPLE(gen);
  gen->destroy = _unur_tdr_free;
  gen->clone = _unur_tdr_clone;
  gen->reinit = _unur_tdr_reinit;
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;
  GEN->guide_factor = PAR->guide_factor; 
  GEN->c_T = PAR->c_T;                
  GEN->darsfactor = PAR->darsfactor;  
  GEN->darsrule = PAR->darsrule;      
  GEN->max_ivs = _unur_max(2*PAR->n_starting_cpoints,PAR->max_ivs);  
#ifdef UNUR_ENABLE_INFO
  GEN->max_ivs_info = PAR->max_ivs;   
#endif
  GEN->max_ratio = PAR->max_ratio;    
  GEN->bound_for_adding = PAR->bound_for_adding;
  if ( (gen->distr->set & UNUR_DISTR_SET_CENTER) ||
       (gen->distr->set & UNUR_DISTR_SET_MODE) ) {
    GEN->center = unur_distr_cont_get_center(gen->distr);
    GEN->center = _unur_max(GEN->center,DISTR.BD_LEFT);
    GEN->center = _unur_min(GEN->center,DISTR.BD_RIGHT);
    gen->set |= TDR_SET_CENTER;
  }
  else {
    GEN->center = 0.;
    gen->variant &= ~TDR_VARFLAG_USECENTER;
  }
  if ( !(gen->distr->set & UNUR_DISTR_SET_MODE)
       || (DISTR.mode < DISTR.BD_LEFT)
       || (DISTR.mode > DISTR.BD_RIGHT))
    gen->variant = gen->variant & (~TDR_VARFLAG_USEMODE);
  GEN->n_starting_cpoints = PAR->n_starting_cpoints;
  if (PAR->starting_cpoints) {
    GEN->starting_cpoints = _unur_xmalloc( PAR->n_starting_cpoints * sizeof(double) );
    memcpy( GEN->starting_cpoints, PAR->starting_cpoints, PAR->n_starting_cpoints * sizeof(double) );
  }
  else {
    GEN->starting_cpoints = NULL;
  }
  GEN->percentiles = NULL;
  if (gen->set & TDR_SET_N_PERCENTILES)
    unur_tdr_chg_reinit_percentiles( gen, PAR->n_percentiles, PAR->percentiles );
  GEN->retry_ncpoints = PAR->retry_ncpoints;   
  GEN->Umin = 0.;
  GEN->Umax = 1.;
  if (!(gen->set & TDR_SET_USE_DARS) && !PAR->starting_cpoints)
    gen->variant |= TDR_VARFLAG_USEDARS;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_tdr_info;
#endif
  return gen;
} 
int
_unur_tdr_reinit( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv,*next;
  double *bak_cpoints;
  int bak_n_cpoints;
  int i;
  int n_trials;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );
  n_trials = 1;
  if (gen->set & TDR_SET_N_PERCENTILES) {
    if (GEN->starting_cpoints==NULL || (GEN->n_starting_cpoints != GEN->n_percentiles)) {
      GEN->n_starting_cpoints = GEN->n_percentiles;
      GEN->starting_cpoints = _unur_xrealloc( GEN->starting_cpoints, GEN->n_percentiles * sizeof(double));
    }
    for (i=0; i<GEN->n_percentiles; i++) {
      GEN->starting_cpoints[i] = unur_tdr_eval_invcdfhat( gen, GEN->percentiles[i], NULL, NULL, NULL );
      if (!_unur_isfinite(GEN->starting_cpoints[i])) 
	n_trials = 2;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TDR_DEBUG_REINIT)
    _unur_tdr_debug_reinit_start(gen);
#endif
  bak_n_cpoints = GEN->n_starting_cpoints;
  bak_cpoints = GEN->starting_cpoints;
  for (;; ++n_trials) {
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
    GEN->iv = NULL;
    GEN->n_ivs = 0;
    GEN->Atotal = 0.;
    GEN->Asqueeze = 0.;
    if (n_trials > 2) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points for reinit");
      GEN->n_starting_cpoints = bak_n_cpoints;
      GEN->starting_cpoints = bak_cpoints;
      return UNUR_FAILURE;
    }
    if (n_trials > 1) {
      GEN->n_starting_cpoints = GEN->retry_ncpoints;
      GEN->starting_cpoints = NULL;
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug & TDR_DEBUG_REINIT)
	_unur_tdr_debug_reinit_retry(gen);
#endif
    }
    if (_unur_tdr_make_gen( gen ) != UNUR_SUCCESS)
      continue;
    if (GEN->Atotal <= 0.)
      continue;
    break;
  }
  if (n_trials > 1) {
    GEN->n_starting_cpoints = bak_n_cpoints;
    GEN->starting_cpoints = bak_cpoints;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TDR_DEBUG_REINIT)
    if (gen->debug) _unur_tdr_debug_reinit_finished(gen);
#endif
  SAMPLE = _unur_tdr_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_tdr_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_tdr_gen*)clone->datap)
  struct unur_gen *clone;
  struct unur_tdr_interval *iv,*next, *clone_iv, *clone_prev;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = next) {
    clone_iv = _unur_xmalloc( sizeof(struct unur_tdr_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tdr_interval) );
    if (clone_prev == NULL) {
      CLONE->iv = clone_iv;
      clone_iv->prev = NULL;
    }
    else {
      clone_prev->next = clone_iv;
      clone_iv->prev = clone_prev;
    }
    next = iv->next;
    clone_prev = clone_iv;
  }
  if (clone_iv) clone_iv->next = NULL;
  if (GEN->starting_cpoints) {
    CLONE->starting_cpoints = _unur_xmalloc( GEN->n_starting_cpoints * sizeof(double) );
    memcpy( CLONE->starting_cpoints, GEN->starting_cpoints, GEN->n_starting_cpoints * sizeof(double) );
  }
  if (GEN->percentiles) {
    CLONE->percentiles = _unur_xmalloc( GEN->n_percentiles * sizeof(double) );
    memcpy( CLONE->percentiles, GEN->percentiles, GEN->n_percentiles * sizeof(double) );
  }
  CLONE->guide = NULL;
  _unur_tdr_make_guide_table(clone);
  return clone;
#undef CLONE
} 
void
_unur_tdr_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_TDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  SAMPLE = NULL;   
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_tdr_debug_free(gen);
#endif
  {
    struct unur_tdr_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }
  if (GEN->starting_cpoints) 
    free (GEN->starting_cpoints);
  if (GEN->percentiles) 
    free (GEN->percentiles);
  if (GEN->guide)  free(GEN->guide);
  _unur_generic_free(gen);
} 
int
_unur_tdr_starting_cpoints( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, fx, fx_last;
  int use_center, use_mode, is_mode, was_mode;
  int i, is_increasing;
  double extra_cpoint;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  use_mode = (gen->variant & TDR_VARFLAG_USEMODE) ? TRUE : FALSE;
  use_center = (!use_mode && (gen->variant & TDR_VARFLAG_USECENTER)) ? TRUE : FALSE;
  extra_cpoint = use_mode ? DISTR.mode : (use_center ? GEN->center : 0. );
  GEN->n_ivs = 0;
  if (!GEN->starting_cpoints) {
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT  - GEN->center);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT - GEN->center);
    diff_angle = (right_angle-left_angle) / (GEN->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   
  x = DISTR.BD_LEFT;
  if (use_mode && DISTR.mode <= x) {
    is_mode = TRUE;
    use_mode = FALSE;  
    is_increasing = FALSE;
  }
  else if (use_center && GEN->center <= x) {
    is_mode = FALSE;
    use_center = FALSE;     
    is_increasing = TRUE;   
  }
  else {
    is_mode = FALSE;
    is_increasing = TRUE;
  }
  fx = fx_last = _unur_FP_is_minus_infinity(x) ? 0. : PDF(x);
  iv = GEN->iv = _unur_tdr_interval_new( gen, x, fx, is_mode );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  
  iv->prev = NULL;
  for( i=0; i<=GEN->n_starting_cpoints; i++ ) {
    was_mode = is_mode;
    if (i < GEN->n_starting_cpoints) {
      if (GEN->starting_cpoints) {   
	x = GEN->starting_cpoints[i];
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	angle += diff_angle;
	x = tan( angle ) + GEN->center;
      }
    }
    else {
      x = DISTR.BD_RIGHT;
    }
    if ((use_mode || use_center) && x >= extra_cpoint) {
      is_mode = use_mode;              
      use_center = use_mode = FALSE;   
      if (x>extra_cpoint) {
	x = extra_cpoint;     
	--i;              
	if (!GEN->starting_cpoints)
	  angle -= diff_angle; 
      }
    }
    else
      is_mode = FALSE;
    fx = _unur_FP_is_infinity(x) ? 0. : PDF(x);
    if (!is_increasing && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal!");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (is_mode && fx < fx_last * (1.-DBL_EPSILON)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"mode -> ignore");
      continue;
    }
    if (was_mode && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"mode");
      return UNUR_ERR_GEN_DATA;
    }
    if (fx <= 0. && fx_last <= 0.) {
      if (is_increasing) {
	if (i<GEN->n_starting_cpoints) {
	  iv->x = x;  
	  continue;   
	}
      }
      else
	break;
    }
    iv->next = _unur_tdr_interval_new( gen, x, fx, is_mode );
    if (iv->next == NULL) return UNUR_ERR_GEN_DATA;  
    iv->next->prev = iv;
    iv = iv->next;
    if (is_increasing && fx < fx_last)
      is_increasing = 0;
    fx_last = fx;
  }
  iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
  iv->Acum = INFINITY;
  iv->ip = iv->x;
  iv->fip = iv->fx;
  iv->next = NULL;         
  --(GEN->n_ivs);           
  return UNUR_SUCCESS;
} 
int
_unur_tdr_starting_intervals( struct unur_gen *gen )
{
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    
    return _unur_tdr_gw_starting_intervals(gen);
  case TDR_VARIANT_PS:    
  case TDR_VARIANT_IA:    
    return _unur_tdr_ps_starting_intervals(gen);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
} 
int
_unur_tdr_run_dars( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv;
  double Atot, Asqueezetot;    
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;
  Atot = 0.;            
  Asqueezetot = 0.;     
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }
  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    
    return _unur_tdr_gw_dars(gen);
  case TDR_VARIANT_PS:    
  case TDR_VARIANT_IA:    
    return _unur_tdr_ps_dars(gen);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
} 
struct unur_tdr_interval *
_unur_tdr_interval_new( struct unur_gen *gen, double x, double fx, int is_mode )
{
  struct unur_tdr_interval *iv;
  double dfx;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);
  if (fx<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return NULL;
  }
  if (_unur_FP_is_infinity(fx)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return NULL;
  }
  iv = _unur_xmalloc( sizeof(struct unur_tdr_interval) );
  iv->next = NULL; 
  ++(GEN->n_ivs);   
  COOKIE_SET(iv,CK_TDR_IV);
  iv->Acum = iv->Ahat = iv->Ahatr = iv->Asqueeze = 0.;
  iv->ip = iv->fip = iv->sq = 0.;
  iv->x = x;              
  iv->fx = fx;            
  if (fx<=0.) {           
    iv->Tfx = -INFINITY;  
    iv->dTfx = INFINITY;  
    return iv;
  }
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    iv->Tfx = log(fx);
    if (is_mode) { 
      iv->dTfx = 0.; break; 
    }
    if (_unur_cont_have_dlogPDF(gen->distr)) {
      iv->dTfx = dlogPDF(x); break; 
    }
    else {
      dfx = dPDF(x);
      if (_unur_iszero(dfx))
	iv->dTfx = 0.;
      else
	iv->dTfx = (1./fx * dfx);   
    }
    break;
  case TDR_VAR_T_SQRT:
    iv->Tfx = -1./sqrt(fx);
    if (is_mode) { 
      iv->dTfx = 0.; break; }
    if (_unur_cont_have_dlogPDF(gen->distr)) {
      iv->dTfx = -0.5 * iv->Tfx * dlogPDF(x);
      break;
    }
    else {
      dfx = dPDF(x);
      if (_unur_iszero(dfx))
	iv->dTfx = 0.;
      else
	iv->dTfx = (dfx<0.) ? -exp( -M_LN2 - 1.5*log(fx) + log(-dfx))
	  : exp( -M_LN2 - 1.5*log(fx) + log(dfx));
    }
    break;
  case TDR_VAR_T_POW:
    break;
  }
  if ( !(iv->dTfx > -INFINITY))
    iv->dTfx = INFINITY;
  return iv;
} 
int
_unur_tdr_tangent_intersection_point( struct unur_gen *gen, struct unur_tdr_interval *iv, double *ipt )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  if ( iv->dTfx > 1.e+140 ) {
    *ipt = iv->x;        
    return UNUR_SUCCESS; 
  }
  if ( iv->next->dTfx < -1.e+140 || _unur_FP_is_infinity(iv->next->dTfx)) {
    *ipt = iv->next->x;   
    return UNUR_SUCCESS; 
  }
  if ( _unur_FP_less( iv->dTfx, iv->next->dTfx ) ) {
    if ( fabs(iv->dTfx) < DBL_EPSILON * fabs(iv->next->dTfx) ) {
      *ipt = iv->x;        
      iv->dTfx = INFINITY;
      return UNUR_SUCCESS; 
    }
    else if ( fabs(iv->next->dTfx) < DBL_EPSILON * fabs(iv->dTfx) ) {
      *ipt = iv->next->x;   
      iv->next->dTfx = INFINITY;
      return UNUR_SUCCESS; 
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). PDF not T-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  if (_unur_FP_approx(iv->dTfx, iv->next->dTfx)) {
    *ipt = 0.5 * (iv->x + iv->next->x);
    return UNUR_SUCCESS;
  }
  *ipt = ( (iv->next->Tfx - iv->Tfx - iv->next->dTfx * iv->next->x + iv->dTfx * iv->x) / 
	   (iv->dTfx - iv->next->dTfx) );
  if (_unur_FP_less(*ipt, iv->x) || _unur_FP_greater(*ipt, iv->next->x))
    *ipt = 0.5 * (iv->x + iv->next->x);
  return UNUR_SUCCESS;
} 
double
_unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
{
  double area = 0.;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
  if (_unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x))
    return 0.;
  if (_unur_FP_same(x, iv->x))
    return 0.;
  if ( _unur_FP_is_infinity(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   
    return INFINITY;
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:   
    if (!_unur_iszero(slope)) {                         
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	area = iv->fx / slope;
      else {
	double t = slope * (x - iv->x);
	if (fabs(t) > 1.e-6) {
	  if (t > MAXLOG / 10.) {
	    double xdiff = (x>iv->x) ? x - iv->x : iv->x - x;
	    area = exp( log(iv->fx) + log(xdiff) + t - log(t) );
	  }
	  else {
	    area = iv->fx * (x - iv->x) * ( exp(t) - 1. ) / t;
	  }
	}
	else if (fabs(t) > 1.e-8)
	  area = iv->fx * (x - iv->x) * (1. + t/2. + t*t/6.);
	else
	  area = iv->fx * (x - iv->x) * (1. + t/2.);
      }
    }
    else { 
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	return INFINITY;
      else
	area = iv->fx * (x - iv->x);
    }
    break;
  case TDR_VAR_T_SQRT:   
    if (!_unur_iszero(slope)) {
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	area = 1. / ( iv->Tfx * slope );
      else {
	double hx = iv->Tfx + slope * (x - iv->x);
	if (hx>=0.)
	  return INFINITY; 
	else
	  area = (x - iv->x) / ( iv->Tfx * hx );
      }
    }
    else { 
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	return INFINITY;
      else
	area = iv->fx * (x - iv->x);
    }
    break;
  case TDR_VAR_T_POW:   
    break;
  }
  return ( (area<0.) ? -area : area );
} 
double
_unur_tdr_interval_xxarea( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
{
  double ev = 0.;
  double hx,u;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
  if (_unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x))
    return 0.;
  if (_unur_FP_same(x, iv->x))
    return 0.;
  if ( _unur_FP_is_infinity(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   
    return INFINITY;
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:    
    if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x)) {
      ev = iv->fx / (slope*slope) * (1-slope*iv->x);
    }
    else {
      u = (x-iv->x) * slope;
      if (fabs(u) > 1.e-6) {
	ev = iv->fx / (slope*slope) * (exp(u)*(slope*x-1.) - slope*iv->x + 1.);
      }
      else {
	ev = 0.5 * (x+iv->x);
	if (fabs(u) > 0) {
	  ev += 1./6. * (2.*x+iv->x) * u;
	  ev += 1./24. * (3.*x+iv->x) * u * u;
	}
	ev *= iv->fx * (x-iv->x);
      }
    }
    break;
  case TDR_VAR_T_SQRT:    
    if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
      return INFINITY;
    hx = iv->Tfx + slope * (x - iv->x);
    if (hx >= 0.)
      return INFINITY; 
    u = (x-iv->x) * slope / iv->Tfx;
    if (fabs(u) > 1.e-6) {
      ev = ( iv->x / (slope * iv->Tfx) - x / (slope * hx)
	     + log( hx / iv->Tfx ) / (slope*slope) );
    }
    else {
      ev = 0.5 * (x+iv->x);
      if (fabs(u) > 0) {
	ev -= 1./3. * (2.*x+iv->x) * u;
	ev += 1./4. * (3.*x+iv->x) * u * u;
      }
      ev *= iv->fx * (x-iv->x);
    }
    break;
  case TDR_VAR_T_POW:    
    break;
  }
  return ((x>iv->x) ? ev : -ev);
} 
double
_unur_tdr_eval_intervalhat( struct unur_gen *gen, struct unur_tdr_interval *iv, double x )
{
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
  if ( _unur_FP_is_minus_infinity(iv->Tfx) || _unur_FP_is_infinity(iv->dTfx) )
    return INFINITY;
  if ( _unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x) ||
       _unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x) )
    return 0.;
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    return (iv->fx * exp( iv->dTfx * (x - iv->x) ));
  case TDR_VAR_T_SQRT:
    {
      double hx = iv->Tfx + iv->dTfx * (x - iv->x);
      return ((hx<0.) ? 1./(hx*hx) : INFINITY);
    }
  case TDR_VAR_T_POW:
    return INFINITY;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }
} 
int
_unur_tdr_make_guide_table( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv;
  double Acum, Asqueezecum, Astep;
  int max_guide_size;
  int j;
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  if (!GEN->guide) {
    max_guide_size = (GEN->guide_factor > 0.) ? ((int)(GEN->max_ivs * GEN->guide_factor)) : 1;
    if (max_guide_size <= 0) max_guide_size = 1;   
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tdr_interval*) );
  }
  Acum = 0.;            
  Asqueezecum = 0.;     
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;
  GEN->guide_size = (int)(GEN->n_ivs * GEN->guide_factor);
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      iv = iv->next;
    if( iv->next == NULL ) {   
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = iv;
    Acum += Astep;
  }
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = iv;
  return UNUR_SUCCESS;
} 
