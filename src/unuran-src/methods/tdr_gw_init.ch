/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

int
_unur_tdr_gw_starting_intervals( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv, *iv_new, *iv_tmp; 
  double x,fx;              
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  for( iv=GEN->iv; iv->next != NULL; ) {
    switch (_unur_tdr_gw_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:      
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:      
      break;
    case UNUR_ERR_SILENT:   
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN->n_ivs);
      if (iv->next==NULL) {
	iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
	iv->Acum = UNUR_INFINITY;
      }
      else
	iv->next->prev = iv;
      continue;
    default:     
      return UNUR_ERR_GEN_CONDITION;
    }
    x = _unur_arcmean(iv->x,iv->next->x);  
    fx = PDF(x);
    if (GEN->n_ivs >= GEN->max_ivs) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  
    if (fx <= 0.) {
      if (iv->fx <= 0.) {
	iv_new->next = iv->next;
	free(iv); 
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	iv_new->prev = NULL;
	iv = iv_new;
      }
      else if (iv->next->fx <= 0.) {
	free(iv->next);
	--(GEN->n_ivs);	
	iv->next = iv_new;
	iv_new->prev = iv;
      }
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	free(iv_new);
	return UNUR_ERR_GEN_CONDITION;
      }
    }
    else {
      iv_new->prev = iv;
      iv_new->next = iv->next;
      iv->next->prev = iv_new;
      iv->next = iv_new;
    }
  }
  return UNUR_SUCCESS;
} 
int
_unur_tdr_gw_dars( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv, *iv_next;
  double Alimit;               
  double x0, x1;               
  double xsp, fxsp;            
  double xAhatl, xAhatr, xAsqueeze;
  int rule;                    
  int n_splitted;              
  int splitted;                
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_ivs < GEN->max_ivs) ) {
    if (GEN->n_ivs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_ivs );
    else
      Alimit = 0.; 
    n_splitted = 0;
    for (iv = GEN->iv; iv->next != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
      if (GEN->n_ivs >= GEN->max_ivs)
	break;
      if ((iv->Ahat - iv->Asqueeze) <= Alimit) 
	continue;  
      iv_next = iv->next;
      x0 = iv->x;
      x1 = iv->next->x;
      for (rule = GEN->darsrule; rule <= 3; rule++) {
	switch (rule) {
	case 1:   
	  if ( _unur_FP_is_minus_infinity(x0) ||
	       _unur_FP_is_infinity(x1) ||
	       _unur_FP_approx(x0,x1) )
	    continue;  
	  xAhatl = _unur_tdr_interval_xxarea( gen, iv, iv->dTfx, iv->ip);
	  xAhatr = _unur_tdr_interval_xxarea( gen, iv->next, iv->next->dTfx, iv->ip);
	  if (iv->Asqueeze > 0.)
	    xAsqueeze = (iv->Tfx > iv->next->Tfx)
	      ? _unur_tdr_interval_xxarea( gen, iv, iv->sq, x1)
	      : _unur_tdr_interval_xxarea( gen, iv->next, iv->sq, x0);
	  else  
	    xAsqueeze = 0.;
	  if (! (_unur_isfinite(xAhatl) && _unur_isfinite(xAhatr) && _unur_isfinite(xAsqueeze))
	      || _unur_FP_equal(iv->Ahat,iv->Asqueeze) )
	    continue;  
	  xsp = (xAhatl+xAhatr-xAsqueeze) / (iv->Ahat - iv->Asqueeze);
	  break;
	case 2:   
	  xsp = _unur_arcmean(x0,x1);
	  break;
	case 3:  
	  if (_unur_FP_is_minus_infinity(x0) || _unur_FP_is_infinity(x1))
	    continue;  
	  xsp = 0.5 * (x0 + x1);
	  break;
	default:   
	  continue;
	}
	fxsp = PDF(xsp);
	splitted = _unur_tdr_gw_interval_split(gen, iv, xsp, fxsp);
	if (splitted==UNUR_SUCCESS || splitted==UNUR_ERR_INF) {
	  if (splitted==UNUR_SUCCESS) ++n_splitted;
	  if (iv->next != iv_next)
	    iv = iv->next;
	  break;
	}
	else if (splitted!=UNUR_ERR_SILENT) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  return UNUR_ERR_GEN_CONDITION;
	}
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
_unur_tdr_gw_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
{
  double Ahatl;    
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  CHECK_NULL(iv->next,UNUR_ERR_NULL);  COOKIE_CHECK(iv->next,CK_TDR_IV,UNUR_ERR_COOKIE); 
  if ( _unur_tdr_tangent_intersection_point(gen,iv,&(iv->ip))!=UNUR_SUCCESS )
    return UNUR_ERR_GEN_CONDITION;
  if (iv->Tfx > -UNUR_INFINITY && iv->next->Tfx > -UNUR_INFINITY) {
    if (_unur_FP_approx(iv->x, iv->next->x) )
      return UNUR_ERR_SILENT;   
    iv->sq = (iv->next->Tfx - iv->Tfx) / (iv->next->x - iv->x);
    if ( ( (iv->sq > iv->dTfx       && (!_unur_FP_approx(iv->sq,iv->dTfx)) ) || 
	   (iv->sq < iv->next->dTfx && (!_unur_FP_approx(iv->sq,iv->next->dTfx)) ) )
	 && iv->next->dTfx < UNUR_INFINITY ) {
      if ( !_unur_iszero(iv->sq) && !_unur_iszero(iv->dTfx) && !_unur_iszero(iv->next->dTfx) ) {
      	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. PDF not T-concave!");
      	return UNUR_ERR_GEN_CONDITION;
      }
    }
    iv->Asqueeze = (iv->Tfx > iv->next->Tfx)
      ? _unur_tdr_interval_area( gen, iv, iv->sq, iv->next->x)
      : _unur_tdr_interval_area( gen, iv->next, iv->sq, iv->x);
    if (!_unur_isfinite(iv->Asqueeze))
      iv->Asqueeze = 0.;
  }
  else {  
    iv->sq = 0.;
    iv->Asqueeze = 0.;
  }
  Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->ip);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, iv->ip);
  if (! (_unur_isfinite(Ahatl) && _unur_isfinite(iv->Ahatr)) )
    return UNUR_ERR_INF;
  iv->Ahat = iv->Ahatr + Ahatl;
  if ( iv->Asqueeze > iv->Ahat && !(_unur_FP_approx(iv->Asqueeze, iv->Ahat)) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"A(squeeze) > A(hat). PDF not T-concave!");
    return UNUR_ERR_GEN_CONDITION; 
  }
  return UNUR_SUCCESS;
} 
int
_unur_tdr_gw_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv_oldl, double x, double fx )
{
  struct unur_tdr_interval *iv_newr;  
  struct unur_tdr_interval iv_bak;    
  int success, success_r;
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_oldl,UNUR_ERR_NULL);  COOKIE_CHECK(iv_oldl,CK_TDR_IV,UNUR_ERR_COOKIE);
  if (!_unur_isfinite(x)) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not finite (skipped)");
    return UNUR_ERR_SILENT;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_gw_debug_split_start( gen,iv_oldl,x,fx );
#endif
  if (x < iv_oldl->x || x > iv_oldl->next->x) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }
  if ( (GEN->n_ivs * (iv_oldl->Ahat - iv_oldl->Asqueeze) / (GEN->Atotal - GEN->Asqueeze))
       < GEN->bound_for_adding)
    return UNUR_ERR_SILENT;
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }
  memcpy(&iv_bak, iv_oldl, sizeof(struct unur_tdr_interval));
  if (fx <= 0.) {
    if (iv_oldl->fx <= 0.) {
      iv_oldl->x = x;
    }
    else if (iv_oldl->next->fx <= 0.) {
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    success = _unur_tdr_gw_interval_parameter(gen, iv_oldl);
    iv_newr = NULL;
  }
  else {
    iv_newr = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_newr == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    iv_newr->prev = iv_oldl;
    iv_newr->next = iv_oldl->next;
    iv_oldl->next->prev = iv_newr;
    iv_oldl->next = iv_newr;
    success   = _unur_tdr_gw_interval_parameter(gen, iv_oldl);
    success_r = _unur_tdr_gw_interval_parameter(gen, iv_newr);
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  if (success!=UNUR_SUCCESS) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
    memcpy(iv_oldl, &iv_bak, sizeof(struct unur_tdr_interval));
    if (iv_oldl->next)
      iv_oldl->next->prev = iv_oldl;
    if (iv_newr) {
      --(GEN->n_ivs); 
      free( iv_newr );
    }
  return success;
  }
  GEN->Atotal   = ( GEN->Atotal - iv_bak.Ahat
		   + iv_oldl->Ahat + ((iv_newr) ? iv_newr->Ahat : 0.) );
  GEN->Asqueeze = ( GEN->Asqueeze - iv_bak.Asqueeze
		   + iv_oldl->Asqueeze + ((iv_newr) ? iv_newr->Asqueeze : 0. ) );
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TDR_DEBUG_SPLIT)
    _unur_tdr_gw_debug_split_stop( gen,iv_oldl,iv_newr );
#endif
  if (GEN->Atotal <= 1.e10 * DBL_MIN) {
    _unur_error(gen->genid,UNUR_ERR_ROUNDOFF,"error below hat (almost) 0");
    return UNUR_ERR_ROUNDOFF;
  }
  return UNUR_SUCCESS;
} 
