/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

int
_unur_tdr_ps_starting_intervals( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv, *iv_new, *iv_tmp; 
  double x,fx;              
  double lb, flb;           
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  iv = GEN->iv;
  lb = iv->x;
  flb = iv->fx;
  if (_unur_FP_is_infinity(iv->dTfx)) {
    GEN->iv = iv->next;
    GEN->iv->prev = NULL;
    free (iv);
    --(GEN->n_ivs);
    iv = GEN->iv;
  }
  iv->ip = lb;
  iv->fip = flb;
  while (iv) {
    if (iv->next == NULL) {
      if (!_unur_FP_is_infinity(iv->dTfx)) {
	iv->next = iv_new = _unur_tdr_interval_new( gen, iv->x, 0., FALSE );
	if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  
	iv_new->prev = iv;
	iv_new->ip = iv->x;
	iv_new->fip = iv->fx;
	iv->next->Asqueeze = iv->next->Ahat = iv->next->Ahatr = 0.;
	iv->Acum = UNUR_INFINITY;
	iv->next-> sq = 0.;
      }
      else
	break;
    }
    switch (_unur_tdr_ps_interval_parameter(gen, iv)) {
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
    if (GEN->n_ivs >= GEN->max_ivs) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (iv->Ahatr >= UNUR_INFINITY) {
      CHECK_NULL(iv->next,0);
      x = _unur_arcmean(iv->x,iv->next->ip);
      fx = PDF(x);
      iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  
      if (fx <= 0.) {
	if (iv->next->fx > 0.) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	  free(iv_new);
	  return UNUR_ERR_GEN_CONDITION;
	}
	free(iv->next);
	--(GEN->n_ivs);
	iv->next = iv_new;
	iv_new->prev = iv;
      }
      else {
	iv_new->prev = iv;
	iv_new->next = iv->next;
	iv->next->prev = iv_new;
	iv->next = iv_new;
      }
    }
    else {
      x = _unur_arcmean(iv->ip,iv->x);
      fx = PDF(x);
      iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  
      if (fx <= 0.) {
	if (iv->fx > 0.) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	  free(iv_new);
	  return UNUR_ERR_GEN_CONDITION;
	}
	iv_new->next = iv->next;
	iv_new->prev = iv->prev;
	iv_new->ip = iv->ip;
	iv_new->fip = iv->fip;
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	free(iv);
	iv = iv_new;
      }
      else {
	if (iv->prev) {
	  iv_tmp = iv->prev;
	  iv_new->prev = iv->prev;
	  iv_new->next = iv;
	  iv->prev->next = iv_new;
	  iv->prev = iv_new;
	  iv_new->ip = iv->ip;
	  iv = iv_tmp;
	}
	else { 
	  iv_new->ip = iv->ip;
	  iv_new->fip = iv->fip;
	  iv_new->prev = NULL;
	  iv_new->next = iv;
	  iv->prev = iv_new;
	  GEN->iv = iv_new;
	  iv = iv_new;
	}
      }
    }
  }
  return UNUR_SUCCESS;
} 
int
_unur_tdr_ps_dars( struct unur_gen *gen )
{
  struct unur_tdr_interval *iv, *iv_next;
  double Adiff;                
  double Alimit;               
  double squeeze_gw;           
  double Asqueeze_gw;          
  double Ahat_gw;              
  double x0, x1;               
  double xsp, fxsp;            
  double xAhatl, xAhatr, xAsqueeze_gw;
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
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
      if (GEN->n_ivs >= GEN->max_ivs)
	break;
      if (iv==GEN->iv) {
	x0 = iv->ip;       
	x1 = iv->x;        
	Adiff = iv->Ahat - iv->Ahatr; 
	Adiff *= 1. - iv->sq;         
      }
      else {
	x0 = iv->prev->x;       
	x1 = iv->x;             
	Adiff = ( (1. - iv->prev->sq) * iv->prev->Ahatr +
		  (1. - iv->sq) * (iv->Ahat - iv->Ahatr) );
      }
      if (Adiff <= Alimit) 
	continue;  
      for (rule = GEN->darsrule; rule <= 3; rule++) {
	switch (rule) {
	case 1:   
	  if ( _unur_FP_is_minus_infinity(x0) ||
	       _unur_FP_is_infinity(x1) ||
	       _unur_FP_approx(x0,x1) )
	    continue;  
	  xAhatl = ( (iv->prev == NULL) ? 0. : 
		     _unur_tdr_interval_xxarea( gen, iv->prev, iv->prev->dTfx, iv->ip) );
	  xAhatr = _unur_tdr_interval_xxarea( gen, iv, iv->dTfx, iv->ip);
	  Ahat_gw = ( ((iv->prev == NULL) ? 0. : iv->prev->Ahatr)
		      + (iv->Ahat - iv->Ahatr) );
	  if (iv->Asqueeze > 0. && iv->prev != NULL) {
	    squeeze_gw = (iv->Tfx - iv->prev->Tfx) / (iv->x - iv->prev->x);
	    xAsqueeze_gw = (iv->Tfx > iv->prev->Tfx)
	      ? _unur_tdr_interval_xxarea( gen, iv, squeeze_gw, iv->prev->x)
	      : _unur_tdr_interval_xxarea( gen, iv->prev, squeeze_gw, iv->x);
	    Asqueeze_gw = (iv->Tfx > iv->prev->Tfx)
	      ? _unur_tdr_interval_area( gen, iv, squeeze_gw, iv->prev->x)
	      : _unur_tdr_interval_area( gen, iv->prev, squeeze_gw, iv->x);
	    if (!_unur_isfinite(Asqueeze_gw))  Asqueeze_gw = 0.;
	  }
	  else { 
	    xAsqueeze_gw = 0.;
	    Asqueeze_gw = 0.;
	  }
	  if (! (_unur_isfinite(xAhatl) && _unur_isfinite(xAhatr) && _unur_isfinite(xAsqueeze_gw))
	      || _unur_FP_equal(Ahat_gw,Asqueeze_gw) )
	    continue;  
	  xsp = (xAhatl+xAhatr-xAsqueeze_gw) / (Ahat_gw - Asqueeze_gw);
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
	iv_next = iv->next;
	splitted = _unur_tdr_ps_interval_split(gen, ((xsp<iv->ip)?iv->prev:iv), xsp, fxsp);
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
_unur_tdr_ps_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
{
  double Ahatl;    
  double hxl, hxr; 
  double sq;       
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  if (_unur_tdr_tangent_intersection_point(gen,iv,&(iv->next->ip))!=UNUR_SUCCESS)
    return UNUR_ERR_GEN_CONDITION;
  iv->next->fip = _unur_FP_is_infinity(iv->next->ip) ? 0. : PDF(iv->next->ip);
  Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->ip);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->next->ip);
  if (! (_unur_isfinite(Ahatl) && _unur_isfinite(iv->Ahatr)) )
    return UNUR_ERR_INF;
  iv->Ahat = iv->Ahatr + Ahatl;
  hxl = _unur_tdr_eval_intervalhat(gen,iv,iv->ip);
  if (_unur_FP_greater(iv->fip, hxl) ) {
    if ( (iv->fip < 1.e-50) || _unur_FP_approx(iv->fip, hxl)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) might be < PDF(x)");
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) < PDF(x)");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  iv->sq = (_unur_FP_is_infinity(hxl) || hxl <= 0.) ? 0. : iv->fip / hxl;
  hxr = _unur_tdr_eval_intervalhat(gen,iv,iv->next->ip);
  if (_unur_FP_greater(iv->next->fip, hxr)) {
    if ((iv->next->fip < 1.e-50) || _unur_FP_approx(iv->next->fip, hxr)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) might be < PDF(x)");
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) < PDF(x)");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  sq = (_unur_FP_is_infinity(hxr) || hxr <= 0.) ? 0. : iv->next->fip / hxr;
  if (iv->sq > sq) iv->sq = sq;
  iv->Asqueeze = iv->Ahat * iv->sq;
  return UNUR_SUCCESS;
} 
int
_unur_tdr_ps_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv, double x, double fx )
{
  struct unur_tdr_interval *oldl, *oldr;  
  struct unur_tdr_interval *iv_new;       
  struct unur_tdr_interval oldl_bak, oldr_bak; 
  int success, success_r;
  CHECK_NULL(gen,UNUR_ERR_NULL); COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);  COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
  if (!_unur_isfinite(x)) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not finite (skipped)");
    return UNUR_ERR_SILENT;
  }
  if ( (GEN->n_ivs * (iv->Ahat - iv->Asqueeze) / (GEN->Atotal - GEN->Asqueeze))
       < GEN->bound_for_adding)
    return UNUR_ERR_SILENT;
  if (x < iv->ip || x > iv->next->ip) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }
  if (x < iv->x) {
    oldl = iv->prev;
    oldr = iv;
  }
  else {
    oldl = iv;
    oldr = iv->next;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_ps_debug_split_start( gen,oldl,oldr,x,fx );
#endif
  if (oldl) memcpy(&oldl_bak, oldl, sizeof(struct unur_tdr_interval));
  memcpy(&oldr_bak, oldr, sizeof(struct unur_tdr_interval));
  if (fx <= 0.) {
    if (oldr->fip <= 0. && oldl==NULL) {
      oldr->ip = x;
      oldr->fip = 0.;
    }
    else if (oldr->fip <= 0. && oldr->next==NULL) {
      oldr->x = x;
      oldr->ip = x;
      oldr->fip = 0.;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = NULL;
  }
  else {
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return -1;
    }
    iv_new->prev = oldl;
    iv_new->next = oldr;
    oldr->prev = iv_new;
    if (oldl) oldl->next = iv_new;
  }
  success = UNUR_SUCCESS;
  if (oldl) {
   success = _unur_tdr_ps_interval_parameter(gen, oldl);
  }
  if( iv_new ) {
    if (!oldl) {
      iv_new->ip = oldr->ip;
      iv_new->fip = oldr->fip;
    }
    success_r = _unur_tdr_ps_interval_parameter(gen, iv_new);
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  if ( oldr->next ) {
    success_r = _unur_tdr_ps_interval_parameter(gen, oldr);
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  if (success!=UNUR_SUCCESS) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
    if (oldl) memcpy(oldl, &oldl_bak, sizeof(struct unur_tdr_interval));
    memcpy(oldr, &oldr_bak, sizeof(struct unur_tdr_interval));
    oldr->prev = oldl;
    if (oldl) oldl->next = oldr;
    if (iv_new) {
      --(GEN->n_ivs); 
      free( iv_new );
    }
  return success;
  }
  if (oldl == NULL && iv_new)
    GEN->iv = iv_new;
  GEN->Atotal   = ( GEN->Atotal + (oldr->Ahat - oldr_bak.Ahat)
		   + ((oldl) ? (oldl->Ahat - oldl_bak.Ahat) : 0.)
		   + ((iv_new) ? iv_new->Ahat : 0.) );
  GEN->Asqueeze = ( GEN->Asqueeze + (oldr->Asqueeze - oldr_bak.Asqueeze)
		   + ((oldl) ? (oldl->Asqueeze - oldl_bak.Asqueeze) : 0.)
		   + ((iv_new) ? iv_new->Asqueeze : 0.) );
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_ps_debug_split_stop( gen,oldl,iv_new,oldr );
#endif
  if (GEN->Atotal <= 1.e10 * DBL_MIN) {
    _unur_error(gen->genid,UNUR_ERR_ROUNDOFF,"error below hat (almost) 0");
    return UNUR_ERR_ROUNDOFF;
  }
  return UNUR_SUCCESS;
} 
