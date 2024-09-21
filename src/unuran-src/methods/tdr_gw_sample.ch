/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_tdr_gw_sample( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_tdr_interval *iv, *pt;
  double U, V;                 
  double X;                    
  double fx, sqx, hx;          
  double Tsqx, Thx;            
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_INFINITY);
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return UNUR_INFINITY;
  } 
  urng = gen->urng;
  while (1) {
    U = GEN->Umin + _unur_call_urng(urng) * (GEN->Umax - GEN->Umin);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }
    U -= iv->Acum;    
    if (-U < iv->Ahatr) { 
      pt = iv->next;
    }
    else {                
      pt = iv;
      U += iv->Ahat;
    }
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      if (_unur_iszero(pt->dTfx))
	X = pt->x + U / pt->fx;
      else
	{
	  double t = pt->dTfx * U / pt->fx;
	  if (fabs(t) > 1.e-6)
	    X = pt->x + log(t + 1.) * U / (pt->fx * t);
	  else if (fabs(t) > 1.e-8)
	    X = pt->x + U / pt->fx * (1 - t/2. + t*t/3.);
	  else
	    X = pt->x + U / pt->fx * (1 - t/2.);
	}
      hx = pt->fx * exp(pt->dTfx*(X - pt->x));       
      V = _unur_call_urng(urng) * hx;  
      if (V <= iv->fx && V <= iv->next->fx)
	return X;
      sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(X - iv->x)) : 0.;     
      if (V <= sqx)
	return X;
      break;
    case TDR_VAR_T_SQRT:
      if (_unur_iszero(pt->dTfx))
	X = pt->x + U /pt->fx;
      else {
	X = pt->x + (pt->Tfx*pt->Tfx*U) / (1.-pt->Tfx*pt->dTfx*U);  
      }
      Thx = pt->Tfx + pt->dTfx * (X - pt->x);      
      hx = 1./(Thx*Thx);
      V = _unur_call_urng(urng) * hx;  
      if (V <= iv->fx && V <= iv->next->fx)
	return X;
      Tsqx = (iv->Asqueeze > 0.) ? (iv->Tfx + iv->sq * (X - iv->x)) : -UNUR_INFINITY;  
      sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
      if (V <= sqx)
	return X;
      break;
    case TDR_VAR_T_POW:
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_INFINITY;
    } 
    fx = PDF(X);
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_gw_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
    if (V <= fx)
      return X;
    urng = gen->urng_aux;
  }
} 
double
_unur_tdr_gw_sample_check( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_tdr_interval *iv, *pt;
  double U, V;                 
  double X;                    
  double fx, sqx, hx;          
#ifdef UNUR_ENABLE_LOGGING
  int error = 0;               
#endif
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_INFINITY);
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return UNUR_INFINITY;
  } 
  urng = gen->urng;
  while (1) {
    U = GEN->Umin + _unur_call_urng(urng) * (GEN->Umax - GEN->Umin);
    X = _unur_tdr_gw_eval_invcdfhat( gen, U, &hx, &fx, &sqx, &iv, &pt );
    if (X < DISTR.BD_LEFT || X > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
#ifdef UNUR_ENABLE_LOGGING
      error = 1;
#endif
    }
    if (_unur_FP_greater(fx,hx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
#ifdef UNUR_ENABLE_LOGGING
      error = 1;
#endif
    }
    if (_unur_FP_less(fx,sqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
#ifdef UNUR_ENABLE_LOGGING
      error = 1;
#endif
    }
#ifdef UNUR_ENABLE_LOGGING
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_gw_debug_sample( gen, iv, pt, X, fx, hx, sqx ); 
#endif
    V = _unur_call_urng(urng) * hx;  
    if (V <= iv->fx && V <= iv->next->fx)
      return X;
    if (V <= sqx)
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_gw_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
    if (V <= fx)
      return X;
    urng = gen->urng_aux;
  }
} 
double
_unur_tdr_gw_eval_invcdfhat( const struct unur_gen *gen, double U, 
			     double *hx, double *fx, double *sqx,
			     struct unur_tdr_interval **ivl,
			     struct unur_tdr_interval **cpt )
{ 
  struct unur_tdr_interval *iv, *pt;        
  double X;                                   
  double Tsqx, Thx;                         
  double t;                                 
  iv =  GEN->guide[(int) (U * GEN->guide_size)];
  U *= GEN->Atotal;
  while (iv->Acum < U) {
    iv = iv->next;
  }
  U -= iv->Acum;
  if (-U < iv->Ahatr) { 
    pt = iv->next;
  }
  else {                
    pt = iv;
    U += iv->Ahat;
  }
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    if (_unur_iszero(pt->dTfx))
      X = pt->x + U / pt->fx;
    else {
      t = pt->dTfx * U / pt->fx;
      if (fabs(t) > 1.e-6)
	X = pt->x + log(t + 1.) * U / (pt->fx * t);
      else if (fabs(t) > 1.e-8)
	X = pt->x + U / pt->fx * (1 - t/2. + t*t/3.);
      else
	X = pt->x + U / pt->fx * (1 - t/2.);
    }
    break;
  case TDR_VAR_T_SQRT:
    if (_unur_iszero(pt->dTfx))
      X = pt->x + U / (pt->fx);
    else {
      X = pt->x + (pt->Tfx*(pt->Tfx)*U) / (1.-(pt->Tfx)*(pt->dTfx)*U);  
    }
    break;
  case TDR_VAR_T_POW:
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    X = UNUR_INFINITY;
  } 
  if (hx != NULL) {
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      *hx = pt->fx * exp(pt->dTfx*(X - pt->x));
      break;
    case TDR_VAR_T_SQRT:
      Thx = pt->Tfx + pt->dTfx * (X - pt->x);      
      *hx = 1./(Thx*Thx);
      break;
    case TDR_VAR_T_POW:
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      *hx = UNUR_INFINITY;
    }
  }
  if (fx != NULL) {
    *fx = PDF(X);
  }
  if (sqx != NULL) {
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      *sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(X - iv->x)) : 0.;
      break;
    case TDR_VAR_T_SQRT:
      if (iv->Asqueeze > 0.) {
	Tsqx = iv->Tfx + iv->sq * (X - iv->x);
	*sqx = 1./(Tsqx*Tsqx);
      }
      else
	*sqx = 0.;
      break;
    case TDR_VAR_T_POW:
    default:  
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      *sqx = 0.;
    }
  }
  if (ivl) *ivl = iv;
  if (cpt) *cpt = pt;
  return X;
} 
int
_unur_tdr_gw_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv,
			  double x, double fx )
{
  int result;
  if (! (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) ) {
    GEN->max_ivs = GEN->n_ivs;
    return UNUR_SUCCESS;
  }
  result = _unur_tdr_gw_interval_split(gen, iv, x, fx);
  if (result!=UNUR_SUCCESS && result!=UNUR_ERR_SILENT && result!=UNUR_ERR_INF) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
    if (gen->variant & TDR_VARFLAG_PEDANTIC || result == UNUR_ERR_ROUNDOFF) {
      SAMPLE = _unur_sample_cont_error;
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  _unur_tdr_make_guide_table(gen);
  return UNUR_SUCCESS;
} 
