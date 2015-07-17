/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_tdr_ps_sample( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_tdr_interval *iv;
  double U, V;                 
  double X;                    
  double fx;                   
  double Thx;                  
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
    U -= iv->Acum - iv->Ahatr;    
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      if (_unur_iszero(iv->dTfx))
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;
    case TDR_VAR_T_SQRT:
      if (_unur_iszero(iv->dTfx))
	X = iv->x + U / iv->fx;
      else {
	X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
      }
      break;
    case TDR_VAR_T_POW:
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_INFINITY;
    } 
    V = _unur_call_urng(urng);
    if (V <= iv->sq)
      	return X;
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      V *= iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);      
      V *= 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      return UNUR_INFINITY;
    } 
    fx = PDF(X);
    if (V <= fx)
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_ps_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
    urng = gen->urng_aux;
  }
} 
double
_unur_tdr_ps_sample_check( struct unur_gen *gen )
{
  UNUR_URNG *urng;             
  struct unur_tdr_interval *iv;
  double U, V;                 
  double X;                    
  double fx, sqx, hx;          
  int squeeze_rejection = FALSE; 
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
    X = _unur_tdr_ps_eval_invcdfhat( gen, U, &hx, &fx, &sqx, &iv );
    V = _unur_call_urng(urng);
    if (V <= iv->sq)
      squeeze_rejection = TRUE;
    V *= hx;
    if (_unur_FP_less(X, DISTR.BD_LEFT) || _unur_FP_greater(X, DISTR.BD_RIGHT) ) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
#ifdef UNUR_ENABLE_LOGGING
      error = 1;
#endif
    }
    if (_unur_FP_greater(fx, hx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
#ifdef UNUR_ENABLE_LOGGING
      error = 1;
#endif
    }
    if (_unur_FP_less(fx, sqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
#ifdef UNUR_ENABLE_LOGGING
      error = 1;
#endif
    }
#ifdef UNUR_ENABLE_LOGGING
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_ps_debug_sample( gen, iv, X, fx, hx, sqx ); 
#endif
    if (squeeze_rejection)
      return X;
    if (V <= fx)
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_ps_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
    urng = gen->urng_aux;
  }
} 
double
_unur_tdr_ps_eval_invcdfhat( const struct unur_gen *gen, double U,
			     double *hx, double *fx, double *sqx,
			     struct unur_tdr_interval **ivl )
{
  struct unur_tdr_interval *iv;             
  double X;                                   
  double Thx;                               
  double t;                                 
  iv =  GEN->guide[(int) (U * GEN->guide_size)];
  U *= GEN->Atotal;
  while (iv->Acum < U) {
    iv = iv->next;
  }
  U -= iv->Acum - iv->Ahatr;
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    if (_unur_iszero(iv->dTfx))
      X = iv->x + U / iv->fx;
    else {
      t = iv->dTfx * U / iv->fx;
      if (fabs(t) > 1.e-6)
	X = iv->x + log(t + 1.) * U / (iv->fx * t);
      else if (fabs(t) > 1.e-8)
	X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
      else
	X = iv->x + U / iv->fx * (1 - t/2.);
    }
    break;
  case TDR_VAR_T_SQRT:
    if (_unur_iszero(iv->dTfx))
      X = iv->x + U / iv->fx;
    else {
      X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
    }
    break;
  case TDR_VAR_T_POW:
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_INFINITY;
  } 
  if (hx != NULL) 
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      *hx = iv->fx * exp(iv->dTfx*(X - iv->x));
      break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);      
      *hx = 1./(Thx*Thx);
      break;
    case TDR_VAR_T_POW:
    default:
      *hx = UNUR_INFINITY;
    } 
  if (fx != NULL) {
    *fx = PDF(X);
  }
  if (sqx != NULL && hx != NULL) {
    *sqx = *hx * iv->sq;
  }
  if (ivl) *ivl = iv;
  return X;
} 
int
_unur_tdr_ps_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv,
			  double x, double fx )
{
  int result;
  if (! (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) ) {
    GEN->max_ivs = GEN->n_ivs;
    return UNUR_SUCCESS;
  }
  result = _unur_tdr_ps_interval_split(gen, iv, x, fx);
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
