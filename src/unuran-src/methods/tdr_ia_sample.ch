/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_tdr_ia_sample( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_tdr_interval *iv;
  int use_ia;
  double U, V, X;
  double fx, hx, Thx;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_INFINITY);
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return UNUR_INFINITY;
  } 
  urng = gen->urng;
  while (1) {
    U = _unur_call_urng(urng);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }
    U -= iv->Acum;    
    if (U >= - iv->sq * iv->Ahat) {
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    U += iv->Ahatr;
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
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; 
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
      }
      break;
    case TDR_VAR_T_POW:
      return 1.;
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;
    } 
    if (use_ia)
      return X;
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);      
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      return 1.;
    } 
    urng = gen->urng_aux;
    V = _unur_call_urng(urng);
    V = (iv->sq + (1 - iv->sq) * V) * hx;
    fx = PDF(X);
    if (V <= fx)
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_ps_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
  }
} 
double
_unur_tdr_ia_sample_check( struct unur_gen *gen )
{
  UNUR_URNG *urng;             
  struct unur_tdr_interval *iv;
  int use_ia;
  double U, V, X;
  double fx, hx, Thx, sqx;
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
    U = _unur_call_urng(urng);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }
    U -= iv->Acum;    
    if (U >= - iv->sq * iv->Ahat) {
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    U += iv->Ahatr;
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
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; 
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
      }
      break;
    case TDR_VAR_T_POW:
      return 1.;
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;
    } 
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);      
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      return 1.;
    } 
    fx = PDF(X);
    sqx = iv->sq*hx;
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
    if (use_ia)
      return X;
    urng = gen->urng_aux;
    V = _unur_call_urng(urng);
    V = (iv->sq + (1 - iv->sq) * V) * hx;
    if (V <= fx)
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_ps_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
  }
} 
