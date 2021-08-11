/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_tabl_rh_sample( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_tabl_interval *iv;
  double U,X,fx,V;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_INFINITY);
  urng = gen->urng;
  while(1) {
    U = GEN->Umin + _unur_call_urng(urng) * (GEN->Umax - GEN->Umin);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U)
      iv = iv->next;
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_INFINITY);
    U = (iv->xmax >= iv->xmin) ? (iv->Acum - U) : (U - iv->Acum + iv->Ahat);
    X = iv->xmax + U * (iv->xmin - iv->xmax)/iv->Ahat;
    V = _unur_call_urng(urng) * iv->fmax;  
    if (V <= iv->fmin)
      return X;
    fx = PDF(X);
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tabl_improve_hat( gen, iv, X, fx ) != UNUR_SUCCESS)
	   && (gen->variant & TABL_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
    if (V <= fx)
      return X;
    urng = gen->urng_aux;
  }
} 
double
_unur_tabl_rh_sample_check( struct unur_gen *gen )
{ 
  UNUR_URNG *urng;             
  struct unur_tabl_interval *iv;
  double U,X,fx,V;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_INFINITY);
  urng = gen->urng;
  while(1) {
    U = GEN->Umin + _unur_call_urng(urng) * (GEN->Umax - GEN->Umin);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U)
      iv = iv->next;
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_INFINITY);
    U = (iv->xmax >= iv->xmin) ? (iv->Acum - U) : (U - iv->Acum + iv->Ahat);
    X = iv->xmax + U * (iv->xmin - iv->xmax)/iv->Ahat;
    V = _unur_call_urng(urng) * iv->fmax;  
    fx = PDF(X);
    if (_unur_FP_greater(fx,iv->fmax))
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
    if (_unur_FP_less(fx,iv->fmin))
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");
    if (V <= iv->fmin)
      return X;
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tabl_improve_hat( gen, iv, X, fx ) != UNUR_SUCCESS)
	   && (gen->variant & TABL_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }
    if (V <= fx)
      return X;
    urng = gen->urng_aux;
  }
} 
double
_unur_tabl_ia_sample( struct unur_gen *gen )
{ 
  struct unur_tabl_interval *iv;
  double U,X,fx;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_INFINITY);
  while(1) {
    U = _unur_call_urng(gen->urng);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U)
      iv = iv->next;
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_INFINITY);
    U = (iv->xmax <= iv->xmin) ? (iv->Acum - U) : (iv->Ahat + U - iv->Acum);
    if( U < iv->Asqueeze ) {
      return( iv->xmax + (iv->Asqueeze-U) * (iv->xmin - iv->xmax)/iv->Asqueeze );
    }
    else {
      X = iv->xmax + (U-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      fx = PDF(X);
      if (GEN->n_ivs < GEN->max_ivs) {
	if ( (_unur_tabl_improve_hat( gen, iv, X, fx ) != UNUR_SUCCESS)
	     && (gen->variant & TABL_VARFLAG_PEDANTIC) )
	  return UNUR_INFINITY;
      }
      U = _unur_call_urng(gen->urng);
      if (fx >= U * (iv->fmax - iv->fmin) + iv->fmin)
	return X;
    }
  }
} 
double
_unur_tabl_ia_sample_check( struct unur_gen *gen )
{ 
  struct unur_tabl_interval *iv;
  double U,X,fx;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_INFINITY);
  while(1) {
    U = _unur_call_urng(gen->urng);
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U)
      iv = iv->next;
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_INFINITY);
    U = (iv->xmax <= iv->xmin) ? (iv->Acum - U) : (iv->Ahat + U - iv->Acum);
    if( U <= iv->Asqueeze ) {
      X = iv->xmax + (iv->Asqueeze-U) * (iv->xmin - iv->xmax)/iv->Asqueeze;
      fx = PDF(X);
      if (_unur_FP_greater(fx,iv->fmax))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
      if (_unur_FP_less(fx,iv->fmin))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");
      return X;
    }
    else {
      X = iv->xmax + (U-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      fx = PDF(X);
      if (_unur_FP_greater(fx,iv->fmax))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
      if (_unur_FP_less(fx,iv->fmin))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");
      if (GEN->n_ivs < GEN->max_ivs) {
	if ( (_unur_tabl_improve_hat( gen, iv, X, fx ) != UNUR_SUCCESS)
	     && (gen->variant & TABL_VARFLAG_PEDANTIC) )
	  return UNUR_INFINITY;
      }
      U = _unur_call_urng(gen->urng);
      if (fx >= U * (iv->fmax - iv->fmin) + iv->fmin)
	return X;
    }
  }
} 
int
_unur_tabl_improve_hat( struct unur_gen *gen, struct unur_tabl_interval *iv,
			double x, double fx)
{
  int result;
  if (! (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) ) {
    GEN->max_ivs = GEN->n_ivs;
    return UNUR_SUCCESS;
  }
  result = _unur_tabl_split_interval( gen, iv, x, fx,(gen->variant & TABL_VARMASK_SPLIT));
  if (! (result == UNUR_SUCCESS || result == UNUR_ERR_SILENT) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
    SAMPLE = _unur_sample_cont_error;
    return UNUR_ERR_GEN_CONDITION;
  }
  if ( _unur_tabl_make_guide_table(gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
    SAMPLE = _unur_sample_cont_error;
    return UNUR_ERR_GEN_CONDITION;
  }
  return UNUR_SUCCESS;
} 
