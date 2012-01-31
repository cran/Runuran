/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double 
_unur_ninv_sample_newton( struct unur_gen *gen )
{
  return _unur_ninv_newton( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
}
double
_unur_ninv_sample_regula( struct unur_gen *gen )
{
  return _unur_ninv_regula( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
} 
double
_unur_ninv_sample_bisect( struct unur_gen *gen )
{
  return _unur_ninv_bisect( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
} 
double
unur_ninv_eval_approxinvcdf( const struct unur_gen *gen, double u )
{ 
  double x;
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);
  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  
  }
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    x = _unur_ninv_newton(gen,u);
    break;
  case NINV_VARFLAG_BISECT:
    x = _unur_ninv_bisect(gen,u);
    break;
  case NINV_VARFLAG_REGULA:
  default:
    x = _unur_ninv_regula(gen,u);
    break;
  }
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];
  return x;
} 
