/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_pinv_sample( struct unur_gen *gen )
{ 
  double U,X;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_INFINITY);
  U = _unur_call_urng(gen->urng);
  X = _unur_pinv_eval_approxinvcdf(gen,U);
  if (X<DISTR.trunc[0]) return DISTR.trunc[0];
  if (X>DISTR.trunc[1]) return DISTR.trunc[1];
  return X;
} 
double
_unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
{
  int i;
  double x,un;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_INFINITY);
  un = u * GEN->Umax;
  i = GEN->guide[(int)(u * GEN->guide_size)];
  while (GEN->iv[i+1].cdfi < un)
    i++;
  un -= GEN->iv[i].cdfi;
  x = _unur_pinv_newton_eval(un, GEN->iv[i].ui, GEN->iv[i].zi, GEN->order);
  return (GEN->iv)[i].xi + x;
} 
double
unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
{
  double x;
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return UNUR_INFINITY;
  }
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_INFINITY);
  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  
  }
  x = _unur_pinv_eval_approxinvcdf(gen,u);
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];
  return x;
} 
double
unur_pinv_eval_approxcdf( const struct unur_gen *gen, double x )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return UNUR_INFINITY;
  }
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_INFINITY);
  if ( (gen->variant & PINV_VARIANT_PDF) && GEN->aCDF == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"'keepcdf' not set");
    return UNUR_INFINITY;
  }
  if (x <= DISTR.domain[0]) return 0.;
  if (x >= DISTR.domain[1]) return 1.;
  if (gen->variant & PINV_VARIANT_PDF) {
    return _unur_lobatto_eval_CDF(GEN->aCDF,x);
  }
  else {
    return (CDF(x));
  }
} 
int
unur_pinv_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE )
{ 
  _unur_check_NULL(GENTYPE, gen, UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);
  unur_test_u_error(gen, max_error, MAE, 1.e-20, samplesize, 
		     FALSE, FALSE, FALSE, NULL);
  return UNUR_SUCCESS;
} 
