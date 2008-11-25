/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_pinv_sample( struct unur_gen *gen )
{ 
  double U,X;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);
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
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);
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
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY;
  }
  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);
  if ( u<0. || u>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }
  if (u<=0.) return DISTR.domain[0];
  if (u>=1.) return DISTR.domain[1];
  x = _unur_pinv_eval_approxinvcdf(gen,u);
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];
  return x;
} 
int
unur_pinv_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE )
{ 
  double score;
  _unur_check_NULL(GENTYPE, gen, UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);
  score = unur_test_inverror(gen, max_error, MAE, 1.e-20, samplesize, 
			     FALSE, FALSE, FALSE, NULL);
  return UNUR_SUCCESS;
} 
