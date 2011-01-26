/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include "tdr_gw_sample.ch"
#include "tdr_ia_sample.ch"
#include "tdr_ps_sample.ch"
double
unur_tdr_eval_invcdfhat( const struct unur_gen *gen, double u,
			 double *hx, double *fx, double *sqx )
{ 
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_TDR ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 
  if ( u<0. || u>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }
  if (u<=0.) return DISTR.domain[0];
  if (u>=1.) return DISTR.domain[1];
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    
    return _unur_tdr_gw_eval_invcdfhat(gen,u,hx,fx,sqx,NULL,NULL);
  case TDR_VARIANT_IA:    
  case TDR_VARIANT_PS:    
    return _unur_tdr_ps_eval_invcdfhat(gen,u,hx,fx,sqx,NULL);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }
} 
