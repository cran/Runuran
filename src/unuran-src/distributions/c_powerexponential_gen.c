/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include "unur_distributions_source.h"
inline static int powerexponential_epd_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define tau    (DISTR.params[0])        
int 
_unur_stdgen_powerexponential_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    { 
      double d_tau = (par) ? par->distr->data.cont.params[0] : tau;
      if (d_tau < 1.) {
	_unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
	return UNUR_ERR_GEN_CONDITION;
      }
    }
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_powerexponential_epd );
    return powerexponential_epd_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS (2)
#define s    GEN->gen_param[0]
#define sm1  GEN->gen_param[1]
inline static int
powerexponential_epd_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  s = 1. / tau;
  sm1 = 1. - s;
  return UNUR_SUCCESS;
} 
double 
_unur_stdgen_sample_powerexponential_epd( struct unur_gen *gen )
{
  double U,u1,V,X,y;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  do {
    U = 2. * uniform() - 1.;                                  
    u1 = fabs(U);                                             
    V = uniform();                                            
    if (u1 <= sm1)
      X = u1;
    else {                       
      y = tau * (1. - u1);                                         
      X = sm1 - s * log(y);
      V *= y;
    }
  } while (log(V) > -exp(log(X)*tau));               
  if (U > 0.)
    X = -X;
  return X;
} 
#undef GEN_N_PARAMS
#undef s
#undef sm1
