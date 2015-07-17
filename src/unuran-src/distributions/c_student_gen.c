/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include "unur_distributions_source.h"
inline static int student_trouo_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define nu (DISTR.params[0])    
int 
_unur_stdgen_student_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_student_tpol );
    return UNUR_SUCCESS;
  case 2:  
    if (par!=NULL && par->distr->data.cont.params[0] < 1.) {   
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return UNUR_ERR_GEN_CONDITION;
    }
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_student_trouo );
    return student_trouo_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
double
_unur_stdgen_sample_student_tpol( struct unur_gen *gen )
{
  double u,v,w;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  do {
    u = 2. * uniform() - 1.;
    v = 2. * uniform() - 1.;
    w = u * u + v * v;
  } while (w > 1.);
  return(u * sqrt( nu * ( exp(- 2. / nu * log(w)) - 1.) / w));
} 
#define GEN_N_PARAMS (6)
#define c       (GEN->gen_param[0])
#define e       (GEN->gen_param[1])
#define p       (GEN->gen_param[2])
#define q       (GEN->gen_param[3])
#define r       (GEN->gen_param[4])
#define vm      (GEN->gen_param[5])
inline static int
student_trouo_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (nu < 1.) {
    _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
    return UNUR_ERR_GEN_CONDITION;
  }
  r = 1. / nu;
  p = 1. / (1. + r);
  q = -0.25 * (nu + 1.);
  c = 4. * pow(p, q);
  e = 16. / c;
  vm = (nu>1.0) ? sqrt(p+p) * pow( (1.-r)*p, 0.25*(nu-1.) ) : 1.;
  return UNUR_SUCCESS;
} 
double
_unur_stdgen_sample_student_trouo( struct unur_gen *gen )
{
  double tru,u,v;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    u = uniform();
    v = uniform();
    v = vm * (v + v - 1.);
    tru = v / u;
    if ( c * u <= 5. - tru * tru) 
      break;
    if (nu >= 3.) 
      if (u * (tru * tru + 3.) >= e) 
	continue;        
    if ( u <= pow(1. + tru * tru * r, q))
      break;
  }
  return tru;
} 
#undef GEN_N_PARAMS
#undef c
#undef e
#undef p
#undef q
#undef r
#undef vm
