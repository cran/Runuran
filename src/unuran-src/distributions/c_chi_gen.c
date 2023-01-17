/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include "unur_distributions_source.h"
inline static int chi_chru_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define nu (DISTR.params[0])    
int 
_unur_stdgen_chi_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    { 
      double d_nu = (par) ? par->distr->data.cont.params[0] : nu;
      if (d_nu < 1.) {
	_unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
	return UNUR_ERR_GEN_CONDITION;
      }
    }
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_chi_chru );
    return chi_chru_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS (4)
#define b       (GEN->gen_param[0])
#define vm      (GEN->gen_param[1])
#define vp      (GEN->gen_param[2])
#define vd      (GEN->gen_param[3])
inline static int
chi_chru_init( struct unur_gen *gen )
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
  if (_unur_isone(nu))
    return UNUR_SUCCESS;
  b = sqrt(nu - 1.);
  vm = - 0.6065306597 * (1. - 0.25 / (b * b + 1.));
  vm = (-b > vm) ? -b : vm;
  vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
  vd = vp - vm;
  return UNUR_SUCCESS;
} 
double
_unur_stdgen_sample_chi_chru( struct unur_gen *gen )
{
  double u,v,z,zz,r;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  if (_unur_isone(nu)) {
    while (1) {
      u = uniform();
      v = uniform() * 0.857763884960707;
      z = v / u;
      if (z < 0) continue;
      zz = z * z;
      r = 2.5 - zz;
      if (z < 0.)
	r = r + zz * z / (3. * z);
      if (u < r * 0.3894003915)
	break;
      if (zz > (1.036961043 / u + 1.4))
	continue;
      if (2 * log(u) < (- zz * 0.5 ))
	break;
    }
  }
  else { 
    while (1) {
      u = uniform();
      v = uniform() * vd + vm;
      z = v / u;
      if (z < -b)
	continue;
      zz = z * z;
      r = 2.5 - zz;
      if (z < 0.0)
	r = r + zz * z / (3.0 * (z + b));
      if (u < r * 0.3894003915) {
	z += b;
	break;
      }
      if (zz > (1.036961043 / u + 1.4))
	continue;
      if (2. * log(u) < (log(1.0 + z / b) * b * b - zz * 0.5 - z * b)) {
	z += b;
	break;
      }
    }
  }
  return z;
} 
#undef GEN_N_PARAMS
#undef b 
#undef vm
#undef vp
#undef vd
