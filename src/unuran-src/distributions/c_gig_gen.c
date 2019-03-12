/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include "unur_distributions_source.h"
inline static int gig_gigru_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define theta  (DISTR.params[0])    
#define omega  (DISTR.params[1])    
#define eta    (DISTR.params[2])    
int 
_unur_stdgen_gig_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    if (par!=NULL && par->distr->data.cont.params[0] <= 0.) {    
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return UNUR_ERR_GEN_CONDITION;
    }
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_gig_gigru );
    return gig_gigru_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS (10)
#define m       (GEN->gen_param[0])
#define linvmax (GEN->gen_param[1])
#define vminus  (GEN->gen_param[2])
#define vdiff   (GEN->gen_param[3])
#define b2      (GEN->gen_param[4])
#define hm12    (GEN->gen_param[5])
#define a       (GEN->gen_param[6])
#define d       (GEN->gen_param[7])
#define e       (GEN->gen_param[8])
#define c       (GEN->gen_param[9])
static const double drittel = 0.3333333333333333;                    
static const double pdrittel = 0.037037037037037;                    
inline static int
gig_gigru_init( struct unur_gen *gen )
{
  double r,s,t,p,q,xeta,fi,fak,yy1,yy2,max,invy1,invy2,vplus,hm1,xm,ym;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (theta <= 0) {
    _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
    return UNUR_ERR_GEN_CONDITION;
  }
  if (theta<=1. && omega<=1.) {
    e = omega * omega;
    d = theta + 1.;
    ym = (-d + sqrt(d*d + e))/omega;
    d = theta - 1.;
    xm = (d + sqrt(d*d + e))/omega;
    d = 0.5 * d;
    e = -0.25 * omega;
    r = xm + 1./xm;
    s = xm*ym;
    a = exp(-0.5*theta*log(s) + 0.5*log(xm/ym) - e*(r - ym - 1.0/ym));
    c = -d * log(xm) - e * r;
    hm12 = b2 = vdiff = vminus = linvmax = m = 0.;
  }
  else {
    hm1 = theta - 1.;
    hm12 = hm1 * 0.5;
    b2 = omega * 0.25;
    m = (hm1 + sqrt(hm1*hm1 + omega*omega))/omega;         
    max = exp(hm12 * log(m) - b2 * (m + (1./m)));          
    linvmax = log(1.0/max);
    r = (6.*m + 2.*theta*m - omega*m*m + omega)/(4.*m*m);
    s = (1. + theta - omega*m)/(2.*m*m);
    t = omega/(-4.*m*m);
    p = (3.*s - r*r) * drittel;
    q = (2.*r*r*r) * pdrittel - (r*s) * drittel + t;
    xeta = sqrt(-(p*p*p)*pdrittel);
    fi = acos(-q/(2.*xeta));
    fak = 2.*exp(log(xeta)*drittel);
    yy1 = fak * cos(fi*drittel) - r*drittel;
    yy2 = fak * cos(fi*drittel + 2.*drittel*M_PI) - r*drittel;
    invy1 = 1./yy1;
    invy2 = 1./yy2;
    vplus = exp(linvmax + log(invy1) + hm12*log(invy1 + m)
		- b2*(invy1 + m + 1./(invy1 + m)));
    vminus = -exp(linvmax + log(-invy2) + hm12 * log(invy2 + m)
		  - b2*(invy2 + m + 1./(invy2 + m)));
    vdiff = vplus - vminus;
    c = e = d = a = 0.;
  }
  return UNUR_SUCCESS;
} 
double
_unur_stdgen_sample_gig_gigru( struct unur_gen *gen )
{
  double U,V,X,Z;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  if (theta<=1. && omega<=1.) {
    do {
      U = uniform();                                        
      V = uniform();                                        
      X = a*(V/U);
    }                                         
    while (((log(U)) > (d*log(X) + e*(X + 1./X) + c)));
  }
  else {
    do {
      do {
	U = uniform();                                      
	V = vminus + uniform() * vdiff;                   
	Z = V/U;
      } while (Z < (-m));
      X = Z + m;
    }                                         
    while ((log(U) > (linvmax + hm12 * log(X) - b2 * (X + 1./X))));
  }
  return ((DISTR.n_params==2) ? X : eta * X );
} 
#undef GEN_N_PARAMS
#undef m
#undef linvmax
#undef vminus
#undef vdiff
#undef b2
#undef hm12
#undef a
#undef d
#undef e
#undef c
