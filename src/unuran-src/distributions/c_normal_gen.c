/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include "unur_distributions_source.h"
inline static int normal_bm_init( struct unur_gen *gen );
inline static int normal_pol_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define mu    (DISTR.params[0])   
#define sigma (DISTR.params[1])   
int 
_unur_stdgen_normal_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 1:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_bm );
    return normal_bm_init( gen );
  case 2:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_pol );
    return normal_pol_init( gen );
  case 3:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_kr );
    return UNUR_SUCCESS;
  case 0:    
  case 4:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_acr );
    return UNUR_SUCCESS;
  case 5:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_nquo );
    return UNUR_SUCCESS;
  case 6:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_quo );
    return UNUR_SUCCESS;
  case 7:    
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_leva );
    return UNUR_SUCCESS;
  case 99:   
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_sum );
    return UNUR_SUCCESS;
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS (1)
#define Xstore  GEN->gen_param[0]
#define flag    GEN->flag
int
normal_bm_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  Xstore = 0.;
  flag = 1;
  return UNUR_SUCCESS;
} 
double
_unur_stdgen_sample_normal_bm( struct unur_gen *gen )
{
  double X;
  double u,v,s;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  do {
    flag = -flag;
    if (flag > 0) { X = Xstore; break; }
    u = uniform();
    v = uniform();
    s = sqrt(-2.0 * log(u));
    Xstore = s * sin(2 * M_PI * v);
    X = s * cos(2 * M_PI * v);
  } while(0);
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
#undef GEN_N_PARAMS
#undef Xstore
#undef flag
#define GEN_N_PARAMS (1)
#define Xstore  GEN->gen_param[0]
#define flag    GEN->flag
int
normal_pol_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  Xstore = 0.;
  flag = 1;
  return UNUR_SUCCESS;
} 
double
_unur_stdgen_sample_normal_pol( struct unur_gen *gen )
{
  double X;
  double s,x,y,tmp;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  do {
    flag = -flag;
    if (flag > 0) { X = Xstore; break; }
    while(1) {
      x = 2. * uniform() - 1.;
      y = 2. * uniform() - 1.;
      s = x*x + y*y;
      if( s < 1. ) {
	tmp = sqrt( -2. * log(s) / s );
	Xstore = y * tmp;
	X = x * tmp;
	break;
      }
    }
  } while(0);
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
#undef GEN_N_PARAMS
#undef Xstore
#undef flag
double
_unur_stdgen_sample_normal_nquo( struct unur_gen *gen )
{
  double X;
  double u,v;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    u = uniform();
    if (_unur_iszero(u)) u = 1.;
    v = (uniform() - 0.5) * 0.857763885 * 2;
    X = v/u;
    if (X*X <= -4. * log(u)) 
      break;
  }
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
double
_unur_stdgen_sample_normal_quo( struct unur_gen *gen )
{
  double X;
  double r,w;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    r = uniform();
    X = (2.101083837941101 * uniform() - 1.050541918970551) / sqrt(r);
    w = X * X;
    if (4. - 4.186837275258269 * r < w) {
      if (1.5/r - 0.920558458320164 < w) 
	continue;
      if (-3.*log(r) < w )
	continue;
    }
    break;
  }
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
double
_unur_stdgen_sample_normal_leva( struct unur_gen *gen )
{
#define S    0.449871
#define T   -0.386595
#define A    0.19600
#define B    0.25472
#define RA   0.27597
#define RB   0.27846
  double X;
  double u,v,x,y,q;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    u = uniform();
    v = uniform();
    v = 1.7156 * (v - 0.5);
    x = u - S;
    y = fabs(v) - T;
    q = x * x + y * (A * y - B * x);
    X = v/u;
    if( q < RA ) break;
    if( q > RB ) continue;
    if (v*v > -4.*log(u)*u*u) continue;
    break;
  }
#undef S
#undef T  
#undef A  
#undef B  
#undef RA 
#undef RB 
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
double
_unur_stdgen_sample_normal_kr( struct unur_gen *gen )
{
#define XI 2.216035867166471
#define PIhochK 0.3989422804
  double U, V, W, X;
  double t, z;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  U = uniform();
  if (U < 0.884070402298758) {
    V = uniform();
    X = XI * (1.131131635444180 * U + V - 1.);
  }
  else if (U >= 0.973310954173898) {
    do {
      V = uniform();
      W = uniform();
      if (_unur_iszero(W)) { t=0.; continue; }
      t = XI * XI/2. - log(W);
    } while ( (V*V*t) > (XI*XI/2.) );
    X = (U < 0.986655477086949) ? pow(2*t,0.5) : -pow(2*t,0.5);
  }
  else if (U>=0.958720824790463) {
    do {
      V = uniform();
      W = uniform();
      z = V - W;
      t = XI - 0.630834801921960 * _unur_min(V,W);
    } while (_unur_max(V,W) > 0.755591531667601 &&
	     0.034240503750111 * fabs(z) > (PIhochK * exp(t*t/(-2.)) - 0.180025191068563*(XI-fabs(t))) );
    X = (z<0) ? t : -t;
  }
  else if (U>=0.911312780288703) {
    do {
      V = uniform();
      W = uniform();
      z = V - W;
      t = 0.479727404222441 + 1.105473661022070 * _unur_min(V,W);
    } while (_unur_max(V,W) > 0.872834976671790 &&
	     0.049264496373128*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
    X = (z<0) ? t : -t;
  }
  else {
    do {
      V = uniform();
      W = uniform(); 
      z = V - W;
      t = 0.479727404222441 - 0.595507138015940 * _unur_min(V,W);
      if (t<=0.) continue;
    } while (_unur_max(V,W)>0.805777924423817 &&
	     0.053377549506886*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
    X = (z<0) ? t : -t;
  }
#undef XI
#undef PIhochK 
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
double
_unur_stdgen_sample_normal_acr( struct unur_gen *gen )
{
#define c1 1.448242853
#define c2 3.307147487
#define c3 1.46754004
#define d1 1.036467755
#define d2 5.295844968
#define d3 3.631288474
#define hm 0.483941449
#define zm 0.107981933
#define hp 4.132731354
#define zp 18.52161694
#define phln 0.4515827053
#define hm1 0.516058551
#define hp1 3.132731354
#define hzm 0.375959516
#define hzmp 0.591923442
#define as 0.8853395638
#define bs 0.2452635696
#define cs 0.2770276848
#define b  0.5029324303
#define x0 0.4571828819
#define ym 0.187308492 
#define s  0.7270572718 
#define t  0.03895759111
  double X;
  double rn,x,y,z;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  do {
    y = uniform();
    if (y>hm1) {
      X = hp*y-hp1; break; }
    else if (y<zm) {  
      rn = zp*y-1;
      X = (rn>0) ? (1+rn) : (-1+rn);
      break;
    } 
    else if (y<hm) {  
      rn = uniform();
      rn = rn-1+rn;
      z = (rn>0) ? 2-rn : -2-rn;
      if ((c1-y)*(c3+fabs(z))<c2) {
	X = z; break; }
      else {  
	x = rn*rn;
	if ((y+d1)*(d3+x)<d2) {
	  X = rn; break; }
	else if (hzmp-y<exp(-(z*z+phln)/2)) {
	  X = z; break; }
	else if (y+hzm<exp(-(x+phln)/2)) {
	  X = rn; break; }
      }
    }
    while (1) {
      x = uniform();
      y = ym * uniform();
      z = x0 - s*x - y;
      if (z>0) 
	rn = 2+y/x;
      else {
	x = 1-x;
	y = ym-y;
	rn = -(2+y/x);
      }
      if ((y-as+x)*(cs+x)+bs<0) {
	X = rn; break; }
      else if (y<x+t)
	if (rn*rn<4*(b-log(x))) {
	  X = rn; break; }
    }
  } while(0);
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
double 
_unur_stdgen_sample_normal_sum( struct unur_gen *gen )
{
  double X;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  X = ( uniform() + uniform() + uniform() + uniform() + uniform() + uniform() +
	uniform() + uniform() + uniform() + uniform() + uniform() + uniform()
	- 6 );
  return ((DISTR.n_params==0) ? X : mu + sigma * X );
} 
#undef mu
#undef sigma
