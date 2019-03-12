/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include "unur_distributions_source.h"
inline static int beta_bc_init( struct unur_gen *gen );
inline static int beta_bb_init( struct unur_gen *gen );
inline static int beta_b00_init( struct unur_gen *gen );
inline static int beta_b01_init( struct unur_gen *gen );
inline static int beta_b1prs_init( struct unur_gen *gen );
#define PAR       ((struct unur_cstd_par*)par->datap) 
#define GEN       ((struct unur_cstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define uniform()  _unur_call_urng(gen->urng) 
#define p     (DISTR.params[0])   
#define q     (DISTR.params[1])   
#define a     (DISTR.params[2])   
#define b     (DISTR.params[3])   
int 
_unur_stdgen_beta_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    if (gen==NULL) return UNUR_SUCCESS; 
    if (p>1. && q>1.) {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_bb );
      return beta_bb_init( gen );
    }
    else {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_bc );
      return beta_bc_init( gen );
    }
  case 2:  
    if (gen==NULL) return UNUR_SUCCESS;  
    if (p > 1. && q > 1.) {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b1prs );
      return beta_b1prs_init( gen );
    }
    else if (p < 1. && q < 1.) {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b00 );
      return beta_b00_init( gen );
    }
    else if (_unur_isone(p) || _unur_isone(q)) { 
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_binv );
      return UNUR_SUCCESS;
    }
    else { 
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b01 );
      return beta_b01_init( gen );
    }
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_PARAMS (8)
#define am      (GEN->gen_param[0])
#define bm      (GEN->gen_param[1])
#define al      (GEN->gen_param[2])
#define alnam   (GEN->gen_param[3])
#define be      (GEN->gen_param[4])
#define si      (GEN->gen_param[5])
#define rk1     (GEN->gen_param[6])
#define rk2     (GEN->gen_param[7])
inline static int
beta_bc_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  am = (p > q) ? p : q;
  bm = (p < q) ? p : q;
  al = am + bm;
  alnam = al * log(al/am) - 1.386294361;
  be = 1.0 / bm;
  si = 1.0 + am - bm;
  rk1 = si * (0.013888889 + 0.041666667 * bm) / (am * be - 0.77777778);
  rk2 = 0.25 + (0.5 + 0.25 / si) * bm;
  return UNUR_SUCCESS;
} 
double 
_unur_stdgen_sample_beta_bc(  struct unur_gen *gen )
{
  double X;
  double u1,u2,v,w,y,z;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    u1 = uniform();
    u2 = uniform();
    if (u1 < 0.5) {
      y = u1 * u2;
      z = u1 * y;
      if ((0.25 * u2 - y + z) >= rk1) 
	continue;  
      v = be * log(u1 / (1.0 - u1));
      if (v > 80.0) {
	if (alnam < log(z))
	  continue;
	X = (_unur_FP_same(am,p)) ? 1.0 : 0.0;
	break;
      }
      else {
	w = am * exp(v);
	if ((al * (log(al / (bm + w)) + v) - 1.386294361) < log(z))
	  continue;  
	X = (!_unur_FP_same(am,p)) ? bm / (bm + w) : w / (bm + w);
	break;
      }
    }
    else {
      z = u1 * u1 * u2;
      if (z < 0.25) {
	v = be * log(u1 / (1.0 - u1));
	if (v > 80.0) {
	  X = (_unur_FP_same(am,p)) ? 1.0 : 0.0;
	  break;
	}
	w = am * exp(v);
	X = (!_unur_FP_same(am,p)) ? bm / (bm + w) : w / (bm + w);
	break;
      }
      else {
	if (z >= rk2)
	  continue;
	v = be * log(u1 / (1.0 - u1));
	if ( v > 80.0) {
	  if (alnam < log(z))
	    continue;
	  X = (_unur_FP_same(am,p)) ? 1.0 : 0.0;
	  break;
	}
	w = am * exp(v);
	if ((al * (log(al / (bm + w)) + v) - 1.386294361) < log(z))
	  continue;  
	X = (!_unur_FP_same(am,p))? bm / (bm + w) : w / (bm + w);
	break;
      }
    }
  }
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);
} 
#undef GEN_N_PARAMS
#undef am
#undef bm
#undef al
#undef alnam
#undef be
#undef si
#undef rk1
#undef rk2
#define GEN_N_PARAMS (5)
#define am      (GEN->gen_param[0])
#define bm      (GEN->gen_param[1])
#define al      (GEN->gen_param[2])
#define be      (GEN->gen_param[3])
#define ga      (GEN->gen_param[4])
inline static int
beta_bb_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  am = (p < q) ? p : q;
  bm = (p > q) ? p : q;
  al = am + bm;
  be = sqrt((al - 2.0)/(2.0 * p * q - al));
  ga = am + 1.0 / be;
  return UNUR_SUCCESS;
} 
double 
_unur_stdgen_sample_beta_bb(  struct unur_gen *gen )
{
  double X;
  double u1,u2,v,w,z,r,s,t;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    u1 = uniform();
    u2 = uniform();
    v = be * log(u1 / (1.0 - u1));
    w = am * exp(v);
    z = u1 * u1 * u2;
    r = ga * v - 1.386294361;
    s = am + r - w;
    if (s + 2.609437912 < 5.0 * z) {
      t = log(z);
      if (s < t)
	if (r + al * log(al/(bm + w)) < t) 
	  continue;
    }
    X = (_unur_FP_same(am,p)) ? w / (bm + w) : bm / (bm + w);
    break;
  }
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);
} 
#undef GEN_N_PARAMS
#undef am
#undef bm
#undef al
#undef be
#undef ga
#define GEN_N_PARAMS (8)
#define p_      (GEN->gen_param[0])
#define q_      (GEN->gen_param[1])
#define c       (GEN->gen_param[2])
#define t       (GEN->gen_param[3])
#define fp      (GEN->gen_param[4])
#define fq      (GEN->gen_param[5])
#define p1      (GEN->gen_param[6])
#define p2      (GEN->gen_param[7])
inline static int
beta_b00_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  p_ = p - 1.;
  q_ = q - 1.;
  c = (q * q_) / (p * p_);                            
  t = _unur_FP_same(c,1.) ? 0.5 : (1. - sqrt(c))/(1. - c);   
  fp = exp(p_ * log(t));
  fq = exp(q_ * log(1. - t));                        
  p1 = t/p;                                           
  p2 = (1. - t)/q + p1;                              
  return UNUR_SUCCESS;
} 
double 
_unur_stdgen_sample_beta_b00(  struct unur_gen *gen )
{
  double U, V, X, Z;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    U = uniform() * p2;
    if (U <= p1) {  
      Z = exp(log(U/p1) / p);
      X = t * Z;
      V = uniform() * fq;
      if (V <= 1. - q_ * X)
	break;
      if (V <= 1. + (fq - 1.) * Z) {
	if (log(V) <= q_ * log(1. - X))
	  break;
      }
    }
    else {          
      Z = exp( log( (U-p1)/(p2-p1) ) / q);
      X = 1. - (1. - t)*Z;
      V = uniform() * fp;
      if (V <= 1.0 - p_*(1. - X))
	break;
      if (V <= 1.0 + (fp - 1.) * Z) {
	if (log(V) <= p_ * log(X))  
	  break;
      }
    }
  }
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);
} 
#undef GEN_N_PARAMS
#undef p_
#undef q_
#undef c
#undef t
#undef fp
#undef fq
#undef p1
#undef p2
#define GEN_N_PARAMS (11)
#define pint    (GEN->gen_param[0])
#define qint    (GEN->gen_param[1])
#define p_      (GEN->gen_param[2])
#define q_      (GEN->gen_param[3])
#define t       (GEN->gen_param[4])
#define fp      (GEN->gen_param[5])
#define fq      (GEN->gen_param[6])
#define ml      (GEN->gen_param[7])
#define mu      (GEN->gen_param[8])
#define p1      (GEN->gen_param[9])
#define p2      (GEN->gen_param[10])
inline static int
beta_b01_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (p>q) {
    pint = q;
    qint = p;
  }
  else {
    pint = p;
    qint = q;
  }
  p_ = pint - 1.;
  q_ = qint - 1.;
  t = p_/(pint - qint);                   
  fq = exp((q_ - 1.) * log(1. - t));
  fp = pint - (pint + q_) * t;
  t -= (t - (1. - fp) * (1. - t) * fq / qint) / (1. - fp*fq);
  fp = exp(p_ * log(t));
  fq = exp(q_ * log(1. - t));                         
  if (q_ <= 1.0) {
    ml = (1. - fq) / t;                               
    mu = q_ * t;                                      
  }
  else {
    ml = q_;
    mu = 1. - fq;
  }
  p1 = t/pint;                                           
  p2 = fq * (1. - t)/qint + p1;                          
  return UNUR_SUCCESS;
} 
double 
_unur_stdgen_sample_beta_b01(  struct unur_gen *gen )
{
  double U, V, X, Z;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    U = uniform() * p2;
    if (U <= p1) {    
      Z = exp( log(U/p1) / pint);
      X = t * Z;
      V = uniform();
      if (V <= 1. - ml * X)
	break;
      if (V <= 1. - mu * Z)
	if (log(V) <= q_ * log(1. - X))
	  break;
    }
    else {             
      Z = exp( log((U-p1)/(p2-p1)) / qint);
      X = 1. - (1. - t) * Z;
      V = uniform() * fp;
      if (V <= 1. - p_ * (1. - X))
	break;
      if (V <= 1. + (fp - 1.) * Z)
	if (log(V) <= p_ * log(X))
	  break;
    }
  }
  if (p>q)
    X = 1. - X;
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);
} 
#undef GEN_N_PARAMS
#undef pint
#undef qint
#undef p_
#undef q_
#undef t 
#undef fp
#undef fq
#undef ml
#undef mu
#undef p1
#undef p2
#define GEN_N_PARAMS (22)
#define p_      (GEN->gen_param[0])
#define q_      (GEN->gen_param[1])
#define s       (GEN->gen_param[2])
#define m       (GEN->gen_param[3])
#define D       (GEN->gen_param[4])
#define Dl      (GEN->gen_param[5])
#define x1      (GEN->gen_param[6])
#define x2      (GEN->gen_param[7])
#define x4      (GEN->gen_param[8])
#define x5      (GEN->gen_param[9])
#define f1      (GEN->gen_param[10])
#define f2      (GEN->gen_param[11])
#define f4      (GEN->gen_param[12])
#define f5      (GEN->gen_param[13])
#define ll      (GEN->gen_param[14])
#define lr      (GEN->gen_param[15])
#define z2      (GEN->gen_param[16])
#define z4      (GEN->gen_param[17])
#define p1      (GEN->gen_param[18])
#define p2      (GEN->gen_param[19])
#define p3      (GEN->gen_param[20])
#define p4      (GEN->gen_param[21])
inline static int
beta_b1prs_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  p_ = p - 1.0;
  q_ = q - 1.0;
  s = p_ + q_;
  m = p_ / s;
  if (p_ > 1. || q_ > 1.)
    D = sqrt(m * (1. - m) / (s - 1.));
  if (p_ <= 1.) {
    x2 = Dl = m * 0.5;
    x1 = z2 = f1 = ll = 0.;
  }
  else {
    x2 = m - D;
    x1 = x2 - D;
    z2 = x2 * (1. - (1. - x2)/(s * D));
    if (x1 <= 0. || (s - 6.) * x2 - p_ + 3. > 0.) {
      x1 = z2;  x2 = (x1 + m) * 0.5;
      Dl = m - x2;
    }
    else {
      Dl = D;
    }
    f1 = exp( p_ * log(x1/m) + q_ * log((1. - x1)/(1. - m)) );
    ll = x1 * (1.0 - x1) / (s * (m - x1));            
  }
  f2 = exp( p_ * log(x2/m) + q_ * log((1. - x2)/(1. - m)) );
  if (q_ <= 1.) {
    D = (1. - m) * 0.5;
    x4 = 1. - D;
    x5 = z4 = 1.;
    f5 = lr = 0.;
  }
  else {
    x4 = m + D;
    x5 = x4 + D;
    z4 = x4 * (1. + (1. - x4)/(s * D));
    if (x5 >= 1. || (s - 6.) * x4 - p_ + 3. < 0.) {
      x5 = z4;
      x4 = (m + x5) * 0.5;
      D = x4 - m;
    }
    f5 = exp( p_ * log(x5/m) + q_ * log((1. - x5)/(1. - m)) );
    lr = x5 * (1. - x5) / (s * (x5 - m));            
  }
  f4 = exp( p_ * log(x4/m) + q_ * log((1. - x4)/(1. - m)) );
  p1 = f2 * (Dl + Dl);                                
  p2 = f4 * (D  + D) + p1;                            
  p3 = f1 * ll       + p2;                            
  p4 = f5 * lr       + p3;                            
  return UNUR_SUCCESS;
} 
double 
_unur_stdgen_sample_beta_b1prs(  struct unur_gen *gen )
{
  double U, V, W, X, Y;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  while (1) {
    U = uniform() * p4;
    if (U <= p1) {
      W = U/Dl - f2;
      if (W <= 0.) {
	X = m - U/f2;
	break;
      }
      if (W <= f1) {
	X = x2 - W/f1 * Dl;
	break;
      }
      U = uniform();
      V = Dl * U;
      X = x2 - V;
      Y = x2 + V;
      if (W * (x2 - z2) <= f2 * (X - z2))
	break;
      V = f2 + f2 - W;
      if (V < 1.) {
	if (V <= f2 + (1. - f2) * U) {
	  X = Y;
	  break;
	}
	if (V <= exp( p_ * log(Y/m) + q_ * log((1. - Y)/(1. - m)) ) ) {
	  X = Y;
	  break;
	}
      }
    }
    else 
      if (U <= p2) {
	U -= p1;
	W = U/D - f4;
	if (W <= 0.) {
	  X = m + U/f4;
	  break;
	}
	if (W <= f5) {
	  X = x4 + W/f5 * D;
	  break;
	}
	U = uniform();
	V = D * U;
	X = x4 + V;
	Y = x4 - V;
	if (W * (z4 - x4) <= f4 * (z4 - X))
	  break;
	V = f4 + f4 - W;
	if (V < 1.) {
	  if (V <= f4 + (1.0 - f4) * U) {
	    X = Y;
	    break;
	  }
	  if (V <= exp( p_ * log(Y/m) + q_ * log((1. - Y)/(1. - m)) ) ) {
	    X = Y;
	    break;
	  }
	}
      }
      else
	if (U <= p3) {                                    
	  U = (U - p2)/(p3 - p2);
	  Y = log(U);
	  X = x1 + ll * Y;
	  if (X <= 0.)                                    
	    continue; 
	  W = U * uniform();
	  if (W <= 1. + Y)
	    break;
	  W *= f1;
	}
	else {                                            
	  U = (U - p3)/(p4 - p3);
	  Y = log(U);
	  X = x5 - lr * Y;
	  if (X >= 1.)                                    
	    continue;
	  W = U * uniform();
	  if (W <= 1. + Y)
	    break;
	  W *= f5;
	}
    if (log(W) <= p_ * log(X/m) + q_ * log((1. - X)/(1. - m)))
      break;
  }
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);
} 
#undef GEN_N_PARAMS
#undef p_
#undef q_
#undef s 
#undef m 
#undef D 
#undef Dl
#undef x1
#undef x2
#undef x4
#undef x5
#undef f1
#undef f2
#undef f4
#undef f5
#undef ll
#undef lr
#undef z2
#undef z4
#undef p1
#undef p2
#undef p3
#undef p4
double 
_unur_stdgen_sample_beta_binv(  struct unur_gen *gen )
{
  double X;
  CHECK_NULL(gen,UNUR_INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_INFINITY);
  if (_unur_isone(p) && _unur_isone(q)) {
    X = uniform();
  }
  else if (_unur_isone(p)) {
    X = 1. - pow(1.-uniform(), 1/q);
  }
  else { 
    X = pow(uniform(), 1/p);
  }
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);
} 
