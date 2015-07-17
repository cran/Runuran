/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>   
#include <methods/dstd_struct.h>
#include <methods/x_gen_source.h>
#include <distr/distr_source.h>
#include "unur_distributions_source.h"
#include "unur_distributions.h"
inline static int poisson_pdtabl_init( struct unur_gen *gen );
inline static int poisson_pdac_init( struct unur_gen *gen );
inline static int poisson_pprsc_init( struct unur_gen *gen );
#define PAR       ((struct unur_dstd_par*)par->datap) 
#define GEN       ((struct unur_dstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define uniform()  _unur_call_urng(gen->urng) 
#define theta  (DISTR.params[0])    
int 
_unur_stdgen_poisson_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
    if (gen==NULL) return UNUR_SUCCESS; 
    if (theta < 10.) {
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pdtabl );
      return poisson_pdtabl_init( gen );
    }
    else { 
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pdac );
      return poisson_pdac_init( gen );
    }
  case 2:  
    if (gen==NULL) return UNUR_SUCCESS; 
    if (theta < 10.) {
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pdtabl );
      return poisson_pdtabl_init( gen );
    }
    else { 
      _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_poisson_pprsc );
      return poisson_pprsc_init( gen );
    }
  default: 
    return UNUR_FAILURE;
  }
} 
#define GEN_N_IPARAMS (2)
#define GEN_N_PARAMS  (39)
#define m    (GEN->gen_iparam[0])
#define ll   (GEN->gen_iparam[1])
#define p0   (GEN->gen_param[0])
#define q    (GEN->gen_param[1])
#define p    (GEN->gen_param[2])
#define pp   ((GEN->gen_param)+3)  
int
poisson_pdtabl_init( struct unur_gen *gen )
{
  int i;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (GEN->gen_iparam == NULL || GEN->n_gen_iparam != GEN_N_IPARAMS) {
    GEN->n_gen_iparam = GEN_N_IPARAMS;
    GEN->gen_iparam = _unur_xrealloc(GEN->gen_iparam, GEN->n_gen_iparam * sizeof(int));
  }
  m = (theta > 1.) ? ((int) theta) : 1;
  ll = 0;
  p0 = q = p = exp(-theta);
  for (i=0; i<36; i++) pp[i]=0.;
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_poisson_pdtabl( struct unur_gen *gen )
{
  double U;
  int K,i;
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  while (1) {
    U = uniform();              
    K = 0;
    if (U <= p0) 
      return K;
    if (ll != 0) {               
      i = (U > 0.458) ? _unur_min(ll,m) : 1;
      for (K = i; K <=ll; K++)
	if (U <= pp[K])
	  return K;
      if (ll == 35) continue;
    }
    for (K = ll +1; K <= 35; K++) {
      p *= theta / (double)K;
      q += p;
      pp[K] = q;
      if (U <= q) {
	ll = K;
	return K;
      }
    }
    ll = 35;
  }
} 
#undef GEN_N_IPARAMS
#undef GEN_N_PARAMS
#undef m 
#undef ll
#undef p0
#undef q 
#undef p 
#undef pp
#define GEN_N_IPARAMS (1)
#define GEN_N_PARAMS  (10)
#define l     (GEN->gen_iparam[0])
#define s     (GEN->gen_param[0])
#define d     (GEN->gen_param[1])
#define omega (GEN->gen_param[2])
#define b1    (GEN->gen_param[3])
#define b2    (GEN->gen_param[4])
#define c     (GEN->gen_param[5])
#define c0    (GEN->gen_param[6])
#define c1    (GEN->gen_param[7])
#define c2    (GEN->gen_param[8])
#define c3    (GEN->gen_param[9])
#define NORMAL  gen->gen_aux    
int
poisson_pdac_init( struct unur_gen *gen )
{
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (GEN->gen_iparam == NULL || GEN->n_gen_iparam != GEN_N_IPARAMS) {
    GEN->n_gen_iparam = GEN_N_IPARAMS;
    GEN->gen_iparam = _unur_xrealloc(GEN->gen_iparam, GEN->n_gen_iparam * sizeof(int));
  }
  if (NORMAL==NULL) {
    struct unur_distr *distr = unur_distr_normal(NULL,0);
    struct unur_par *par = unur_cstd_new( distr );
    NORMAL = (par) ? _unur_init(par) : NULL;
    _unur_check_NULL( NULL, NORMAL, UNUR_ERR_NULL );
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
    _unur_distr_free( distr );
  }
  s = sqrt(theta);
  d = 6. * theta * theta;
  l = (int)(theta - 1.1484);
  omega = 0.3989423 / s;
  b1 = 0.416666666667e-1 / theta;
  b2 = 0.3 * b1 * b1;
  c3 = 0.1428571 * b1 * b2;
  c2 = b2 - 15.0 * c3;
  c1 = b1 - 6.0 * b2 + 45.0 * c3;
  c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
  c = 0.1069 / theta;
  return UNUR_SUCCESS;
} 
#define  a0  -0.5000000002
#define  a1   0.3333333343
#define  a2  -0.2499998565
#define  a3   0.1999997049
#define  a4  -0.1666848753
#define  a5   0.1428833286
#define  a6  -0.1241963125
#define  a7   0.1101687109
#define  a8  -0.1142650302
#define  a9   0.1055093006
int
_unur_stdgen_sample_poisson_pdac( struct unur_gen *gen )
{
  static const int fac[] = {1,1,2,6,24,120,720,5040,40320,362880};
  double t,g,theta_k;
  double gx,gy,px,py,x,xx,delta,v;
  int sign;
  double E, U;
  int K;
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  t = _unur_sample_cont(NORMAL);
  g = theta + s * t;
  if (g >= 0.) {
    K = (int) g;
    if (K >= l) 
      return K;
    U = uniform();
    theta_k = theta - K;
    if (d * U >= theta_k * theta_k * theta_k)
      return K;
    if (K < 10) {
      px = -theta;
      py = exp(K * log(theta)) / fac[K];
    }
    else {  
      delta = 0.83333333333e-1 / (double)K;
      delta = delta - 4.8 * delta*delta*delta * (1.-1./(3.5*K*K));
      v = (theta_k) / (double)K;
      if (fabs(v) > 0.25)
	px = K * log(1. + v) - theta_k - delta;
      else {
	px = K * v * v;
	px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
		  a4)*v+a3)*v+a2)*v+a1)*v+a0;
	px -= delta;
      }
      py = 0.3989422804 / sqrt((double)K);
    }
    x = (0.5 - theta_k) / s;
    xx = x * x;
    gx = -0.5 * xx;
    gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
    if (gy * (1.0 - U)  <= py * exp(px - gx))
      return K;
  }
  while (1) {
    do {
      E = - log(uniform());
      U = uniform();
      U = U + U - 1.;
      sign = (U < 0.) ? -1 : 1;
      t = 1.8 + E * sign;
    } while (t <= -0.6744);
    K = (int)(theta + s * t);
    theta_k = theta - K;
    if (K < 10) {
      px = -theta;
      py = exp(K * log(theta)) / fac[K];
    }
    else { 
      delta = 0.83333333333e-1 / (double)K;
      delta = delta - 4.8*delta*delta*delta*(1.0-1.0/(3.5*K*K));
      v = (theta_k) / (double)K;
      if (fabs(v) > 0.25)
	px = K * log(1. + v) - theta_k - delta;
      else {
	px = K * v * v;
	px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
		  a4)*v+a3)*v+a2)*v+a1)*v+a0;
	px -= delta;
      }
      py = 0.3989422804 / sqrt((double)K);
    }
    x = (0.5 - theta_k) / s;
    xx = x * x;
    gx = -0.5 * xx;
    gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
    if (c * sign * U <= py * exp(px + E) - gy * exp(gx + E)) 
      return K;
  }
} 
#undef  a0
#undef  a1
#undef  a2
#undef  a3
#undef  a4
#undef  a5
#undef  a6
#undef  a7
#undef  a8
#undef  a9
#undef GEN_N_IPARAMS
#undef GEN_N_PARAMS
#undef l
#undef s
#undef d
#undef omega
#undef b1
#undef b2
#undef c 
#undef c0
#undef c1
#undef c2
#undef c3
#undef NORMAL
inline static double f(int k, double l_nu, double c_pm)
{
  return  exp(k * l_nu - _unur_SF_ln_factorial(k) - c_pm);
}
#define GEN_N_IPARAMS (5)
#define GEN_N_PARAMS  (20)
#define m       (GEN->gen_iparam[0])
#define k2      (GEN->gen_iparam[1])
#define k4      (GEN->gen_iparam[2])
#define k1      (GEN->gen_iparam[3])
#define k5      (GEN->gen_iparam[4])
#define dl      (GEN->gen_param[0])
#define dr      (GEN->gen_param[1])
#define r1      (GEN->gen_param[2])
#define r2      (GEN->gen_param[3])
#define r4      (GEN->gen_param[4])
#define r5      (GEN->gen_param[5])
#define ll      (GEN->gen_param[6])
#define lr      (GEN->gen_param[7])
#define l_theta (GEN->gen_param[8])
#define c_pm    (GEN->gen_param[9])
#define f2      (GEN->gen_param[10])
#define f4      (GEN->gen_param[11])
#define f1      (GEN->gen_param[12])
#define f5      (GEN->gen_param[13])
#define p1      (GEN->gen_param[14])
#define p2      (GEN->gen_param[15])
#define p3      (GEN->gen_param[16])
#define p4      (GEN->gen_param[17])
#define p5      (GEN->gen_param[18])
#define p6      (GEN->gen_param[19])
int
poisson_pprsc_init( struct unur_gen *gen )
{
  double Ds;
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL || GEN->n_gen_param != GEN_N_PARAMS) {
    GEN->n_gen_param = GEN_N_PARAMS;
    GEN->gen_param = _unur_xrealloc(GEN->gen_param, GEN->n_gen_param * sizeof(double));
  }
  if (GEN->gen_iparam == NULL || GEN->n_gen_iparam != GEN_N_IPARAMS) {
    GEN->n_gen_iparam = GEN_N_IPARAMS;
    GEN->gen_iparam = _unur_xrealloc(GEN->gen_iparam, GEN->n_gen_iparam * sizeof(int));
  }
  Ds = sqrt(theta + 0.25);
  m  = (int) theta;
  k2 = (int) (theta + 0.5 - Ds);
  k4 = (int) (theta - 0.5 + Ds);
  k1 = k2 + k2 - m + 1;
  k5 = k4 + k4 - m;
  dl = (double) (k2 - k1);
  dr = (double) (k5 - k4);
  r1 = theta / (double) k1;
  r2 = theta / (double) k2;
  r4 = theta / (double)(k4 + 1);
  r5 = theta / (double)(k5 + 1);
  ll =  log(r1);                                   
  lr = -log(r5);                                   
  l_theta = log(theta);
  c_pm = m * l_theta - _unur_SF_ln_factorial(m);
  f2 = f(k2, l_theta, c_pm);
  f4 = f(k4, l_theta, c_pm);
  f1 = f(k1, l_theta, c_pm);
  f5 = f(k5, l_theta, c_pm);
  p1 = f2 * (dl + 1.);                            
  p2 = f2 * dl        + p1;                       
  p3 = f4 * (dr + 1.) + p2;                       
  p4 = f4 * dr        + p3;                       
  p5 = f1 / ll        + p4;                       
  p6 = f5 / lr        + p5;                       
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_poisson_pprsc( struct unur_gen *gen )
{
  int    Dk, X, Y;
  double U, V, W;
  while (1) {
    U = uniform() * p6;
    if (U < p2) {
      V = U - p1;
      if (V < 0.)
	return (k2 + (int)(U/f2));
      W = V / dl;
      if (W < f1)
	return (k1 + (int)(V/f1));
      Dk = (int)(dl * uniform()) + 1;
      if (W <= f2 - Dk * (f2 - f2/r2))
	return (k2 - Dk);                            
      if ((V = f2 + f2 - W) < 1.) {
	Y = k2 + Dk;
	if (V <= f2 + Dk * (1. - f2)/(dl + 1.))
	  return Y;                                  
	if (V <= f(Y, l_theta, c_pm))  
	  return Y;
      }
      X = k2 - Dk;
    }
    else if (U < p4) {
      (V = U - p3);
      if (V < 0.)
	return (k4 - (int)((U - p2)/f4));
      W = V / dr;
      if (W < f5 )
	return (k5 - (int)(V/f5));
      Dk = (int)(dr * uniform()) + 1;
      if (W <= f4 - Dk * (f4 - f4*r4))
	return(k4 + Dk);                             
      if ((V = f4 + f4 - W) < 1.0) {
	Y = k4 - Dk;
	if (V <= f4 + Dk * (1.0 - f4)/ dr)
	  return Y;                                 
	if (V <= f(Y, l_theta, c_pm))
	  return Y;       
      }
      X = k4 + Dk;
    }
    else {
      W = uniform();
      if (U < p5) {
	Dk = (int)(1. - log(W)/ll);
	X = k1 - Dk;
	if (X < 0)
	  continue;           
	W *= (U - p4) * ll;                          
	if (W <= f1 - Dk * (f1 - f1/r1))
	  return X; 
      }
      else {
	Dk = (int)(1. - log(W)/lr);
	X  = k5 + Dk;                                
	W *= (U - p5) * lr;                          
	if (W <= f5 - Dk * (f5 - f5*r5))
	  return X; 
      }
    }
    if (log(W) <= X * l_theta - _unur_SF_ln_factorial(X) - c_pm)
      return X;
  }
} 
#undef GEN_N_IPARAMS
#undef GEN_N_PARAMS
#undef m 
#undef k2
#undef k4
#undef k1
#undef k5
#undef dl
#undef dr
#undef r1
#undef r2
#undef r4
#undef r5
#undef ll
#undef lr
#undef l_theta
#undef c_pm
#undef f2
#undef f4
#undef f1
#undef f5
#undef p1
#undef p2
#undef p3
#undef p4
#undef p5
#undef p6
