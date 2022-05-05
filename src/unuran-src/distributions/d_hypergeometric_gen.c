/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>   
#include <methods/dstd_struct.h>
#include "unur_distributions_source.h"
#define PAR       ((struct unur_dstd_par*)par->datap) 
#define GEN       ((struct unur_dstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define uniform()  _unur_call_urng(gen->urng) 
#define par_N  (DISTR.params[0])
#define par_M  (DISTR.params[1])
#define par_n  (DISTR.params[2])
inline static int hypergeometric_hruec_init( struct unur_gen *gen );
static int _unur_stdgen_sample_hypergeometric_hruec( struct unur_gen *gen );
inline static int h_util(int N_, int M_, int n_, int k_);
int 
_unur_stdgen_hypergeometric_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
     _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_hypergeometric_hruec );
     return hypergeometric_hruec_init( gen );
  default: 
    return UNUR_FAILURE;
  }
} 
#define flogfak(k) (_unur_SF_ln_factorial(k))
#define delta(k) (flogfak(k)+flogfak(Mc-k)+flogfak(nc-k)+flogfak(NMn+k))
#define GEN_N_IPARAMS (9)
#define GEN_N_PARAMS  (8)
#define N       (GEN->gen_iparam[0])
#define M       (GEN->gen_iparam[1])
#define n       (GEN->gen_iparam[2])
#define b       (GEN->gen_iparam[3])
#define m       (GEN->gen_iparam[4])
#define NMn     (GEN->gen_iparam[5])
#define Mc      (GEN->gen_iparam[6])
#define nc      (GEN->gen_iparam[7])
#define N_half  (GEN->gen_iparam[8])
#define NMnp    (GEN->gen_param[0])
#define Np      (GEN->gen_param[1])
#define Mp      (GEN->gen_param[2])
#define np      (GEN->gen_param[3])
#define g       (GEN->gen_param[4])
#define a       (GEN->gen_param[5])
#define h       (GEN->gen_param[6])
#define p0      (GEN->gen_param[7])
int
hypergeometric_hruec_init( struct unur_gen *gen )
{
  int k1,bh;
  double x,p,q,c,my;
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
  N = (int) par_N;
  M = (int) par_M;
  n = (int) par_n;
  N_half = N/2;                      
  Mc = (M<=N_half) ? M : N-M;        
  nc = (n<=N_half) ? n : N-n;        
  Np = (double) N;
  Mp = (double) Mc;
  np = (double) nc;
  NMn = N - Mc - nc;
  NMnp = Np - Mp - np;
  p = Mp / Np;
  q = 1.0 - p;
  my = np * p;
  bh = _unur_min(nc,Mc);
  m = (int) ((np+1.0)*(Mp+1.0)/(Np+2.0));       
  if (m < 5) {
    c = my + 10.0*sqrt(my*q*(1.0-np/Np));
    b = _unur_min(bh,(int)c);                   
    p0 = exp(flogfak(N-Mc)+flogfak(N-nc)-flogfak(NMn)-flogfak(N));
    h = a = g = 0.;
  }
  else {
    a = my+0.5;
    c = sqrt(2.0*a*q*(1.0-np/Np));
    b = _unur_min(bh,(int)(a+7.0*c));           
    g = delta(m);
    k1 = (int)(a-c);
    x = (a-k1-1.0)/(a-k1);
    if((np-k1)*(p-(double)k1/Np)*x*x > (k1+1)*(q-(np-k1-1.0)/Np))
      k1++;
    h = (a-k1)*exp(0.5*(g-delta(k1))+M_LN2);    
    p0 = 0.;
  }
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_hypergeometric_hruec( struct unur_gen *gen )
{
  int k,i;
  double x,u,f,lf;
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  if (m<5) {                                     
    double pk;
    k = 0;
    pk = p0;
    u = uniform();
    while (u>pk) {
      ++k;
      if (k>b) {
	u = uniform();
	k = 0;
	pk = p0;
      }
      else {
	u -= pk;
	pk *= ((Mp-k+1.0)*(np-k+1.0)) / ((double)k*(NMnp+k));
      }
    }
    return (h_util(N,M,n,k));
  }
  for (;;) {                                    
    do {
      u = uniform();
      x = a + h*(uniform()-0.5) / u;
    } while (x < 0 || ((k=(int)x) > b));        
    if (m <= 20 || abs(m-k) <= 15) {           
      f = 1.0;
      if (m<k) {
	for (i=m+1;i<=k;i++)
	  f *= ((Mp-i+1.0)*(np-i+1.0)) / ((double)i*(NMnp+i));
	if (u*u <= f) break;                    
      }
      else {
	for (i=k+1;i<=m;i++)
	  f *= ((Mp-i+1.0)*(np-i+1.0)) / ((double)i*(NMnp+i));
	if (u*u*f <= 1.0) break;                
      }
    }
    else {
      lf = g - delta(k);                        
      if ( u * (4.0 - u) - 3.0 <= lf) break;    
      if (u*(u-lf) <= 1.0)                      
	if (2.0*log(u) <= lf) break;            
    }
  }
  return (h_util(N,M,n,k));
} 
#undef GEN_N_IPARAMS
#undef GEN_N_PARAMS
#undef N
#undef M
#undef n
#undef b    
#undef m   
#undef NMn 
#undef Mc
#undef nc
#undef N_half
#undef NMnp
#undef Np  
#undef Mp  
#undef np  
#undef g   
#undef a   
#undef h   
#undef p0  
#undef delta
#undef flogfak
#undef par_N
#undef par_M
#undef par_n
int
h_util(int N_, int M_, int n_, int k)
{
 int N_half;
 N_half = N_/2;
 if (n_ <= N_half)
	 {
	  if (M_ <= N_half) return(k);   
	  else return(n_ - k);           
	 }                            
 else
	 {
	  if (M_ <= N_half) return(M_ - k); 
	  else return(M_ - N_ + n_ + k);       
	 }                            
} 
