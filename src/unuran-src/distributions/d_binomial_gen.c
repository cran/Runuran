/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/cstd.h>   
#include <methods/dstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
inline static int binomial_bruec_init( struct unur_gen *gen );
static int _unur_stdgen_sample_binomial_bruec( struct unur_gen *gen );
#define PAR       ((struct unur_dstd_par*)par->datap) 
#define GEN       ((struct unur_dstd_gen*)gen->datap) 
#define DISTR     gen->distr->data.discr 
#define uniform()  _unur_call_urng(gen->urng) 
#define MAX_gen_params   11    
#define MAX_gen_iparams   2    
#define n  (DISTR.params[0])
#define p  (DISTR.params[1])
int 
_unur_stdgen_binomial_init( struct unur_par *par, struct unur_gen *gen )
{
  switch ((par) ? par->variant : gen->variant) {
  case 0:  
  case 1:  
     if (gen==NULL) return UNUR_SUCCESS; 
     _unur_dstd_set_sampling_routine( par,gen,_unur_stdgen_sample_binomial_bruec );
     return binomial_bruec_init( gen );
  case UNUR_STDGEN_INVERSION:   
  default: 
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
} 
#define flogfak(k) _unur_sf_ln_factorial(k)
#define b       (GEN->gen_iparam[0])
#define m       (GEN->gen_iparam[1])
#define par     (GEN->gen_param[0])
#define q1      (GEN->gen_param[1])
#define np      (GEN->gen_param[3])
#define a       (GEN->gen_param[4])
#define h       (GEN->gen_param[5])
#define g       (GEN->gen_param[6])
#define r       (GEN->gen_param[7])
#define t       (GEN->gen_param[8])
#define r1      (GEN->gen_param[9])
#define p0      (GEN->gen_param[10])
int
binomial_bruec_init( struct unur_gen *gen )
{
  int bh,k1;
  double c,x; 
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);
  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
    GEN->n_gen_iparam = MAX_gen_iparams;
    GEN->gen_iparam = _unur_xmalloc(GEN->n_gen_param * sizeof(int));
  }
  par = _unur_min(p,1.0-p);
  q1 = 1.0 - par;
  np = n*par;                                
  if (np < 5) {
    p0 = exp(n*log(q1));                     
    bh = (int)(np + 10.0*sqrt(np*q1));
    b = _unur_min(n,bh);                     
  }
  else {                                     
    m = (int)(np + par);                     
    a = np + 0.5;                            
    c = sqrt(2.0 * a * q1);
    r = par/q1;
    t = (n+1) * r;
    r1 = log(r);
    bh = (int)(a + 7.0*c);
    b = _unur_min(n,bh);                     
    g = flogfak(m) + flogfak(n-m);           
    k1 = (int)(a-c);
    x = (a-k1-1.0)/(a-k1);
    if((n-k1)*par*x*x > (k1+1)*q1)
      k1++;                                  
    h = (a-k1) * exp(.5*((k1-m)*r1+g-flogfak(k1)-flogfak(n-k1))+M_LN2);
  }
  return UNUR_SUCCESS;
} 
int
_unur_stdgen_sample_binomial_bruec( struct unur_gen *gen )
{
  int i,k;
  double u,f,x,lf;
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  if (np<5) {
    double pk;
    k = 0;
    pk = p0;
    u = uniform();
    while (u>pk) {
      ++k;
      if (k>b) {
	u = uniform();
	k = 0;
	pk=p0;
      }
      else {
	u -= pk;
	pk=(double)(((n-k+1)*par*pk)/(k*q1));
      }
    }
    return ((p>0.5) ? n-k:k);
  }
  for (;;) {
    do {
      u = uniform();
      x = a+h*(uniform()-0.5)/u;
    } while (x < 0 || ((k=(int)x) > b));    
    if ((labs(m-k)<=15) && ((k>29)||(n-k>29)) ) {
      f = 1.0;                              
      if (m<k) {
	for (i=m;i<k;) 
	  f *= t / (double)++i-r;           
	if (u*u <= f) break;                
      }
      else {
	for (i=k;i<m;)
	  f *= t / (double)++i-r;           
	if (u*u*f <= 1.0)
	  break;                            
      }
    }
    else {
      lf = (k-m)*r1+g-flogfak(k)-flogfak(n-k);       
      if ( u * (4.0 - u) - 3.0 <= lf) break;         
      if (u*(u-lf) <= 1.0)                           
	if (2.0*log(u) <= lf) break;                 
    }
  }
  return((p > 0.5) ? n-k : k);
} 
#undef b       
#undef m       
#undef b       
#undef par     
#undef q1      
#undef np      
#undef al      
#undef h       
#undef g       
#undef r       
#undef t       
#undef r1      
#undef p0      
