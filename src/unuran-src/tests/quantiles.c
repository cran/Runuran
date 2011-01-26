/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"
static char test_name[] = "Quantiles";
int
unur_test_quartiles( UNUR_GEN *gen, double *q0 ,double *q1, double *q2, double *q3, double *q4,
		     int samplesize, int verbosity, FILE *out )
{
  double x = 0.;
  int n;
  double  p = 0.5; 
  double  h[5];    
  int     pos[5];  
  double  dpos[5]; 
  int     i, j;    
  int     sgnd;    
  double  tmp;
  double  d;       
  _unur_check_NULL(test_name, gen, UNUR_ERR_NULL);
  if (! ( ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute quartiles for distribution");
    return UNUR_ERR_GENERIC;
  }
  if (samplesize < 10) 
    samplesize = 10;
  for (n=0; n<samplesize; n++) {
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = _unur_sample_discr(gen); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(gen); break;
    }
    if ( n==0 ){  
      h[n]    = x;
      pos[n]  = n;
      dpos[0] = 0.;
      dpos[1] = p/2.;
      dpos[2] = 4.*p;
      dpos[3] = 2. + 2.*p;
      dpos[4] = 4.;
    }                  
    else if ( n<4 ){   
      h[n]   = x;
      pos[n] = n;
    }                  
    else if ( n==4 ){
      h[n] = x;
      pos[n] = n;
      for (i=0; i<4; i++){
	for (j=0; j<4-i; j++){
	  if ( h[j] > h[j+1] ){  
	    tmp    = h[j];
	    h[j]   = h[j+1];
	    h[j+1] = tmp;
	  }
	}
      }  
    }      
    else{  
      if ( x < h[0] )
	h[0] = x;
      if ( x > h[4] )
	h[4] = x;
      for (i=1; i<4; i++){
	if ( x < h[i] ){
	  pos[i]++;
	}
      } 
      pos[4]++; 
      dpos[1] = (double) n * p/2.;
      dpos[2] = (double) n * p;
      dpos[3] = (double) n * (1.+p)/2.;
      dpos[4] = (double) n;
      for (i=1; i<4; i++){
	d = dpos[i] - pos[i];  
        if ( (d >=  1.  &&  pos[i+1]-pos[i] >  1) ||
             (d <= -1.  &&  pos[i-1]-pos[i] < -1)   ){
	  sgnd = (d < 0.) ? -1 : 1;
	  d = (double) sgnd;
          tmp = h[i] + d/(pos[i+1]-pos[i-1]) * 
	       ( (pos[i]-pos[i-1]+d) * (h[i+1]-h[i])/(pos[i+1]-pos[i]) +
                 (pos[i+1]-pos[i]-d) * (h[i]-h[i-1])/(pos[i]-pos[i-1]) );
	  if ( tmp <= h[i-1] || tmp >= h[i+1]){
	    tmp = h[i] + d*(h[i]-h[i+sgnd])/(pos[i]-pos[i+sgnd]);
	  }
	  h[i]    = tmp;      
	  pos[i] += sgnd;     
	} 
      } 
    } 
  }
  *q0 = h[0];
  *q1 = h[1];
  *q2 = h[2];
  *q3 = h[3];
  *q4 = h[4];
  if (verbosity) {
    fprintf(out,"\nQuartiles:\n");
    fprintf(out,"\tmin = \t%6.5g\n",*q0);
    fprintf(out,"\t25%% =\t%6.5g\n",*q1);
    fprintf(out,"\t50%% =\t%6.5g\n",*q2);
    fprintf(out,"\t75%% =\t%6.5g\n",*q3);
    fprintf(out,"\tmax = \t%6.5g\n",*q4);
  }
  return UNUR_SUCCESS;
} 
