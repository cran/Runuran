/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "fmax_source.h"
double
_unur_util_find_max( struct unur_funct_generic fs,  
                     double interval_min,   
		     double interval_max,   
		     double guess_max       
		   )
{
#define MAX_SRCH (100)
  int i;
  double x[3];      
  double fx[3];  
  double max;   
  double max_l; 
  double max_u; 
  double step;
  int unbound_left; 
  int unbound_right; 
  max = (_unur_FP_is_infinity(guess_max)) ? 0. : guess_max;
  if ( _unur_FP_is_minus_infinity(interval_min) &&
       _unur_FP_is_infinity(interval_max) ){
    unbound_left = 1;
    unbound_right = 1;
    x[1]  = max;
    fx[1] = fs.f(x[1], fs.params);    
    max_l = max - 100.0;
    max_u = max + 100.0;
  }
  else if ( ! _unur_FP_is_minus_infinity(interval_min) &&
              _unur_FP_is_infinity(interval_max) ){
    unbound_left = 0;
    unbound_right = 1;
    if ( max >= interval_min ){
      x[1]  = max;
      fx[1] = fs.f(x[1], fs.params);
      max_l = interval_min;
      max_u = 2 * max - interval_min;
    }
    else{
      x[1] = interval_min + 100.0;
      fx[1] = fs.f(x[1], fs.params);
      max_l = interval_min;
      max_u = x[1] + 100.0;
    }
  }
  else if ( _unur_FP_is_minus_infinity(interval_min) &&
          ! _unur_FP_is_infinity(interval_max) ){
    unbound_left = 1;
    unbound_right = 0;
    if ( max <= interval_max ){
      x[1]  = max;
      fx[1] = fs.f(x[1], fs.params);
      max_l = interval_max - 2 * max;
      max_u = interval_max;
    }
    else{
      x[1] = interval_max - 100.0;
      fx[1] = fs.f(x[1], fs.params);
      max_l = x[1] - 100.0;
      max_u = interval_max;      
    }
  }
  else {  
    unbound_left = 0;
    unbound_right = 0;
    if ( max >= interval_min && max <= interval_max ){
      x[1]  = max;
      fx[1] = fs.f(x[1], fs.params);
    }
    else{
      x[1] = interval_min/2.0 + interval_max/2.0;
      fx[1] = fs.f(x[1], fs.params);
    }
    max_l = interval_min;
    max_u = interval_max;
  }
  max = x[1];  
  step = pow(x[1]-max_l, 1.0/MAX_SRCH);
  i = 0;  
  while (i <= MAX_SRCH && _unur_FP_same(0.0, fx[1]) ){
    x[1]  = max - pow(step, (double)i);
    fx[1] = fs.f(x[1], fs.params);
    i++;
  }
  if( _unur_FP_same(0.0, fx[1]) ){
    step = pow(max_u-x[1], 1.0/MAX_SRCH);
    i = 0;
    while (i <= MAX_SRCH && _unur_FP_same(0.0, fx[1]) ){
      x[1]  = max + pow(step, (double)i);
      fx[1] = fs.f(x[1], fs.params);
      i++;
    }
  }
  if( _unur_FP_same(fx[1], 0.0) )
     return INFINITY; 
  if ( unbound_left ){
    x[2] = x[1];       fx[2] = fx[1];
    x[1] = x[2] - 1.0; fx[1] = fs.f(x[1], fs.params);
    x[0] = x[2] - 2.0; fx[0] = fs.f(x[0], fs.params);
  }
  else if ( unbound_right ){
    x[0] = x[1];       fx[0] = fx[1];
    x[1] = x[0] + 1.0; fx[1] = fs.f(x[1], fs.params);
    x[2] = x[0] + 2.0; fx[2] = fs.f(x[2], fs.params);
  }
  else{      
    x[0] = interval_min;  fx[0] = fs.f(x[0], fs.params);
    x[2] = interval_max;  fx[2] = fs.f(x[2], fs.params);
    if ( _unur_FP_same(x[1], interval_min) ||
         _unur_FP_same(x[1], interval_max) ){
      x[1]  = interval_min/2.0 + interval_max/2.0; 
      fx[1] = fs.f(x[1], fs.params);
    }
  }
  step = 1.0;
  if ( unbound_right ){
    while(fx[0] <= fx[1] && fx[1] <= fx[2]){ 
      step *= 2.0;
      x[0]  = x[1]; fx[0] = fx[1];
      x[1]  = x[2]; fx[1] = fx[2];
      x[2] += step; fx[2] = fs.f(x[2], fs.params);   
    }
  }
  step = 1.0;  
  if ( unbound_left ){
    while(fx[0] >= fx[1] && fx[1] >= fx[2]){ 
      step *= 2.0;
      x[2]  = x[1]; fx[2] = fx[1];
      x[1]  = x[0]; fx[1] = fx[0];
      x[0] -= step; fx[0] = fs.f(x[0], fs.params);
    }
  }
  max = _unur_util_brent( fs, x[0], x[2], x[1], FLT_MIN );
  if (!(_unur_FP_is_infinity( max )) ){
  }
  else {
    return INFINITY; 
  }
  return max; 
#undef MAX_SRCH
} 
double
_unur_util_brent(            
     struct unur_funct_generic fs, 
     double a,                     
     double b,                     
     double c,                     
     double tol)                   
{
#define SQRT_EPSILON  (1.e-7)           
#define MAXIT         (1000)            
#define f(x) ( (-1) * ((fs.f)(x, fs.params)) )          
  int i;
  double x,v,w;                         
  double fx;                            
  double fv;                            
  double fw;                            
  const double r = (3.-sqrt(5.0))/2;    
  CHECK_NULL(fs.f,INFINITY);
  if ( tol < 0. || b <= a || c <= a || b <= c) {
    _unur_error("CMAX",UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }
  v = c;  fv = f(v);                    
  x = v;  w = v;
  fx=fv;  fw=fv;
  for(i=0; i < MAXIT; i++) {            
    double range = b-a;                 
    double middle_range = (a+b)/2;
    double tol_act =                    
		SQRT_EPSILON*fabs(x) + tol/3;
    double new_step;                    
    if( fabs(x-middle_range) + range/2 <= 2*tol_act )
      return x;                         
    new_step = r * ( x<middle_range ? b-x : a-x );
    if( fabs(x-w) >= tol_act ) {        
      register double p;                
      register double q;                
      register double t;
      t = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*t;
      q = 2*(q-t);
      if( q>(double)0 )                 
	p = -p;                         
      else                              
	q = -q;
      if( fabs(p) < fabs(new_step*q) && 
	  p > q*(a-x+2*tol_act)      && 
	  p < q*(b-x-2*tol_act)  )      
	new_step = p/q;                 
    }
    if( fabs(new_step) < tol_act ) {    
      if( new_step > (double)0 )        
	new_step = tol_act;
      else
	new_step = -tol_act;
    }
    {                            
      register double t = x + new_step;	
      register double ft = f(t);
      if( ft <= fx ) {                  
	if( t < x )   
	  b = x;                        
	else                            
	  a = x;
	v = w;  w = x;  x = t;          
	fv=fw;  fw=fx;  fx=ft;
      }
      else {                            
	if( t < x )
	  a = t;                        
	else
	  b = t;
        if( ft <= fw || _unur_FP_same(w,x) ) {
	  v = w;  w = t;
	  fv=fw;  fw=ft;
        }
        else if( ft<=fv || _unur_FP_same(v,x) || _unur_FP_same(v,w) ) {
	  v = t;
	  fv=ft;
        }
      }
    }                   
  }                
  return INFINITY;
#undef f
#undef MAXIT
#undef SQRT_EPSILON
} 
