/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_empk_par {
  const UNUR_GEN *kerngen;  
  UNUR_GEN *kernel;    
  double  alpha;       
  double  beta;        
  double  smoothing;   
  double  kernvar;     
};
struct unur_empk_gen {
  double *observ;      
  int     n_observ;    
  UNUR_GEN *kerngen;   
  double  smoothing;   
  double  kernvar;     
  double  bwidth;      
  double  bwidth_opt;  
  double  mean_observ; 
  double  stddev_observ; 
  double  sconst;      
  double  alpha;       
  double  beta;        
};
