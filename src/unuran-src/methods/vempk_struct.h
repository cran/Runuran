/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_vempk_par {
  double  smoothing;   
};
struct unur_vempk_gen {
  double *observ;      
  int     n_observ;    
  int     dim;         
  UNUR_GEN *kerngen;   
  double  smoothing;   
  double  hopt;        
  double  hact;        
  double  corfac;      
  double *xbar;        
};
