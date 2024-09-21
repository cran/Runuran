/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_hitro_par { 
  double r;                  
  int thinning;              
  int burnin;                
  double adaptive_mult;      
  double vmax;               
  const double *umin, *umax; 
  const double *x0;          
};
struct unur_hitro_gen {
  int dim;                   
  int thinning;              
  double r;                  
  double *state;             
  int    coord;              
  double *direction;         
  double *vu;                
  double *vumin, *vumax;     
  double *x;                 
  const double *center;      
  double adaptive_mult;      
  int burnin;                
  double *x0;                
  double fx0;                
};
