/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_vnrou_par { 
  double r;		    
  double *umin, *umax;      
  double vmax;              
};
struct unur_vnrou_gen { 
  int    dim;               
  double r;		    
  double *umin, *umax;      
  double vmax;              
  const double *center;       
};
