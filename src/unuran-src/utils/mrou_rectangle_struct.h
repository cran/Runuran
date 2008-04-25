/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct MROU_RECTANGLE {
  UNUR_DISTR *distr;        
  int    dim;               
  double r;	            
  int bounding_rectangle;   
  double *umin, *umax;      
  double vmax;              
  const double *center;     
  int aux_dim;              
  const char *genid;        
};
