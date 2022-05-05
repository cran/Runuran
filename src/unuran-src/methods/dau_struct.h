/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_dau_par { 
  double  urn_factor;  
};
struct unur_dau_gen { 
  int     len;         
  int     urn_size;    
  double *qx;          
  int    *jx;          
  double  urn_factor;  
};
