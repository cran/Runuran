/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_ninv_par { 
  int     max_iter;          
  double  x_resolution;      
  double  u_resolution;      
  double  s[2];              
  int     table_on;          
  int     table_size;        
};
struct unur_ninv_gen { 
  int     max_iter;          
  double  x_resolution;      
  double  u_resolution;      
  double *table;             
  double *f_table;	     
  int     table_on;          
  int     table_size;        
  double  Umin, Umax;        
  double  CDFmin, CDFmax;    
  double  s[2];              
  double  CDFs[2];           
};
