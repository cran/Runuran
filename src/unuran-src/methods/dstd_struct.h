/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_dstd_par { 
  int dummy;              
};
struct unur_dstd_gen { 
  double *gen_param;      
  int     n_gen_param;    
  int    *gen_iparam;     
  int     n_gen_iparam;   
  double  Umin;           
  double  Umax;           
  int  is_inversion;           
  const char *sample_routine_name; 
};
