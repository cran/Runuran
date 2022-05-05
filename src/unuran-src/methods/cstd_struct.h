/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_cstd_par { 
  int dummy;              
};
struct unur_cstd_gen { 
  double *gen_param;      
  int     n_gen_param;    
  int    flag;            
  double  Umin;           
  double  Umax;           
  int  is_inversion;           
  const char *sample_routine_name; 
};
