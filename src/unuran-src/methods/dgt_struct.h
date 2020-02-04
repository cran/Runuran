/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_dgt_par { 
  double  guide_factor; 
};
struct unur_dgt_gen { 
  double  sum;          
  double *cumpv;        
  int    *guide_table;  
  int     guide_size;   
  double  guide_factor; 
};
