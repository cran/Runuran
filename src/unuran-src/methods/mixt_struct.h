/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_mixt_par { 
  int n_comp;                   
  const double *prob;           
  struct unur_gen **comp;       
};
struct unur_mixt_gen { 
  int is_inversion;             
};
