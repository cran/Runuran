/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_vmt_par { 
  int dummy;
};
struct unur_vmt_gen { 
  struct unur_gen **marginalgen_list;   
  double *cholesky;         
  int    dim;               
};
