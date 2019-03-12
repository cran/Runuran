/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_norta_par { 
  int dummy;
};
struct unur_norta_gen { 
  int    dim;                          
  double *copula;                      
  struct unur_distr *normaldistr;      
  struct unur_gen **marginalgen_list;  
};
