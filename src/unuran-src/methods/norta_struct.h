/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_norta_par { 
  int dummy;
};
struct unur_norta_gen { 
  int    dim;                          
  double *copula;                      
  struct unur_distr *normaldistr;      
  struct unur_gen **marginalgen_list;  
  UNUR_URNG *marginal_urng;            
  double urng_U[1];                    
};
