/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_gibbs_par { 
  int thinning;             
  int burnin;               
  double  c_T;              
  const double *x0;         
};
struct unur_gibbs_gen {
  int dim;                  
  int thinning;             
  double  c_T;              
  double *state;            
  struct unur_distr *distr_condi; 
  int    coord;             
  double *direction;        
  int burnin;               
  double *x0;               
};
