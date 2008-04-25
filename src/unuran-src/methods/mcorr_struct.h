/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_mcorr_par { 
  int    dim;            
  const double *eigenvalues;   
};
struct unur_mcorr_gen { 
  int    dim;            
  double *H;             
  double *M;             
  double *eigenvalues;   
};
