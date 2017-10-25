/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_lobatto_nodes {
  double x;   
  double u;   
}; 
struct unur_lobatto_table {
  struct unur_lobatto_nodes *values; 
  int n_values;              
  int cur_iv;                
  int size;                  
  UNUR_LOBATTO_FUNCT *funct; 
  struct unur_gen *gen;      
  double tol;                
  UNUR_LOBATTO_ERROR *uerror; 
  double bleft;              
  double bright;             
  double integral;           
}; 
