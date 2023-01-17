/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_itdr_par { 
  double xi;                 
  double cp, ct;             
};
struct unur_itdr_gen { 
  double bx;                 
  double Atot;                    
  double Ap, Ac, At;              
  double cp, xp;             
  double alphap, betap;      
  double by;                 
  double sy;                 
  double ct, xt;             
  double Tfxt, dTfxt;        
  double pole;               
  double bd_right;           
  double sign;               
  double xi;                 
};
