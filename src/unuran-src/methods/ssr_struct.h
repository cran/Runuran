/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_ssr_par { 
  double  Fmode;             
  double  fm;                
  double  um;                
};
struct unur_ssr_gen { 
  double  fm;                
  double  um;                
  double  vl, vr;            
  double  xl, xr;            
  double  al, ar;            
  double  A;                 
  double  Aleft, Ain;        
  double  Fmode;             
};
