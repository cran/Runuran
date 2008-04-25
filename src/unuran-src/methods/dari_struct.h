/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_dari_par { 
  int     squeeze;       
  int     size;          
  double  c_factor;      
};
struct unur_dari_gen { 
  double  vt;            
  double  vc;            
  double  vcr;           
  double  xsq[2];        
  double  y[2];          
  double  ys[2];         
  double  ac[2];         
  double  pm;            
  double  Hat[2];        
  double  c_factor;      
  int     m;             
  int     x[2];          
  int     s[2];          
  int     n[2];          
  int     size;          
  int     squeeze;       
  double *hp;            
  char   *hb;            
};
