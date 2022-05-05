/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_pinv_par { 
  int order;               
  int smooth;              
  double u_resolution;     
  double bleft;            
  double bright;           
  int sleft;               
  int sright;              
  int max_ivs;             
};
struct unur_pinv_interval {
  double *ui;  
  double *zi;  
  double xi;   
  double cdfi; 
#ifdef UNUR_COOKIES
  unsigned cookie;         
#endif
};
struct unur_pinv_gen { 
  int order;               
  int    *guide;            
  int     guide_size;      
  double  Umax;            
  double  u_resolution;    
  int smooth;              
  double  bleft;           
  double  bright;          
  struct unur_pinv_interval *iv; 
  int n_ivs;               
  int max_ivs;             
  double  bleft_par;       
  double  bright_par;      
  double  dleft;           
  double  dright;          
  int sleft;               
  int sright;              
  double area;              
  struct unur_lobatto_table *aCDF; 
};
