/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_tabl_par { 
  const double *slopes; 
  int     n_slopes;     
  double  bleft;        
  double  bright;       
  int     max_ivs;      
  double  max_ratio;    
  const double *cpoints; 
  int     n_cpoints;    
  int     n_stp;        
  double  area_fract;   
  double  darsfactor;   
  double  guide_factor; 
};
struct unur_tabl_interval {
  double  xmax;         
  double  fmax;         
  double  xmin;         
  double  fmin;         
  double  Ahat;         
  double  Asqueeze;     
  double  Acum;         
  struct unur_tabl_interval *next;  
#ifdef UNUR_COOKIES
  unsigned cookie;      
#endif
};
struct unur_tabl_gen { 
  double  Atotal;               
  double  Asqueeze;             
  double  bleft;                
  double  bright;               
  struct unur_tabl_interval **guide; 
  int     guide_size;           
  double  guide_factor;         
  double  Umin, Umax;           
  struct unur_tabl_interval *iv;     
  int     n_ivs;                
  int     max_ivs;              
  double  max_ratio;            
  double  darsfactor;           
#ifdef UNUR_ENABLE_INFO
  int     max_ivs_info;         
#endif
};
