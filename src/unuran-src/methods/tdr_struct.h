/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_tdr_par { 
  double  guide_factor;         
  const double *starting_cpoints; 
  int     n_starting_cpoints;   
  const double *percentiles;    
  int n_percentiles;            
  int retry_ncpoints;           
  int     max_ivs;              
  double  max_ratio;            
  double  bound_for_adding;     
  double  c_T;                             
  double  darsfactor;           
  int     darsrule;             
};
struct unur_tdr_interval {
  double  x;                    
  double  fx;                    
  double  Tfx;                   
  double  dTfx;                 
  double  sq;                   
  double  ip;                   
  double  fip;                  
  double  Acum;                 
  double  Ahat;                 
  double  Ahatr;                
  double  Asqueeze;             
  struct unur_tdr_interval *next; 
  struct unur_tdr_interval *prev; 
#ifdef UNUR_COOKIES
  unsigned cookie;              
#endif
};
struct unur_tdr_gen { 
  double  Atotal;               
  double  Asqueeze;             
  double  c_T;                             
  double  Umin, Umax;           
  struct unur_tdr_interval *iv; 
  int     n_ivs;                
  int     max_ivs;              
  double  max_ratio;            
  double  bound_for_adding;     
  struct unur_tdr_interval **guide; 
  int     guide_size;           
  double  guide_factor;         
  double  center;               
  double *starting_cpoints;     
  int     n_starting_cpoints;   
  double *percentiles;          
  int n_percentiles;            
  int retry_ncpoints;           
  double  darsfactor;           
  int     darsrule;             
#ifdef UNUR_ENABLE_INFO
  int     max_ivs_info;         
#endif
};
