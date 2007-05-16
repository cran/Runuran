/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_tdrgw_par { 
  const double *starting_cpoints; 
  int n_starting_cpoints;         
  const double *percentiles; 
  int n_percentiles;         
  int retry_ncpoints;        
  int max_ivs;               
};
struct unur_tdrgw_interval {
  double  x;              
  double  logfx;          
  double  dlogfx;         
  double  sq;             
  double  Acum;           
  double  logAhat;        
  double  Ahatr_fract;    
  struct unur_tdrgw_interval *next; 
#ifdef DEBUG_STORE_IP 
  double  ip;             
#endif
#ifdef UNUR_COOKIES
  unsigned cookie;        
#endif
};
struct unur_tdrgw_gen { 
  double  Atotal;               
  double  logAmax;              
  struct unur_tdrgw_interval *iv; 
  int     n_ivs;                
  int     max_ivs;              
  double *starting_cpoints;     
  int     n_starting_cpoints;   
  double *percentiles;       
  int n_percentiles;         
  int retry_ncpoints;        
};
