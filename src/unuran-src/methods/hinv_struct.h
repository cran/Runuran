/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_hinv_par { 
  int order;               
  double u_resolution;     
  double  guide_factor;    
  double  bleft;           
  double  bright;          
  const double *stp;       
  int     n_stp;           
  int     max_ivs;         
};
#define UNUR_HINV_MAX_ORDER   (5)
struct unur_hinv_interval {
  double spline[UNUR_HINV_MAX_ORDER+1];   
  double p;                  
  double u;                
  double f;                
  double df;               
  struct unur_hinv_interval *next;  
#ifdef UNUR_COOKIES
  unsigned cookie;         
#endif
};
struct unur_hinv_gen { 
  int order;               
  int N;                   
  double *intervals;       
  int    *guide;            
  int     guide_size;      
  double  guide_factor;    
  double  Umin, Umax;      
  double  CDFmin, CDFmax;  
  double  u_resolution;    
  double  bleft;           
  double  bright;          
  struct unur_hinv_interval *iv; 
  double  tailcutoff_left;  
  double  tailcutoff_right; 
  int     max_ivs;         
  const double *stp;       
  int     n_stp;           
  double  bleft_par;       
  double  bright_par;      
};
