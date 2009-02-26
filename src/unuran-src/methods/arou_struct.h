/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_arou_par { 
  double  guide_factor;         
  double  bound_for_adding;     
  double  max_ratio;            
  int     n_starting_cpoints;   
  const double *starting_cpoints;  
  int     max_segs;             
  double  darsfactor;           
};
struct unur_arou_segment {
  double Acum;                  
  double Ain;                   
  double Aout;                  
  double ltp[2];                
  double dltp[3];               
  double mid[2];                
  double *rtp;                  
  double *drtp;                 
  struct unur_arou_segment *next; 
#ifdef UNUR_COOKIES
  unsigned cookie;              
#endif
};
struct unur_arou_gen { 
  double  Atotal;               
  double  Asqueeze;             
  double  max_ratio;            
  struct unur_arou_segment **guide;  
  int     guide_size;           
  double  guide_factor;         
  struct unur_arou_segment *seg;     
  int     n_segs;               
  int     max_segs;             
  double  darsfactor;           
  double  center;               
#ifdef UNUR_ENABLE_INFO
  int     max_segs_info;        
#endif
};
