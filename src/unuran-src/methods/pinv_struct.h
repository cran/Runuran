/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_pinv_par { 
  int order;               
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
struct unur_pinv_CDFvalues {
  double x;   
  double u;   
}; 
struct unur_pinv_CDFtable {
  struct unur_pinv_CDFvalues *values; 
  int n_values;            
  int cur_iv;              
  int size;                
}; 
struct unur_pinv_gen { 
  int order;               
  int    *guide;            
  int     guide_size;      
  double  Umax;            
  double  u_resolution;    
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
  double logPDFconstant;   
#ifdef PINV_USE_CDFTABLE
  struct unur_pinv_CDFtable *CDFtable; 
#endif
};
