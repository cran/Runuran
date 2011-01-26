/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_urng {
  double (*sampleunif)(void *state);  
  void *state;                        
  unsigned int (*samplearray)(void *state, double *X, int dim);
  void (*sync)(void *state);          
  unsigned long seed;                 
  void (*setseed)(void *state, unsigned long seed);  
  void (*reset)(void *state);         
  void (*nextsub)(void *state);       
  void (*resetsub)(void *state);      
  void (*anti)(void *state, int a);   
  void (*delete)(void *state);        
#ifdef UNUR_COOKIES
  unsigned cookie;            
#endif
};
