/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

typedef double UNUR_SAMPLING_ROUTINE_CONT(struct unur_gen *gen);
typedef int UNUR_SAMPLING_ROUTINE_DISCR(struct unur_gen *gen);
typedef int UNUR_SAMPLING_ROUTINE_CVEC(struct unur_gen *gen, double *vec);
struct unur_par {
  void *datap;                
  size_t s_datap;             
  struct unur_gen* (*init)(struct unur_par *par);
  unsigned method;            
  unsigned variant;           
  unsigned set;               
  UNUR_URNG *urng;            
  UNUR_URNG *urng_aux;        
  const struct unur_distr *distr;  
  int distr_is_privatecopy;   
  unsigned debug;             
#ifdef UNUR_COOKIES
  unsigned cookie;            
#endif
};
struct unur_gen { 
  void *datap;                
  union {
    UNUR_SAMPLING_ROUTINE_CONT  *cont;
    UNUR_SAMPLING_ROUTINE_DISCR *discr;
    UNUR_SAMPLING_ROUTINE_CVEC  *cvec;
    UNUR_SAMPLING_ROUTINE_CVEC  *matr;
  } sample;                   
  UNUR_URNG *urng;            
  UNUR_URNG *urng_aux;        
  struct unur_distr *distr;   
  int distr_is_privatecopy;   
  unsigned method;            
  unsigned variant;           
  unsigned set;               
  unsigned status;            
  char *genid;                
  struct unur_gen *gen_aux;   
  struct unur_gen **gen_aux_list; 
  int n_gen_aux_list;         
  size_t s_datap;             
  unsigned debug;             
  void (*destroy)(struct unur_gen *gen);  
  struct unur_gen* (*clone)(const struct unur_gen *gen ); 
  int (*reinit)(struct unur_gen *gen);  
#ifdef UNUR_ENABLE_INFO
  struct unur_string *infostr; 
  void (*info)(struct unur_gen *gen, int help); 
#endif
#ifdef UNUR_COOKIES
  unsigned cookie;            
#endif
};
