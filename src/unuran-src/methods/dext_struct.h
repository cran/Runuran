/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_dext_par {   
  int (*init)  (UNUR_GEN *gen);      
  int (*sample)(UNUR_GEN *gen);      
};
struct unur_dext_gen { 
  int (*init)  (UNUR_GEN *gen);      
  int (*sample)(UNUR_GEN *gen);      
  void *param;                       
  size_t size_param;                 
};
