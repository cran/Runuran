/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_cext_par {   
  int    (*init)  (UNUR_GEN *gen);   
  double (*sample)(UNUR_GEN *gen);   
};
struct unur_cext_gen { 
  int    (*init)  (UNUR_GEN *gen);   
  double (*sample)(UNUR_GEN *gen);   
  void *param;               
  size_t size_param;         
};
