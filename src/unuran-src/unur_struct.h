/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_STRUCT_H_SEEN
#define UNUR_STRUCT_H_SEEN
struct unur_distr;    
struct unur_par;      
struct unur_gen;      
typedef double UNUR_FUNCT_GENERIC  (double  x, void *params);
typedef double UNUR_FUNCT_VGENERIC (double *x, void *params);
struct unur_funct_generic {
  UNUR_FUNCT_GENERIC *f;
  void *params;
};
struct unur_funct_vgeneric {
  UNUR_FUNCT_VGENERIC *f;
  void *params;
};
#include <utils/slist_struct.h>
#include <utils/string_struct.h>
#include <parser/functparser_struct.h>
#include <urng/urng_struct.h>
#include <distr/distr_struct.h>
#include <methods/x_gen_struct.h>
#endif  
