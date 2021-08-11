/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_METHODS_SOURCE_H_SEEN
#define UNUR_METHODS_SOURCE_H_SEEN
#include <methods/unur_methods.h>
#define _unur_check_par_object( par,type ) \
  if ( (par)->method != UNUR_METH_##type ) { \
    _unur_error(#type,UNUR_ERR_PAR_INVALID,""); \
    return (UNUR_ERR_PAR_INVALID); } \
  COOKIE_CHECK((par),CK_##type##_PAR,UNUR_ERR_COOKIE)
#define _unur_check_gen_object( gen,type,rval ) \
  if ( (gen)->method != UNUR_METH_##type ) { \
    _unur_error((gen)->genid,UNUR_ERR_GEN_INVALID,""); \
    return rval; } \
  COOKIE_CHECK((gen),CK_##type##_GEN,UNUR_ERR_COOKIE)
#endif  
