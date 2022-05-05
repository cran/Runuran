/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_ERRNO_H_SEEN
#define UNUR_ERRNO_H_SEEN
enum { 
  UNUR_SUCCESS            = 0x00,                               
  UNUR_FAILURE            = 0x01,     
  UNUR_ERR_DISTR_SET      = 0x11,     
  UNUR_ERR_DISTR_GET      = 0x12,     
  UNUR_ERR_DISTR_NPARAMS  = 0x13,     
  UNUR_ERR_DISTR_DOMAIN   = 0x14,     
  UNUR_ERR_DISTR_GEN      = 0x15,     
  UNUR_ERR_DISTR_REQUIRED = 0x16,     
  UNUR_ERR_DISTR_UNKNOWN  = 0x17,     
  UNUR_ERR_DISTR_INVALID  = 0x18,     
  UNUR_ERR_DISTR_DATA     = 0x19,     
  UNUR_ERR_DISTR_PROP     = 0x20,     
  UNUR_ERR_PAR_SET        = 0x21,     
  UNUR_ERR_PAR_VARIANT    = 0x22,     
  UNUR_ERR_PAR_INVALID    = 0x23,     
  UNUR_ERR_GEN            = 0x31,     
  UNUR_ERR_GEN_DATA       = 0x32,     
  UNUR_ERR_GEN_CONDITION  = 0x33,     
  UNUR_ERR_GEN_INVALID    = 0x34,     
  UNUR_ERR_GEN_SAMPLING   = 0x35,     
  UNUR_ERR_NO_REINIT      = 0x36,     
  UNUR_ERR_NO_QUANTILE    = 0x37,     
  UNUR_ERR_URNG           = 0x41,     
  UNUR_ERR_URNG_MISS      = 0x42,     
  UNUR_ERR_STR            = 0x51,     
  UNUR_ERR_STR_UNKNOWN    = 0x52,     
  UNUR_ERR_STR_SYNTAX     = 0x53,     
  UNUR_ERR_STR_INVALID    = 0x54,     
  UNUR_ERR_FSTR_SYNTAX    = 0x55,     
  UNUR_ERR_FSTR_DERIV     = 0x56,     
  UNUR_ERR_DOMAIN         = 0x61,     
  UNUR_ERR_ROUNDOFF       = 0x62,     
  UNUR_ERR_MALLOC         = 0x63,     
  UNUR_ERR_NULL           = 0x64,      
  UNUR_ERR_COOKIE         = 0x65,     
  UNUR_ERR_GENERIC        = 0x66,     
  UNUR_ERR_SILENT         = 0x67,     
  UNUR_ERR_INF            = 0x68,     
  UNUR_ERR_NAN            = 0x69,     
  UNUR_ERR_COMPILE        = 0xa0,     
  UNUR_ERR_SHOULD_NOT_HAPPEN = 0xf0  
};
#endif  
