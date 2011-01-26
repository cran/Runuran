/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_SOURCE_H_SEEN
#define UNUR_SOURCE_H_SEEN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#else
#  error "config.h" required
#endif
#include <unuran_config.h>
#ifdef __GNUC__
#  define ATTRIBUTE__FORMAT(a,b)   __attribute__ (( __format__ (printf, (a), (b)) ))
#  define ATTRIBUTE__UNUSED        __attribute__ ((unused))
#  define ATTRIBUTE__MALLOC        __attribute__ ((malloc))
#else
#  define ATTRIBUTE__FORMAT(a,b)
#  define ATTRIBUTE__UNUSED
#  define ATTRIBUTE__MALLOC
#endif
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif
#include <unur_typedefs.h>
#include <unur_struct.h>
#if !HAVE_DECL_LOG1P
#  include <specfunct/unur_specfunct_source.h> 
#endif
#include <urng/urng_source.h>
#include <unur_cookies.h>
#include <utils/debug.h>
#include <utils/debug_source.h>
#include <utils/error.h>
#include <utils/error_source.h>
#include <utils/stream.h>
#include <utils/stream_source.h>
#include <utils/unur_errno.h>
#include <utils/unur_fp_source.h>
#include <utils/unur_fp_const_source.h>
#include <utils/umath.h>
#include <utils/umath_source.h>
#include <utils/unur_math_source.h>
#include <utils/vector_source.h>
#include <utils/string_source.h>
#include <utils/umalloc_source.h>
#include <utils/slist.h>
#endif  
