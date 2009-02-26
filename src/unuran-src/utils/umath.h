/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef MATH_H_SEEN
#define MATH_H_SEEN
#include <math.h>
#ifndef INFINITY
#  ifdef HUGE_VAL
#    define INFINITY  (HUGE_VAL)
#  else
extern const double INFINITY;
#  endif
#endif
#define UNUR_INFINITY  (INFINITY)
#ifndef TRUE
#define TRUE   (1)
#endif
#ifndef FALSE
#define FALSE  (0)
#endif
#endif  
