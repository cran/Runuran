/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "urng_builtin.h"
#define SEED10  (12345L)
#define SEED11  (23456L)
#define SEED12  (34067L)
#define SEED20  (45678L)
#define SEED21  (56789L)
#define SEED22  (67890L)
static unsigned long x10 = SEED10;
static unsigned long x11 = SEED11;
static unsigned long x12 = SEED12;
static unsigned long x20 = SEED20;
static unsigned long x21 = SEED21;
static unsigned long x22 = SEED22;
static unsigned long x10_start = SEED10;
static unsigned long x11_start = SEED11;
static unsigned long x12_start = SEED12;
static unsigned long x20_start = SEED20;
static unsigned long x21_start = SEED21;
static unsigned long x22_start = SEED22;
double
unur_urng_MRG31k3p (void *dummy ATTRIBUTE__UNUSED)
{
# define m1      2147483647
# define m2      2147462579
# define norm    4.656612873077393e-10
# define mask11  511
# define mask12  16777215
# define mask20  65535
  register unsigned long yy1, yy2;  
  yy1 = ( (((x11 & mask11) << 22) + (x11 >> 9))
	  + (((x12 & mask12) << 7)  + (x12 >> 24)) );
  if (yy1 > m1) yy1 -= m1;
  yy1 += x12;
  if (yy1 > m1) yy1 -= m1;
  x12 = x11;  x11 = x10;  x10 = yy1;
  yy1 = ((x20 & mask20) << 15) + 21069 * (x20 >> 16);
  if (yy1 > m2) yy1 -= m2;
  yy2 = ((x22 & mask20) << 15) + 21069 * (x22 >> 16);
  if (yy2 > m2) yy2 -= m2;
  yy2 += x22;
  if (yy2 > m2) yy2 -= m2;
  yy2 += yy1;
  if (yy2 > m2) yy2 -= m2;
  x22 = x21;  x21 = x20;  x20 = yy2;
  if (x10 <= x20)
    return ((x10 - x20 + m1) * norm);
  else 
    return ((x10 - x20) * norm);
} 
void
unur_urng_MRG31k3p_seed (void *dummy ATTRIBUTE__UNUSED, unsigned long seed)
{
  if (seed==0) {
    _unur_error("URNG.mrg31k3p",UNUR_ERR_GENERIC,"seed = 0");
    return;
  }
  x10 = x10_start = seed; 
  x11 = x11_start = seed; 
  x12 = x12_start = seed; 
  x20 = x20_start = seed; 
  x21 = x21_start = seed; 
  x22 = x22_start = seed; 
} 
void
unur_urng_MRG31k3p_reset (void *dummy ATTRIBUTE__UNUSED)
{
  x10 = x10_start;
  x11 = x11_start;
  x12 = x12_start;
  x20 = x20_start;
  x21 = x21_start;
  x22 = x22_start;
} 
