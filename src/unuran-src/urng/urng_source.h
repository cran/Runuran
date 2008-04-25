/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef URNG_SOURCE_H_SEEN
#define URNG_SOURCE_H_SEEN
#ifdef UNUR_URNG_UNURAN
#define _unur_call_urng(urng)    ((urng)->sampleunif((urng)->state))
#else
#error
#error UNUR_URNG changed!
#error
#error Define _unur_call_urng(urng) and _unur_call_reset(urng) in
#error file 'urng_source'
#error
#endif  
#endif  
