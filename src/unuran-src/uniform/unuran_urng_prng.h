/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#  define __BEGIN_DECLS extern "C" {
#  define __END_DECLS }
#else
#  define __BEGIN_DECLS 
#  define __END_DECLS 
#endif
__BEGIN_DECLS
#ifndef UNURAN_URNG_PRNG_H_SEEN
#define UNURAN_URNG_PRNG_H_SEEN
#include <prng.h>
UNUR_URNG *unur_urng_prng_new( const char *prngstr );
UNUR_URNG *unur_urng_prngptr_new( struct prng *urng );
#endif  
__END_DECLS
