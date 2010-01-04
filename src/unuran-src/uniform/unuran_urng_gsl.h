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
#ifndef UNURAN_URNG_GSL_H_SEEN
#define UNURAN_URNG_GSL_H_SEEN
#include <gsl/gsl_rng.h>
UNUR_URNG *unur_urng_gsl_new( const gsl_rng_type *urngtype );
UNUR_URNG *unur_urng_gslptr_new( gsl_rng *urng );
#include <gsl/gsl_qrng.h>
UNUR_URNG *unur_urng_gslqrng_new( const gsl_qrng_type *qrngtype, unsigned int dim );
#endif  
__END_DECLS
