/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/

/* UNU.RAN header files */
#include <unuran.h>

/* R header files */
#ifndef R_NO_REMAP
# define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>

/*****************************************************************************/
/* Create UNU.RAN object for distribution defined in pure C code.            */

typedef SEXP RUNURAN_EXT_FUNCT_INIT
( SEXP sexp_obj, SEXP sexp_params, SEXP sexp_domain,
  UNUR_FUNCT_CONT *cdf, UNUR_FUNCT_CONT *pdf, UNUR_FUNCT_CONT *dpdf, int islog,
  double *mode, double *center, char *name );

RUNURAN_EXT_FUNCT_INIT Runuran_ext_cont_init;

/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for continuous distribution.         */
/*---------------------------------------------------------------------------*/
