/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* define macros for GCC attributes                                          */

#ifdef __GNUC__
#  define ATTRIBUTE__UNUSED        __attribute__ ((unused))
#else
#  define ATTRIBUTE__UNUSED
#endif

/*****************************************************************************/
/* Special packing functions                                                 */

void _Runuran_pack_pinv (struct unur_gen *gen, SEXP sexp_unur);
/*---------------------------------------------------------------------------*/
/* Pack Runuran generator object for method PINV into R list                 */
/*---------------------------------------------------------------------------*/

SEXP _Runuran_sample_pinv (SEXP sexp_data, int n);
/*---------------------------------------------------------------------------*/
/* Sample from generator object: use R data list (packed object)             */
/*---------------------------------------------------------------------------*/

SEXP _Runuran_quantile_pinv (SEXP sexp_data, SEXP sexp_U, SEXP sexp_unur);
/*---------------------------------------------------------------------------*/
/* Evaluate approximate quantile function:  use R data list (packed object)  */
/*---------------------------------------------------------------------------*/
