/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

/* uncomment this line to switch on some consistency checks */
/* #define RUNURAN_DEBUG 1 */

/*---------------------------------------------------------------------------*/

/* R header files */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* UNU.RAN header files */
#include <unuran.h>

/*****************************************************************************/
/*                                                                           */
/*   R <--> C programming interface (used in .Call)                          */
/*   (public part)                                                           */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/* Loading and unloading DLL                                                 */

void R_init_Runuran (DllInfo *info);
/*---------------------------------------------------------------------------*/
/* Initialization routine when loading the DLL.                              */
/*---------------------------------------------------------------------------*/

void R_unload_Runuran (DllInfo *info);
/*---------------------------------------------------------------------------*/
/* Clear memory before unloading the DLL.                                    */
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/* Create and handle UNU.RAN generator objects                               */

SEXP Runuran_init (SEXP sexp_obj, SEXP sexp_distr, SEXP sexp_method);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN generator object.                           */
/*---------------------------------------------------------------------------*/

SEXP Runuran_sample (SEXP sexp_unur, SEXP sexp_n);
/*---------------------------------------------------------------------------*/
/* Sample from UNU.RAN generator object.                                     */
/*---------------------------------------------------------------------------*/

SEXP Runuran_quantile (SEXP sexp_unur, SEXP sexp_U);
/*---------------------------------------------------------------------------*/
/* Quantile for distribution in UNU.RAN generator object.                    */
/*---------------------------------------------------------------------------*/

SEXP Runuran_PDF (SEXP sexp_obj, SEXP sexp_x, SEXP sexp_islog);
/*---------------------------------------------------------------------------*/
/* Evaluate PDF or PMF for UNU.RAN distribution or generator object.         */
/*---------------------------------------------------------------------------*/

SEXP Runuran_CDF (SEXP sexp_obj, SEXP sexp_x);
/*---------------------------------------------------------------------------*/
/* Evaluate CDF for UNU.RAN distribution or generator object.                */
/*---------------------------------------------------------------------------*/

SEXP Runuran_print (SEXP sexp_unur, SEXP sexp_help);
/*---------------------------------------------------------------------------*/
/* Print information about UNU.RAN generator object.                         */
/*---------------------------------------------------------------------------*/

SEXP Runuran_pack (SEXP sexp_unur);
/*---------------------------------------------------------------------------*/
/* Pack Runuran objects into R lists                                         */
/*---------------------------------------------------------------------------*/

SEXP Runuran_performance (SEXP sexp_unur, SEXP sexp_debug);
/*---------------------------------------------------------------------------*/
/* Get some informations about UNU.RAN generator object in an R list.        */
/*---------------------------------------------------------------------------*/

SEXP Runuran_verify_hat (SEXP sexp_unur, SEXP sexp_n);
/*---------------------------------------------------------------------------*/
/* Verify hat of UNU.RAN generator object that implements rejection method.  */
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/* Meta methods                                                              */

SEXP Runuran_mixt (SEXP sexp_obj, SEXP sexp_prob, SEXP sexp_comp, SEXP sexp_inversion);
/*---------------------------------------------------------------------------*/
/* Create UNU.RAN generator object for mixture of distribution.              */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/* Create and handle UNU.RAN distribution objects                            */

SEXP Runuran_cont_init (SEXP sexp_obj, SEXP sexp_env, 
			SEXP sexp_cdf, SEXP sexp_pdf, SEXP sexp_dpdf, SEXP sexp_islog,
			SEXP sexp_mode, SEXP sexp_center, SEXP sexp_domain, 
			SEXP sexp_area, SEXP sexp_name);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for continuous distribution.         */
/*---------------------------------------------------------------------------*/

SEXP Runuran_discr_init (SEXP sexp_obj, SEXP sexp_env,
			 SEXP sexp_cdf, SEXP sexp_pv, SEXP sexp_pmf,
			 SEXP sexp_mode, SEXP sexp_domain,
			 SEXP sexp_sum, SEXP sexp_name);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for discrete distribution.           */
/*---------------------------------------------------------------------------*/

SEXP Runuran_cmv_init (SEXP sexp_obj, SEXP sexp_env, 
		       SEXP sexp_dim, SEXP sexp_pdf, 
		       SEXP sexp_mode, SEXP sexp_center, 
		       SEXP sexp_ll, SEXP sexp_ur, SEXP sexp_name);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for cont. multivariate distribution. */
/*---------------------------------------------------------------------------*/


SEXP Runuran_std_cont (SEXP sexp_obj, SEXP sexp_name, SEXP sexp_params, SEXP sexp_domain);
/*---------------------------------------------------------------------------*/
/* Create UNU.RAN object for special continuous distribution.                */
/*---------------------------------------------------------------------------*/

UNUR_DISTR *_Runuran_get_std_cont( const char *name, const double *params, int n_params );
/*---------------------------------------------------------------------------*/
/* Create UNU.RAN object for distribution 'name'.                            */
/*---------------------------------------------------------------------------*/

SEXP Runuran_std_discr (SEXP sexp_obj, SEXP sexp_name, SEXP sexp_params, SEXP sexp_domain);
/*---------------------------------------------------------------------------*/
/* Create UNU.RAN object for special discrete distribution.                  */
/*---------------------------------------------------------------------------*/

UNUR_DISTR *_Runuran_get_std_discr( const char *name, const double *params, int n_params );
/*---------------------------------------------------------------------------*/
/* Create UNU.RAN object for distribution 'name'.                            */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Registering native routines */
/* Not implemented yet (it seems to be slower) */ 
/* static const R_CallMethodDef Runuran_CallEntries[] = { */
/*     {"Runuran_init", (DL_FUNC) &Runuran_init, 2}, */
/*     {"Runuran_sample", (DL_FUNC) &Runuran_sample, 2}, */
/*     {NULL, NULL, 0} */
/* }; */
/* Remark: This list has to be updated before usage. */
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/*                                                                           */
/*   Internal functions (not used by .Call from R)                           */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* define macros for GCC attributes                                          */

/* #ifdef __GNUC__ */
/* #  define ATTRIBUTE__UNUSED        __attribute__ ((unused)) */
/* #else */
/* #  define ATTRIBUTE__UNUSED */
/* #endif */

/*****************************************************************************/

#define _Runuran_fatal() \
  errorcall_return(R_NilValue,"[UNU.RAN - error] cannot create UNU.RAN distribution object")
/*---------------------------------------------------------------------------*/
/* Handle fatal error: print error message and exit.                         */
/*---------------------------------------------------------------------------*/

void _Runuran_set_error_handler(int status);
/*---------------------------------------------------------------------------*/
/* set status of error handler (on / off).                                   */
/*---------------------------------------------------------------------------*/

void _Runuran_free(SEXP sexp_gen);
/*---------------------------------------------------------------------------*/
/* Free UNU.RAN generator object.                                            */
/*---------------------------------------------------------------------------*/

SEXP _Runuran_tag(void); 
/*---------------------------------------------------------------------------*/
/* Make tag for R generator object [Contains static variable!]               */
/*---------------------------------------------------------------------------*/

void _Runuran_distr_free(SEXP sexp_distr);
/*---------------------------------------------------------------------------*/
/* Free UNU.RAN distribution object.                                         */
/*---------------------------------------------------------------------------*/

SEXP _Runuran_distr_tag(void); 
/*---------------------------------------------------------------------------*/
/* Make tag for R object [Contains static variable!]                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Check pointer to R UNU.RAN generator object.                              */
#define ALLWAYS_CHECK_UNUR_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != _Runuran_tag()) \
      error("[UNU.RAN - error] invalid UNU.RAN object");		 \
  } while (0)

#ifdef RUNURAN_DEBUG
#define CHECK_UNUR_PTR(s) ALLWAYS_CHECK_UNUR_PTR(s)
#else
#define CHECK_UNUR_PTR(s) do {} while(0)
#endif

/*---------------------------------------------------------------------------*/
/* Check pointer to R UNU.RAN distribution object.                           */
#ifdef RUNURAN_DEBUG
#define CHECK_DISTR_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != _Runuran_distr_tag()) \
      error("[UNU.RAN - error] invalid UNU.RAN distribution object");	\
  } while (0)
#else
#define CHECK_DISTR_PTR(s)  do {} while(0)
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


/*****************************************************************************/
/* Auxiliary URNG                                                            */

SEXP Runuran_use_aux_urng (SEXP sexp_unur, SEXP sexp_set);
/*---------------------------------------------------------------------------*/
/* check, set or unset auxiliary URNG for given generator object.            */
/*---------------------------------------------------------------------------*/

SEXP Runuran_set_aux_seed (SEXP sexp_seed);
/*---------------------------------------------------------------------------*/
/* set seed for auxiliary URNG.                                              */
/*---------------------------------------------------------------------------*/
