/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

/* uncomment this line to switch on some consistency checks */
/* #define RUNURAN_DEBUG 1 */

/*---------------------------------------------------------------------------*/

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
/* UNU.RAN objects                                                           */

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

SEXP Runuran_PDF (SEXP sexp_distr, SEXP sexp_x);
/*---------------------------------------------------------------------------*/
/* Evaluate PDF or PMF for UNU.RAN distribution object. [EXPERIMENTAL]       */
/*---------------------------------------------------------------------------*/

SEXP Runuran_print (SEXP sexp_unur, SEXP sexp_help);
/*---------------------------------------------------------------------------*/
/* Print information about UNU.RAN generator object.                         */
/*---------------------------------------------------------------------------*/

SEXP Runuran_pack (SEXP sexp_unur);
/*---------------------------------------------------------------------------*/
/* Pack Runuran objects into R lists                                         */
/*---------------------------------------------------------------------------*/

SEXP Runuran_performance (SEXP sexp_unur);
/*---------------------------------------------------------------------------*/
/* Get some informations about UNU.RAN generator object in an R list.        */
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/* UNU.RAN distribution objects                                              */

SEXP Runuran_discr_init (SEXP sexp_obj, SEXP sexp_env,
			 SEXP sexp_pv, SEXP sexp_pmf,
			 SEXP sexp_mode, SEXP sexp_domain,
			 SEXP sexp_sum, SEXP sexp_name);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for discrete distribution.           */
/*---------------------------------------------------------------------------*/

SEXP Runuran_cont_init (SEXP sexp_obj, SEXP sexp_env, 
			SEXP sexp_cdf, SEXP sexp_pdf, SEXP sexp_dpdf, SEXP sexp_islog,
			SEXP sexp_mode, SEXP sexp_center, SEXP sexp_domain, 
			SEXP sexp_area, SEXP sexp_name);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for continuous distribution.         */
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
/* Internal functions (not used by .Call from R)                             */
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

void _Runuran_distr_free(SEXP sexp_distr);
/*---------------------------------------------------------------------------*/
/* Free UNU.RAN distribution object.                                         */
/*---------------------------------------------------------------------------*/

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
