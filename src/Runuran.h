/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

SEXP Runuran_init (SEXP sexp_distr, SEXP sexp_method);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN generator object.                           */
/*---------------------------------------------------------------------------*/

SEXP Runuran_sample (SEXP sexp_gen, SEXP sexp_n);
/*---------------------------------------------------------------------------*/
/* Sample from UNU.RAN generator object.                                     */
/*---------------------------------------------------------------------------*/

/* Registering native routines */
/* Not implemented yet (it seems to be slower) */ 
/* static const R_CallMethodDef Runuran_CallEntries[] = { */
/*     {"Runuran_init", (DL_FUNC) &Runuran_init, 2}, */
/*     {"Runuran_sample", (DL_FUNC) &Runuran_sample, 2}, */
/*     {NULL, NULL, 0} */
/* }; */

/*---------------------------------------------------------------------------*/

SEXP Runuran_discr_init (SEXP sexp_pv);
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for discrete distribution.           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
