
/*---------------------------------------------------------------------------*/

#include "Runuran.h"
#include "Runuran_ext.h"

/*****************************************************************************/
/*                                                                           */
/*  Continuous Univariate Distributions (CONT)                               */
/*                                                                           */
/*****************************************************************************/

SEXP Runuran_ext_cont_init 
( SEXP sexp_obj, SEXP sexp_params, SEXP sexp_domain,
  UNUR_FUNCT_CONT *cdf, UNUR_FUNCT_CONT *pdf, UNUR_FUNCT_CONT *dpdf, int islog,
  double *mode, double *center, char *name )
/*---------------------------------------------------------------------------*/
/* Create and initialize UNU.RAN object for continuous distribution          */
/* defined in pure C code.                                                   */
/*                                                                           */
/* Parameters:                                                               */
/*   obj    ... S4 class that contains Runuran distribution object           */ 
/*   params ... vector of parameter values (R object)                        */
/*   domain ... domain of distribution (R object)                            */
/*   cdf    ... CDF of distribution                                          */
/*   pdf    ... PDF of distribution                                          */
/*   dpdf   ... derivative of PDF of distribution                            */
/*   islog  ... boolean: TRUE if logarithms of CDF|PDF|dPDF are given        */
/*   mode   ... mode of distribution                                         */
/*   center ... "center" (typical point) of distribution                     */
/*   name   ... name of distribution                                         */
/*---------------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct unur_distr *distr;
  const double *params;
  int n_params;
  const double *domain;
  unsigned int errcode = 0u;

  /* extract parameters */
  if (! (sexp_params && TYPEOF(sexp_params)==REALSXP) )
    Rf_error("[Runuran-Ext] invalid argument 'params'");
  params = REAL(sexp_params);
  n_params = Rf_length(sexp_params);

  /* extract domain of distribution */
  if (! (sexp_domain && TYPEOF(sexp_domain)==REALSXP && Rf_length(sexp_domain)==2) )
    Rf_error("[Runuran-Ext] invalid argument 'domain'");
  domain = REAL(sexp_domain);

  /* create continuous distribution object */
  distr = unur_distr_cont_new();
  if (distr == NULL)
    Rf_error("[Runuran-Ext] cannot create UNU.RAN object");

  /* set CDF, PDF and its derivative */
  if (islog) {
    if (cdf)  errcode |= unur_distr_cont_set_logcdf(distr, cdf);
    if (pdf)  errcode |= unur_distr_cont_set_logpdf(distr, pdf);
    if (dpdf) errcode |= unur_distr_cont_set_dlogpdf(distr, dpdf);
  }
  else {
    if (cdf)  errcode |= unur_distr_cont_set_cdf(distr, cdf);
    if (pdf)  errcode |= unur_distr_cont_set_pdf(distr, pdf);
    if (dpdf) errcode |= unur_distr_cont_set_dpdf(distr, dpdf);
  }

  /* set parameters */
  errcode |= unur_distr_cont_set_pdfparams(distr, params, n_params);

  /* set domain */
  errcode |= unur_distr_cont_set_domain( distr, domain[0], domain[1] );

  /* set mode, center and of distribution */
  if (mode)   errcode |= unur_distr_cont_set_mode(distr, *mode);
  if (center) errcode |= unur_distr_cont_set_center(distr, *center);
  if (name)   errcode |= unur_distr_set_name(distr, name);

  /* check return codes */
  if (errcode) {
    unur_distr_free (distr);
    Rf_error("[Runuran-Ext] cannot create UNU.RAN object");
  }

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, _Runuran_distr_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, _Runuran_distr_free);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_distr);

} /* end of Runuran_ext_cont_init() */

/*---------------------------------------------------------------------------*/
