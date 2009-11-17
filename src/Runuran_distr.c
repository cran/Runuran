/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: Runuran_distr.c                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         R interface to UNU.RAN distribution objects                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2007-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <unuran.h>
#include "Runuran.h"

/*---------------------------------------------------------------------------*/

/* structure for storing pointers to R objects                               */
struct Runuran_distr_discr {
  SEXP env;                 /* R environment                                 */
  SEXP pmf;                 /* PMF of distribution                           */
};

struct Runuran_distr_cont {
  SEXP env;                 /* R environment                                 */
  SEXP cdf;                 /* CDF of distribution                           */
  SEXP pdf;                 /* PDF of distribution                           */
  SEXP dpdf;                /* derivative of PDF of distribution             */
};

struct Runuran_distr_cmv {
  SEXP env;                 /* R environment                                 */
  SEXP pdf;                 /* PDF of distribution                           */
};

/*---------------------------------------------------------------------------*/
/*  Discrete Distributions (DISCR)                                           */

static double _Runuran_discr_eval_pmf( int k, const struct unur_distr *distr );
/* Evaluate PMF function.                                                    */

/*---------------------------------------------------------------------------*/
/*  Continuous Univariate Distributions (CONT)                               */

static double _Runuran_cont_eval_cdf( double x, const struct unur_distr *distr );
/* Evaluate CDF function.                                                    */

static double _Runuran_cont_eval_pdf( double x, const struct unur_distr *distr );
/* Evaluate PDF function.                                                    */

static double _Runuran_cont_eval_dpdf( double x, const struct unur_distr *distr );
/* Evaluate derivative of PDF function.                                      */

/*---------------------------------------------------------------------------*/
/*  Continuous Multivariate Distributions (CMV)                              */

static double _Runuran_cmv_eval_pdf( const double *x, struct unur_distr *distr );
/* Evaluate PDF function.                                                    */

/*---------------------------------------------------------------------------*/

static SEXP _Runuran_distr_tag(void); 
/*---------------------------------------------------------------------------*/
/* Make tag for R object [Contains static variable!]                         */
/*---------------------------------------------------------------------------*/

#define CHECK_UNUR_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != _Runuran_distr_tag()) \
      error("[UNU.RAN - error] invalid UNU.RAN distribution object");	\
  } while (0)
/*---------------------------------------------------------------------------*/
/* Check pointer to R UNU.RAN distribution object.                           */
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/*                                                                           */
/*  Discrete Distributions (DISCR)                                           */
/*                                                                           */
/*****************************************************************************/

SEXP
Runuran_discr_init (SEXP sexp_obj, SEXP sexp_env,
		    SEXP sexp_pv, SEXP sexp_pmf,
		    SEXP sexp_mode, SEXP sexp_domain,
		    SEXP sexp_sum, SEXP sexp_name)
     /*----------------------------------------------------------------------*/
     /* Create and initialize UNU.RAN object for discrete distribution.      */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj    ... S4 class that contains unuran distribution object       */ 
     /*   env    ... R environment                                           */
     /*   pv     ... PV of distribution                                      */
     /*   pmf    ... PMF of distribution                                     */
     /*   mode   ... mode of distribution                                    */
     /*   domain ... domain of distribution                                  */
     /*   sum    ... sum over PV / PMF                                       */
     /*   name   ... name of distribution                                    */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct Runuran_distr_discr *Rdistr;
  struct unur_distr *distr;
  const double *pv;
  int n_pv;
  const double *domain;
  double mode, sum;
  int lb, ub;
  const char *name;
  unsigned int error = 0u;

  /* domain of distribution */
  if (! (sexp_domain && TYPEOF(sexp_domain)==REALSXP && length(sexp_domain)==2) )
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'domain'");
  domain = REAL(sexp_domain);
  lb = (domain[0] < (double) INT_MIN) ? INT_MIN : (int) domain[0]; 
  ub = (domain[1] > (double) INT_MAX) ? INT_MAX : (int) domain[1];
  if (lb >= ub)
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid domain: lb >= ub");

  /* create distribution object */
  distr = unur_distr_discr_new();
  if (distr == NULL) _Runuran_fatal();

  /* set domain */
  error |= unur_distr_discr_set_domain( distr, lb, ub );

  /* set probability vector */
  if (!isNull(sexp_pv)) {
    pv = REAL(AS_NUMERIC(sexp_pv));
    if (ISNAN(pv[0]))
      errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'pv'");
    n_pv = length(sexp_pv);
    error |= unur_distr_discr_set_pv(distr,pv,n_pv);
  }

  /* store pointers to R objects */
  Rdistr = Calloc(1,struct Runuran_distr_discr);
  Rdistr->env = sexp_env;
  Rdistr->pmf = sexp_pmf;

  /* set function pointers */
  error |= unur_distr_set_extobj(distr, Rdistr);
  if (!isNull(sexp_pmf)) {
    error |= unur_distr_discr_set_pmf(distr, _Runuran_discr_eval_pmf);
  }

  /* set mode, center and PDFarea of distribution */
  mode = REAL(AS_NUMERIC(sexp_mode))[0];
  sum = REAL(AS_NUMERIC(sexp_sum))[0];

  if (!ISNAN(mode))
    error |= unur_distr_discr_set_mode(distr, (int) mode);
  if (!ISNAN(sum))
    error |= unur_distr_discr_set_pmfsum(distr, sum);

  /* set name of distribution */
  if (sexp_name && TYPEOF(sexp_name) == STRSXP) {
    name = CHAR(STRING_ELT(sexp_name,0));
    unur_distr_set_name(distr,name);
  }
  /* else we simply ignore the 'name' argument */

  /* check return codes */
  if (error) {
    Free(Rdistr);
    unur_distr_free (distr); 
    _Runuran_fatal();
  } 

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, _Runuran_distr_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, _Runuran_distr_free);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_distr);

} /* end of Runuran_discr_init() */

/*---------------------------------------------------------------------------*/

double
_Runuran_discr_eval_pmf( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* Evaluate PMF function.                                               */
     /*----------------------------------------------------------------------*/
{
  const struct Runuran_distr_discr *Rdistr;
  SEXP R_fcall, arg;
  double y;

  Rdistr = unur_distr_get_extobj(distr);
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = (double)k;
  PROTECT(R_fcall = lang2(Rdistr->pmf, R_NilValue));
  SETCADR(R_fcall, arg);
  y = REAL(eval(R_fcall, Rdistr->env))[0];
  UNPROTECT(2);
  return y;
} /* end of _Runuran_discr_eval_pmf() */


/*****************************************************************************/
/*                                                                           */
/*  Continuous Univariate Distributions (CONT)                               */
/*                                                                           */
/*****************************************************************************/

SEXP
Runuran_cont_init (SEXP sexp_obj, SEXP sexp_env, 
		   SEXP sexp_cdf, SEXP sexp_pdf, SEXP sexp_dpdf, SEXP sexp_islog,
		   SEXP sexp_mode, SEXP sexp_center, SEXP sexp_domain,
		   SEXP sexp_area, SEXP sexp_name)
     /*----------------------------------------------------------------------*/
     /* Create and initialize UNU.RAN object for continuous distribution.    */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj    ... S4 class that contains unuran distribution object       */ 
     /*   env    ... R environment                                           */
     /*   cdf    ... CDF of distribution                                     */
     /*   pdf    ... PDF of distribution                                     */
     /*   dpdf   ... derivative of PDF of distribution                       */
     /*   islog  ... boolean: TRUE if logarithms of CDF|PDF|dPDF are given   */
     /*   mode   ... mode of distribution                                    */
     /*   center ... "center" (typical point) of distribution                */
     /*   domain ... domain of distribution                                  */
     /*   area   ... area below PDF                                          */
     /*   name   ... name of distribution                                    */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct Runuran_distr_cont *Rdistr;
  struct unur_distr *distr;
  const double *domain;
  double mode, center, area;
  const char *name;
  int islog;
  unsigned int error = 0u;

#ifdef RUNURAN_DEBUG
  /*   /\* 'this' must be an S4 class *\/ */
  /*   if (!IS_S4_OBJECT(sexp_this)) */
  /*     errorcall_return(R_NilValue,"[UNU.RAN - error] invalid object"); */

  /* all other variables are tested in the R routine */
  /* TODO: add checks in DEBUGging mode */
#endif

  /* domain of distribution */
  if (! (sexp_domain && TYPEOF(sexp_domain)==REALSXP && length(sexp_domain)==2) )
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'domain'");
  domain = REAL(sexp_domain);

  /* whether we are given logarithm of CDF|PDF|dPDF or not */
  islog = LOGICAL(sexp_islog)[0];

  /* store pointers to R objects */
  Rdistr = Calloc(1,struct Runuran_distr_cont);
  Rdistr->env = sexp_env;
  Rdistr->cdf = sexp_cdf;
  Rdistr->pdf = sexp_pdf;
  Rdistr->dpdf = sexp_dpdf;

  /* create distribution object */
  distr = unur_distr_cont_new();
  if (distr == NULL) _Runuran_fatal();

  /* set domain */
  error |= unur_distr_cont_set_domain( distr, domain[0], domain[1] );

  /* set function pointers */
  error |= unur_distr_set_extobj(distr, Rdistr);
  if (islog) {
    if (!isNull(sexp_cdf))
      error |= unur_distr_cont_set_logcdf(distr, _Runuran_cont_eval_cdf);
    if (!isNull(sexp_pdf))
      error |= unur_distr_cont_set_logpdf(distr, _Runuran_cont_eval_pdf);
    if (!isNull(sexp_dpdf))
      error |= unur_distr_cont_set_dlogpdf(distr, _Runuran_cont_eval_dpdf);
  }
  else {
    if (!isNull(sexp_cdf))
      error |= unur_distr_cont_set_cdf(distr, _Runuran_cont_eval_cdf);
    if (!isNull(sexp_pdf))
      error |= unur_distr_cont_set_pdf(distr, _Runuran_cont_eval_pdf);
    if (!isNull(sexp_dpdf))
      error |= unur_distr_cont_set_dpdf(distr, _Runuran_cont_eval_dpdf);
  }

  /* set mode, center and PDFarea of distribution */
  mode = REAL(AS_NUMERIC(sexp_mode))[0];
  center = REAL(AS_NUMERIC(sexp_center))[0];
  area = REAL(AS_NUMERIC(sexp_area))[0];

  if (!ISNAN(mode))
    error |= unur_distr_cont_set_mode(distr, mode);
  if (!ISNAN(center))
    error |= unur_distr_cont_set_center(distr, center);
  if (!ISNAN(area))
    error |= unur_distr_cont_set_pdfarea(distr, area);

  /* set name of distribution */
  if (sexp_name && TYPEOF(sexp_name) == STRSXP) {
    name = CHAR(STRING_ELT(sexp_name,0));
    unur_distr_set_name(distr,name);
  }
  /* else we simply ignore the 'name' argument */

  /* check return codes */
  if (error) {
    Free(Rdistr);
    unur_distr_free (distr); 
    _Runuran_fatal();
  } 

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, _Runuran_distr_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, _Runuran_distr_free);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_distr);

} /* end of Runuran_cont_init() */

/*---------------------------------------------------------------------------*/

double
_Runuran_cont_eval_cdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* Evaluate CDF function.                                               */
     /*----------------------------------------------------------------------*/
{
  const struct Runuran_distr_cont *Rdistr;
  SEXP R_fcall, arg;
  double y;

  Rdistr = unur_distr_get_extobj(distr);
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = x;
  PROTECT(R_fcall = lang2(Rdistr->cdf, R_NilValue));
  SETCADR(R_fcall, arg);
  y = REAL(eval(R_fcall, Rdistr->env))[0];
  UNPROTECT(2);
  return y;
} /* end of _Runuran_cont_eval_cdf() */

/*---------------------------------------------------------------------------*/

double
_Runuran_cont_eval_pdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* Evaluate PDF function.                                               */
     /*----------------------------------------------------------------------*/
{
  const struct Runuran_distr_cont *Rdistr;
  SEXP R_fcall, arg;
  double y;

  Rdistr = unur_distr_get_extobj(distr);
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = x;
  PROTECT(R_fcall = lang2(Rdistr->pdf, R_NilValue));
  SETCADR(R_fcall, arg);
  y = REAL(eval(R_fcall, Rdistr->env))[0];
  UNPROTECT(2);
  return y;
} /* end of _Runuran_cont_eval_pdf() */

/*---------------------------------------------------------------------------*/

double
_Runuran_cont_eval_dpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* Evaluate derivative of PDF function.                                 */
     /*----------------------------------------------------------------------*/
{
  const struct Runuran_distr_cont *Rdistr;
  SEXP R_fcall, arg;
  double y;

  Rdistr = unur_distr_get_extobj(distr);
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = x;
  PROTECT(R_fcall = lang2(Rdistr->dpdf, R_NilValue));
  SETCADR(R_fcall, arg);
  y = REAL(eval(R_fcall, Rdistr->env))[0];
  UNPROTECT(2);
  return y;
} /* end of _Runuran_cont_eval_dpdf() */


/*****************************************************************************/
/*                                                                           */
/*  Continuous Multivariate Distributions (CMV)                              */
/*                                                                           */
/*****************************************************************************/

SEXP
Runuran_cmv_init (SEXP sexp_obj, SEXP sexp_env, 
		  SEXP sexp_dim, SEXP sexp_pdf,
		  SEXP sexp_mode, SEXP sexp_center,
		  SEXP sexp_ll, SEXP sexp_ur, SEXP sexp_name)
     /*----------------------------------------------------------------------*/
     /* Create and initialize UNU.RAN object for continuous multivariate     */
     /* distribution.                                                        */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj    ... S4 class that contains unuran distribution object       */ 
     /*   env    ... R environment                                           */
     /*   dim    ... dimensions of distribution                              */
     /*   pdf    ... PDF of distribution                                     */
     /*   mode   ... mode of distribution                                    */
     /*   center ... center of distribution                                  */
     /*   ll, ur ... lower left and upper right vertex of rectangular domain */
     /*   name   ... name of distribution                                    */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct Runuran_distr_cmv *Rdistr;
  struct unur_distr *distr;
  const int *dim;
  const double *mode, *center;
  const double *ll, *ur;
  const char *name;
  unsigned int error = 0u;

#ifdef RUNURAN_DEBUG
  /*   /\* 'this' must be an S4 class *\/ */
  /*   if (!IS_S4_OBJECT(sexp_this)) */
  /*     errorcall_return(R_NilValue,"[UNU.RAN - error] invalid object"); */

  /* all other variables are tested in the R routine */
  /* TODO: add checks in DEBUGging mode */
#endif

  /* number of dimensions */
  dim = INTEGER(sexp_dim);

  /* store pointers to R objects */
  Rdistr = Calloc(1,struct Runuran_distr_cmv);
  Rdistr->env = sexp_env;
  Rdistr->pdf = sexp_pdf;

  /* create distribution object */
  distr = unur_distr_cvec_new(dim[0]);
  if (distr == NULL) _Runuran_fatal();

  /* set function pointers */
  error |= unur_distr_set_extobj(distr, Rdistr);
  if (!isNull(sexp_pdf))
    error |= unur_distr_cvec_set_pdf(distr, _Runuran_cmv_eval_pdf);

  /* set domain */
  if (!isNull(sexp_ll) && !isNull(sexp_ur)) {
    ll = REAL(sexp_ll);
    ur = REAL(sexp_ur);
    error |= unur_distr_cvec_set_domain_rect(distr, ll, ur);
  }  
  
  /* set mode */
  if (!isNull(sexp_mode)) {
    mode = REAL(sexp_mode);
    error |= unur_distr_cvec_set_mode(distr, mode);
  }

  /* set center */
  if (!isNull(sexp_center)) {
    center = REAL(sexp_center);
    error |= unur_distr_cvec_set_center(distr, center);
  }

  /* set name of distribution */
  if (sexp_name && TYPEOF(sexp_name) == STRSXP) {
    name = CHAR(STRING_ELT(sexp_name,0));
    unur_distr_set_name(distr,name);
  }
  /* else we simply ignore the 'name' argument */

  /* check return codes */
  if (error) {
    Free(Rdistr);
    unur_distr_free (distr); 
    _Runuran_fatal();
  } 

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, _Runuran_distr_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, _Runuran_distr_free);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_distr);

} /* end of Runuran_cmv_init() */

/*---------------------------------------------------------------------------*/

double
_Runuran_cmv_eval_pdf( const double *x, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* Evaluate PDF function.                                               */
     /*----------------------------------------------------------------------*/
{
  const struct Runuran_distr_cmv *Rdistr;
  SEXP R_fcall, arg;
  double *rarg;
  double y;
  int i, dim;
  
  /* get dimension of distribution */
  dim = unur_distr_get_dim(distr);

  /* pointer to R object */
  Rdistr = unur_distr_get_extobj(distr);

  /* copy x into R object of type "numeric" */
  PROTECT(arg = NEW_NUMERIC(dim));
  rarg = REAL(arg);
  for (i=0; i<dim; i++)
    rarg[i] = x[i];

  /* evaluate PDF */
  PROTECT(R_fcall = lang2(Rdistr->pdf, R_NilValue));
  SETCADR(R_fcall, arg);
  y = REAL(eval(R_fcall, Rdistr->env))[0];
  UNPROTECT(2);

  return y;
} /* end of _Runuran_cmv_eval_pdf() */


/*****************************************************************************/
/*                                                                           */
/* Special distributions                                                     */
/*                                                                           */
/*****************************************************************************/

SEXP
Runuran_std_cont (SEXP sexp_obj, SEXP sexp_name, SEXP sexp_params, SEXP sexp_domain)
     /*----------------------------------------------------------------------*/
     /* Create UNU.RAN object for special continuous distribution.           */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj    ... S4 class that contains unuran distribution object       */ 
     /*   name   ... name of special distribution                            */
     /*   params ... vector of parameter values                              */
     /*   domain ... domain of distribution                                  */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct unur_distr *distr;
  const char *name;
  const double *params;
  int n_params;
  const double *domain;
  unsigned int error = 0u;

  /* name of distribution */
  if (! (sexp_name && TYPEOF(sexp_name) == STRSXP) && length(sexp_name)>2)
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'name'");
  name = CHAR(STRING_ELT(sexp_name,0));

  /* parameters */
  if (! (sexp_params && TYPEOF(sexp_params)==REALSXP) )
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'params'");
  params = REAL(sexp_params);
  n_params = length(sexp_params);

  /* domain of distribution */
  if (! (sexp_domain && TYPEOF(sexp_domain)==REALSXP && length(sexp_domain)==2) )
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'domain'");
  domain = REAL(sexp_domain);

  /* create distribution object */
  distr = _Runuran_get_std_cont( name, params, n_params );
  if (distr == NULL) { _Runuran_fatal(); }
  
  /* set domain */
  error |= unur_distr_cont_set_domain( distr, domain[0], domain[1] );

  /* check return codes */
  if (error) {
    unur_distr_free (distr);
    _Runuran_fatal();
  }

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, _Runuran_distr_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, _Runuran_distr_free);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_distr);
} /* end of Runuran_std_cont() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_std_discr (SEXP sexp_obj, SEXP sexp_name, SEXP sexp_params, SEXP sexp_domain)
     /*----------------------------------------------------------------------*/
     /* Create UNU.RAN object for special discrete distribution.             */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj    ... S4 class that contains unuran distribution object       */ 
     /*   name   ... name of special distribution                            */
     /*   params ... vector of parameter values                              */
     /*   domain ... domain of distribution                                  */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct unur_distr *distr;
  const char *name;
  const double *params;
  int n_params;
  const double *domain;
  int lb, ub;
  unsigned int error = 0u;

  /* name of distribution */
  if (! (sexp_name && TYPEOF(sexp_name) == STRSXP) && length(sexp_name)>2)
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'name'");
  name = CHAR(STRING_ELT(sexp_name,0));

  /* parameters */
  if (! (sexp_params && TYPEOF(sexp_params)==REALSXP && length(sexp_params)>0) )
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'params'");
  params = REAL(sexp_params);
  n_params = length(sexp_params);

  /* domain of distribution */
  if (! (sexp_domain && TYPEOF(sexp_domain)==REALSXP && length(sexp_domain)==2) )
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'domain'");
  domain = REAL(sexp_domain);
  lb = (domain[0] < (double) INT_MIN) ? INT_MIN : (int) domain[0];
  ub = (domain[1] > (double) INT_MAX) ? INT_MAX : (int) domain[1];

  /* create distribution object */
  distr = _Runuran_get_std_discr( name, params, n_params );
  if (distr == NULL) { _Runuran_fatal(); }
  
  /* set domain */
  error |= unur_distr_discr_set_domain( distr, lb, ub );

  /* check return codes */
  if (error) {
    unur_distr_free (distr);
    _Runuran_fatal();
  }

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, _Runuran_distr_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, _Runuran_distr_free);

  /* return pointer to R */
  UNPROTECT(1);
  return (sexp_distr);
} /* end of Runuran_std_discr() */


/*****************************************************************************/
/*                                                                           */
/*  Common Routines                                                          */
/*                                                                           */
/*****************************************************************************/

void
_Runuran_distr_free (SEXP sexp_distr)
     /*----------------------------------------------------------------------*/
     /* Free UNU.RAN distribution object.                                    */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr;
  const void *Rdistr;

#ifdef DEBUG
  /* check pointer */
  CHECK_PTR(sexp_distr);
  Rprintf("Runuran_distr_free called!\n");
#endif

  /* Extract pointer to distribution object */
  distr = R_ExternalPtrAddr(sexp_distr);

  /* free structure that stores R object */
  Rdistr = unur_distr_get_extobj(distr);
  Free(Rdistr);

  /* free distribution object */
  unur_distr_free(distr);

  R_ClearExternalPtr(sexp_distr);

} /* end of _Runuran_distr_free() */

/*---------------------------------------------------------------------------*/

SEXP _Runuran_distr_tag(void) 
     /*----------------------------------------------------------------------*/
     /* Make tag for R UNU.RAN distribution object                           */
     /*                                                                      */
     /* Parameters: none                                                     */
     /*                                                                      */
     /* Return:                                                              */
     /*   tag (R object)                                                     */ 
     /*----------------------------------------------------------------------*/
{
  static SEXP tag = NULL;

  /* make tag for R object */
  if (!tag) tag = install("R_UNURAN_DISTR_TAG");

  return tag;
} /* end _Runuran_distr_tag() */

/*---------------------------------------------------------------------------*/
