/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: Runuran.c                                                         *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         R interface for UNU.RAN                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2007-2011 Wolfgang Hoermann and Josef Leydold             *
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

#include "Runuran.h"

/* internal header files for UNU.RAN */
#include <unur_source.h>

/*****************************************************************************/

SEXP
Runuran_init (SEXP sexp_obj, SEXP sexp_distr, SEXP sexp_method)
     /*----------------------------------------------------------------------*/
     /* Create and initialize UNU.RAN generator object.                      */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj    ... S4 class that contains 'Runuran' generator object       */ 
     /*   distr  ... distribution (string or S4 object)                      */
     /*   method ... method (string)                                         */
     /*                                                                      */
     /* Return:                                                              */
     /*   pointer to UNU.RAN generator object                                */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_gen;
  struct unur_gen *gen;
  const char *method;
  SEXP sexp_is_inversion;

  /* check argument */
  if (!sexp_method || TYPEOF(sexp_method) != STRSXP)
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'method'");
  if (!sexp_distr)
    /* this error should not happen. But we add it just in case. */
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid NULL pointer");

  /* get method string */
  method = CHAR(STRING_ELT(sexp_method,0));

  /* create generator object */
  switch (TYPEOF(sexp_distr)) {
  case STRSXP:
    {
      const char *distr = CHAR(STRING_ELT(sexp_distr,0));
      gen = unur_makegen_ssu( distr, method, NULL );
    }
    break;

  case EXTPTRSXP:
    {
      struct unur_distr *distr = R_ExternalPtrAddr(sexp_distr);
      gen = unur_makegen_dsu( distr, method, NULL );
    }
    break;

  default:
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'distribution'");
  }

  /* 'gen' must not be a NULL pointer */
  if (gen == NULL) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] cannot create UNU.RAN object");
  }

  /* set slot 'inversion' to true when 'gen' implements an inversion method. */
  PROTECT(sexp_is_inversion = NEW_LOGICAL(1));
  LOGICAL(sexp_is_inversion)[0] = unur_gen_is_inversion(gen);
  SET_SLOT(sexp_obj, install("inversion"), sexp_is_inversion);

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_gen = R_MakeExternalPtr(gen, _Runuran_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_gen, _Runuran_free);

  /* return pointer to R */
  UNPROTECT(2);
  return (sexp_gen);

} /* end of Runuran_init() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_sample (SEXP sexp_unur, SEXP sexp_n)
     /*----------------------------------------------------------------------*/
     /* Sample from UNU.RAN generator object.                                */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*   n    ... sample size (positive integer)                            */
     /*                                                                      */
     /* Return:                                                              */
     /*   random sample of size 'n'                                          */
     /*----------------------------------------------------------------------*/
{
  int n;
  SEXP sexp_gen;
  SEXP sexp_data;
  struct unur_gen *gen;

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* Extract and check sample size */
  n = *(INTEGER (AS_INTEGER (sexp_n)));
  if (n<=0) {
    error("sample size 'n' must be positive integer");
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  if (! isNull(sexp_gen)) {
    CHECK_UNUR_PTR(sexp_gen);
    gen = R_ExternalPtrAddr(sexp_gen);
    if (gen != NULL) {
      return _Runuran_sample_unur(gen,n);
    }
  }

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (! isNull(sexp_data)) {
    return _Runuran_sample_data(sexp_data,n);
  }

  /* Neither the UNU.RAN object nor the packed data list exists */
  errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
 
} /* end of Runuran_sample() */

/*---------------------------------------------------------------------------*/

SEXP
_Runuran_sample_unur (struct unur_gen *gen, int n)
     /*----------------------------------------------------------------------*/
     /* Sample from generator object: use UNU.RAN object                     */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   gen ... pointer to UNU.RAN generator object                        */
     /*   n   ... sample size (positive integer)                             */
     /*                                                                      */
     /* Return:                                                              */
     /*   random sample of size 'n'                                          */
     /*----------------------------------------------------------------------*/
{
/*   struct unur_gen *gen; */
  int i,k;
  SEXP sexp_res = R_NilValue;
  double *res;

  /* get state for the R built-in URNG */
  GetRNGstate();

  /* generate random vector of length n */
  switch (unur_distr_get_type(unur_get_distr(gen))) {

  case UNUR_DISTR_CONT:   /* univariate continuous distribution */
  case UNUR_DISTR_CEMP:   /* empirical continuous univariate distribution */
    PROTECT(sexp_res = NEW_NUMERIC(n));
    for (i=0; i<n; i++) {
      REAL(sexp_res)[i] = unur_sample_cont(gen); }
    break;

  case UNUR_DISTR_DISCR:  /* discrete univariate distribution */
    PROTECT(sexp_res = NEW_NUMERIC(n));
    for (i=0; i<n; i++) {
      REAL(sexp_res)[i] = (double) unur_sample_discr(gen); }
    break;

  case UNUR_DISTR_CVEC:   /* continuous mulitvariate distribution */
    {
      int dim = unur_get_dimension(gen);
      double *x = (double*) R_alloc(dim, sizeof(double) );
      PROTECT(sexp_res = allocMatrix(REALSXP, n, dim));
      res = REAL(sexp_res);
      for (i=0; i<n; i++) {
	if (unur_sample_vec(gen,x)!=UNUR_SUCCESS)
	  for (k=0; k<dim; k++) res[i + n*k] = NA_REAL;
	else
	  for (k=0; k<dim; k++) res[i + n*k] = x[k];
      }
    }
    break;

  case UNUR_DISTR_CVEMP:  /* empirical continuous multivariate distribution */
  case UNUR_DISTR_MATR:   /* matrix distribution */
  default:
    error("[UNU.RAN - error] '%s': Distribution type not support",
	  unur_distr_get_name(unur_get_distr(gen)) );
  }

  /* update state for the R built-in URNG */
  PutRNGstate();

  /* return result to R */
  UNPROTECT(1);
  return sexp_res;
 
} /* end of _Runuran_sample_unur() */

/*---------------------------------------------------------------------------*/

SEXP
_Runuran_sample_data (SEXP sexp_data, int n)
     /*----------------------------------------------------------------------*/
     /* Sample from generator object: use R data list (packed object)        */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   data ... data for generation method (R list)                       */
     /*   n    ... sample size (positive integer)                            */
     /*                                                                      */
     /* Return:                                                              */
     /*   random sample of size 'n'                                          */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_res = R_NilValue;
  int mid = INTEGER(VECTOR_ELT(sexp_data,0))[0];   /* method ID */

  /* get state for the R built-in URNG */
  GetRNGstate();

  switch (mid) {
  case UNUR_METH_PINV:
    PROTECT(sexp_res = _Runuran_sample_pinv(sexp_data,n));
    break;
  default:
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* update state for the R built-in URNG */
  PutRNGstate();

  /* return result to R */
  UNPROTECT(1);
  return sexp_res;
 
} /* end of _Runuran_sample_data() */

/*---------------------------------------------------------------------------*/

SEXP 
Runuran_quantile (SEXP sexp_unur, SEXP sexp_U)
     /*----------------------------------------------------------------------*/
     /* Evaluate approximate quantile function when a UNU.RAN object that    */
     /* implements an inversion method.                                      */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*   U    ... u-value (numeric array)                                   */
     /*                                                                      */
     /* Return:                                                              */
     /*   (approximate) quantiles for given 'U' values                       */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_is_inversion;
  SEXP sexp_gen;
  SEXP sexp_data;
  const char *class;               /* class name of 'unur' */ 
  struct unur_gen *gen;

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* check type of U */
  if (TYPEOF(sexp_U)!=REALSXP)
    error("[UNU.RAN - error] argument invalid: 'U' must be numeric");

  /* check class name */
  class = translateChar(STRING_ELT(GET_CLASS(sexp_unur), 0));
  if (strcmp(class,"unuran")) {
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");
  }

  /* check whether UNU.RAN object implements inversion method */
  sexp_is_inversion = GET_SLOT(sexp_unur, install("inversion"));
  if ( ! LOGICAL(sexp_is_inversion)[0]) {
    error("[UNU.RAN - error] invalid UNU.RAN object: inversion method required!\n\
\tUse methods 'HINV', 'NINV', 'PINV'; or 'DGT'");
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  if (! isNull(sexp_gen)) {
    CHECK_UNUR_PTR(sexp_gen);
    gen = R_ExternalPtrAddr(sexp_gen);
    if (gen != NULL) {
      return _Runuran_quantile_unur(gen,sexp_U);
    }
  }

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (! isNull(sexp_data)) {
    return _Runuran_quantile_data(sexp_data,sexp_U,sexp_unur);
  }

  /* Neither the UNU.RAN object nor the packed data list exists */
  errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");

} /* end of Runuran_quantile() */

/*---------------------------------------------------------------------------*/

SEXP
_Runuran_quantile_unur (struct unur_gen *gen, SEXP sexp_U)
     /*----------------------------------------------------------------------*/
     /* Evaluate approximate quantile function: use UNU.RAN object           */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   gen ... pointer to UNU.RAN generator object                        */
     /*   U    ... u-value (numeric array)                                   */
     /*                                                                      */
     /* Return:                                                              */
     /*   (approximate) quantiles for given 'U' values                       */
     /*----------------------------------------------------------------------*/
{
  double *U;
  int n = 1;
  int i;
  SEXP sexp_res = R_NilValue;

  /* Extract U */
  U = REAL(sexp_U);
  n = length(sexp_U);

  /* evaluate inverse CDF */
  PROTECT(sexp_res = NEW_NUMERIC(n));
  for (i=0; i<n; i++) {
    if (ISNAN(U[i]))
      /* if NA or NaN is given then we simply return the same value */
      REAL(sexp_res)[i] = U[i];
    else 
      REAL(sexp_res)[i] = unur_quantile(gen,U[i]); 
  }
  UNPROTECT(1);

  /* return result to R */
  return sexp_res;
 
} /* end of _Runuran_quantile_unur() */

/*---------------------------------------------------------------------------*/

SEXP
_Runuran_quantile_data (SEXP sexp_data, SEXP sexp_U, SEXP sexp_unur)
     /*----------------------------------------------------------------------*/
     /* Evaluate approximate quantile function:  use R data list (packed)    */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   data ... data for generation method (R list)                       */
     /*   U    ... u-value (numeric array)                                   */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*                                                                      */
     /* Return:                                                              */
     /*   (approximate) quantiles for given 'U' values                       */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_res = R_NilValue;
  int mid = INTEGER(VECTOR_ELT(sexp_data,0))[0];   /* method ID */

  switch (mid) {
  case UNUR_METH_PINV:
    return _Runuran_quantile_pinv(sexp_data,sexp_U,sexp_unur);
    break;
  default:
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* return result to R */
  return sexp_res;

} /* end of _Runuran_quantile_data() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_PDF (SEXP sexp_obj, SEXP sexp_x, SEXP sexp_islog)
     /*----------------------------------------------------------------------*/
     /* Evaluate PDF or PMF for UNU.RAN distribution or generator object.    */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj   ... 'Runuran' object (distribution or generator, S4 class)   */ 
     /*   x     ... x-value (numeric array)                                  */
     /*   islog ... boolean: if TRUE then the log-density is return          */
     /*                                                                      */
     /* Return:                                                              */
     /*   PDF or log-PDF in UNU.RAN object for given 'x' values              */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;                 /* S4 class containing distribution object */
  SEXP sexp_gen;                   /* S4 class containing generator object */
  SEXP sexp_res = R_NilValue;      /* R vector for storing result */
  struct unur_distr *distr = NULL; /* UNU.RAN distribution object */
  struct unur_gen *gen = NULL;     /* UNU.RAN generator object */
  const char *class;               /* class name of 'obj' */ 
  double *x;                       /* pointer to array arguments for PDF */
  int islog;                       /* whether we retur the log-density */
  int funct_missing = FALSE;
  int n = 1;
  int i;

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_obj))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* check type of x */
  if (TYPEOF(sexp_x)!=REALSXP && TYPEOF(sexp_x)!=INTSXP)
    error("[UNU.RAN - error] argument invalid: 'x' must be numeric");

  /* get class name */ 
  class = translateChar(STRING_ELT(GET_CLASS(sexp_obj), 0));

  /* which type of Runuran object */
  if (!strcmp(class,"unuran.cont") || !strcmp(class,"unuran.discr") ) {
    /* distribution object */
    sexp_distr = GET_SLOT(sexp_obj, install("distr"));
    CHECK_DISTR_PTR(sexp_distr);
    distr = R_ExternalPtrAddr(sexp_distr);
  }
  else if (!strcmp(class,"unuran")) {
    /* generator object */
    sexp_gen = GET_SLOT(sexp_obj, install("unur"));
    if (! isNull(sexp_gen)) {
      CHECK_UNUR_PTR(sexp_gen);
      gen = R_ExternalPtrAddr(sexp_gen);
      if (gen!=NULL) {
	distr = unur_get_distr(gen);
      }
    }
    if (distr==NULL) {
      SEXP sexp_data = GET_SLOT(sexp_obj, install("data"));
      if (! isNull(sexp_data)) {
	/* the generator object is packed */
	error("[UNU.RAN - error] cannot compute PDF for packed UNU.RAN object");
      }
      else {
	error("[UNU.RAN - error] broken UNU.RAN object");
      }
    }
  }
  else {
    error("[UNU.RAN - error] broken UNU.RAN object");
  }

  if ( ! ((distr->type == UNUR_DISTR_CONT) || (distr->type == UNUR_DISTR_DISCR)) ) {
    error("[UNU.RAN - error] invalid distribution type");
  }
  
  /* extract x */
  sexp_x = PROTECT(AS_NUMERIC(sexp_x));
  x = REAL(sexp_x);
  n = length(sexp_x);

  /* whether we have to return the log-density */
  islog = LOGICAL(sexp_islog)[0];

  /* check object for required functions */
  funct_missing = FALSE;
  if (distr->type == UNUR_DISTR_CONT) {
    if ( (islog  && distr->data.cont.logpdf == NULL) ||
	 (!islog && distr->data.cont.pdf == NULL) ) {
      funct_missing = TRUE;
      warning("[UNU.RAN - error] UNU.RAN object does not contain (log)PDF");
    }
  }
  if (distr->type == UNUR_DISTR_DISCR) {
    if ( (islog)   /* not implemented yet */
	 || distr->data.discr.pmf == NULL) {
      funct_missing = TRUE;
      warning("[UNU.RAN - error] UNU.RAN object does not contain (log)PMF");
    }
  }

  /* allocate memory for result */
  PROTECT(sexp_res = NEW_NUMERIC(n));

  /* evaluate CDF */
  for (i=0; i<n; i++) {

    if (funct_missing) {
      /* function not implemented */
      REAL(sexp_res)[i] = NA_REAL;
      continue;
    }

    if (ISNAN(x[i])) {
      /* if NA or NaN is given then we simply return the same value */
      REAL(sexp_res)[i] = x[i];
      continue;
    }

    /* switch (unur_distr_get_type(distr)) { */
    switch (distr->type) {
    case UNUR_DISTR_CONT:
      /* univariate continuous distribution  --> evaluate PDF */
      REAL(sexp_res)[i] = (islog)
	? unur_distr_cont_eval_logpdf(x[i], distr)
	: unur_distr_cont_eval_pdf(x[i], distr);
      break;

    case UNUR_DISTR_DISCR:
      /* discrete univariate distribution  --> evaluate PMF */
      if (x[i] < INT_MIN || x[i] > INT_MAX) 
	REAL(sexp_res)[i] = 0.;
      else
	REAL(sexp_res)[i] = unur_distr_discr_eval_pmf ((int) x[i], distr);
      /* remark: logPMF yet not implemented */
      break;

    default:
      /* this code should not be reachable */
      error("[UNU.RAN - error] internal error");
    }
  }
  
  /* return result to R */
  UNPROTECT(2);
  return sexp_res;

} /* end of Runuran_PDF() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_CDF (SEXP sexp_obj, SEXP sexp_x)
     /*----------------------------------------------------------------------*/
     /* Evaluate CDF for UNU.RAN distribution or generator object.           */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   obj ... 'Runuran' object (distribution or generator, S4 class)     */ 
     /*   x   ... x-value (numeric array)                                    */
     /*                                                                      */
     /* Return:                                                              */
     /*   CDF in UNU.RAN object for given 'x' values                         */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;                 /* S4 class containing distribution object */
  SEXP sexp_gen;                   /* S4 class containing generator object */
  SEXP sexp_res = R_NilValue;      /* R vector for storing result */
  struct unur_distr *distr = NULL; /* UNU.RAN distribution object */
  struct unur_gen *gen = NULL;     /* UNU.RAN generator object */
  const char *class;               /* class name of 'obj' */ 
  double *x;                       /* pointer to array arguments for PDF */
  int n = 1;
  int i;

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_obj))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* check type of x */
  if (TYPEOF(sexp_x)!=REALSXP && TYPEOF(sexp_x)!=INTSXP)
    error("[UNU.RAN - error] argument invalid: 'x' must be numeric");

  /* get class name */
  class = translateChar(STRING_ELT(GET_CLASS(sexp_obj), 0));

  /* which type of Runuran object */
  if (!strcmp(class,"unuran.cont") || !strcmp(class,"unuran.discr") ) {
    /* distribution object */
    sexp_distr = GET_SLOT(sexp_obj, install("distr"));
    CHECK_DISTR_PTR(sexp_distr);
    distr = R_ExternalPtrAddr(sexp_distr);
  }
  else if (!strcmp(class,"unuran")) {
    /* generator object */
    sexp_gen = GET_SLOT(sexp_obj, install("unur"));
    if (! isNull(sexp_gen)) {
      CHECK_UNUR_PTR(sexp_gen);
      gen = R_ExternalPtrAddr(sexp_gen);
      if (gen!=NULL) {
	distr = unur_get_distr(gen);
      }
    }
    if (distr==NULL) {
      SEXP sexp_data = GET_SLOT(sexp_obj, install("data"));
      if (! isNull(sexp_data)) {
	/* the generator object is packed */
	error("[UNU.RAN - error] cannot compute CDF for packed UNU.RAN object");
      }
      else {
	error("[UNU.RAN - error] broken UNU.RAN object");
      }
    }
  }
  else {
    error("[UNU.RAN - error] broken UNU.RAN object");
  }
  
  if ( ! ((distr->type == UNUR_DISTR_CONT) || (distr->type == UNUR_DISTR_DISCR)) ) {
    error("[UNU.RAN - error] invalid distribution type");
  }

  /* check objects */
  if (distr->type == UNUR_DISTR_DISCR && distr->data.discr.cdf == NULL)
      error("[UNU.RAN - error] UNU.RAN object does not contain CDF");

  if (distr->type == UNUR_DISTR_CONT && distr->data.cont.cdf == NULL) {
    if (gen==NULL)
      error("[UNU.RAN - error] UNU.RAN object does not contain CDF");
    else if (gen->method != UNUR_METH_PINV)
      error("[UNU.RAN - error] function requires method PINV");
  }

  /* extract x */
  sexp_x = PROTECT(AS_NUMERIC(sexp_x));
  x = REAL(sexp_x);
  n = length(sexp_x);

  /* allocate memory for result */
  PROTECT(sexp_res = NEW_NUMERIC(n));

  /* evaluate CDF */
  for (i=0; i<n; i++) {
    if (ISNAN(x[i])) {
      /* if NA or NaN is given then we simply return the same value */
      REAL(sexp_res)[i] = x[i];
      continue;
    }

    switch (distr->type) {
    case UNUR_DISTR_CONT:
      /* univariate continuous distribution */
      if (distr->data.cont.cdf != NULL)
	REAL(sexp_res)[i] = unur_distr_cont_eval_cdf(x[i], distr);
      else
	REAL(sexp_res)[i] = unur_pinv_eval_approxcdf(gen, x[i]);
      break;

    case UNUR_DISTR_DISCR:
      /* discrete univariate distribution */
      if (x[i] < INT_MIN) 
	REAL(sexp_res)[i] = 0.;
      else if (x[i] > INT_MAX) 
	REAL(sexp_res)[i] = 1.;
      else
	REAL(sexp_res)[i] = unur_distr_discr_eval_cdf ((int) x[i], distr);
      break;

    default:
      /* this code should not be reachable */
      error("[UNU.RAN - error] internal error");
    }
  }

  /* return result to R */
  UNPROTECT(2);
  return sexp_res;

} /* end of Runuran_CDF() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_print (SEXP sexp_unur, SEXP sexp_help)
     /*----------------------------------------------------------------------*/
     /* Print information about UNU.RAN generator object.                    */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*   help ... whether to print detailed information (integer / boolean) */
     /*                                                                      */
     /* Return:                                                              */
     /*   info string                                                        */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_gen;
  SEXP sexp_data;
  struct unur_gen *gen = NULL;
  int help;
  const char *info;
  SEXP sexp_info;

  /* These two variables are used for replacing URNG in generator object */
  /* while executing unur_gen_info(). However, this might change the     */
  /* URNG permanently when the useR hits CRTL-C during execution!        */
  /* So it is not used, yet!                                             */
  /* UNUR_URNG *urng_tmp, *urng_aux_tmp; */

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (! isNull(sexp_data)) {
    Rprintf("Object is PACKED !\n\n");
    return R_NilValue;
  }

  /* Extract slot for UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  if (isNull(sexp_gen)) {
    warningcall_immediate(R_NilValue,"[UNU.RAN - warning] empty UNU.RAN object");
    return R_NilValue;
  }
  ALLWAYS_CHECK_UNUR_PTR(sexp_gen);

  /* Extract pointer to UNU.RAN generator */
  gen = R_ExternalPtrAddr(sexp_gen);
  if (gen == NULL) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* Extract help flag */
  help = *(INTEGER (AS_INTEGER (sexp_help)));

  /* Replace URNG by UNU.RAN built-in URNG */
  /* urng_aux_tmp = urng_tmp = unur_urng_builtin(); */
  /* urng_aux_tmp = unur_chg_urng(gen, urng_aux_tmp ); */
  /* urng_tmp = unur_chg_urng(gen, urng_tmp ); */

  /* get state for the R built-in URNG */
  GetRNGstate();

  /* get info string */
  info = unur_gen_info(gen,help);
  if (info==NULL) { 
    /* this should not happen. but we want to protect against NULL. */
    info = ""; 
  }

  /* update state for the R built-in URNG */
  PutRNGstate();

  /* restore URNG in generator object */
  /* urng_aux_tmp = unur_chg_urng(gen, urng_aux_tmp ); */
  /* urng_tmp = unur_chg_urng(gen, urng_tmp ); */
  /* unur_urng_free(urng_tmp); */

  /* create R string */
  PROTECT(sexp_info = mkString(info));
  UNPROTECT(1);
  return sexp_info;

} /* end of Runuran_print() */

/*---------------------------------------------------------------------------*/

void
_Runuran_free (SEXP sexp_gen)
     /*----------------------------------------------------------------------*/
     /* Free UNU.RAN generator object.                                       */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   gen ... pointer to UNU.RAN generator object                        */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

#ifdef RUNURAN_DEBUG
  /* check pointer */
  CHECK_UNUR_PTR(sexp_gen);
  Rprintf("Runuran_free called!\n");
#endif

  /* Extract pointer to generator */
  gen = R_ExternalPtrAddr(sexp_gen);
  if (gen==NULL) {
    /* UNU.RAN object already destroyed: nothing to do */
    return;
  }

  /* free generator object */
  unur_free(gen);
  R_ClearExternalPtr(sexp_gen);

} /* end of _Runuran_free() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_pack (SEXP sexp_unur)
     /*----------------------------------------------------------------------*/
     /* Pack Runuran generator object into R list                            */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*                                                                      */
     /* Return:                                                              */
     /*   data for generation method (R list)                                */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  SEXP sexp_gen;
  SEXP sexp_data;

  /* argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (! isNull(sexp_data)) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] object already packed");
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  CHECK_UNUR_PTR(sexp_gen);
  if (isNull(sexp_gen) || 
      ((gen=R_ExternalPtrAddr(sexp_gen)) == NULL) ) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* call packing subroutine for given generator object */
  switch (unur_get_method(gen)) {
  case UNUR_METH_PINV:
    _Runuran_pack_pinv(gen, sexp_unur);
    break;

  default:
    errorcall_return(R_NilValue,"[UNU.RAN - error] cannot pack UNU.RAN object");
  }

  /* remove UNU.RAN object */
  unur_free(gen);
  R_ClearExternalPtr(sexp_gen);

  /* o.k. */
  return R_NilValue;

} /* end of Runuran_pack() */

/*---------------------------------------------------------------------------*/

SEXP Runuran_set_error_level (SEXP sexp_level)
     /*----------------------------------------------------------------------*/
     /* Set verbosity level of UNU.RAN error handler.                        */
     /* Pack Runuran generator object into R list                            */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   level ... verbosity level (integer)                                */ 
     /*                                                                      */
     /* Return:                                                              */
     /*   old level                                                          */
     /*----------------------------------------------------------------------*/
{
  int level;
  int old_level;
  SEXP sexp_old_level;

  /* Extract and check verbosity level */
  level = *(INTEGER (AS_INTEGER (sexp_level)));
  if (level < 0L || level > 3L) {
    error("verbosity 'level' of UNU.RAN error handler must be 0, 1, 2, or 3");
  }

  /* set new verbosity level */
  old_level = _Runuran_set_error_handler(level);

  /* return old one */
  PROTECT(sexp_old_level = NEW_INTEGER(1));
  INTEGER(sexp_old_level)[0] = old_level;
  UNPROTECT(1);

  return sexp_old_level;

} /* end of Runuran_set_error_level() */

/*---------------------------------------------------------------------------*/

int
_Runuran_set_error_handler(int level)
     /*----------------------------------------------------------------------*/
     /* set verbosity level of UNU.RAN error handler.                        */
     /*                                                                      */
     /* Parameter:                                                           */
     /*   level ... 0 = suppress all warnings / error messages               */
     /*             1 = print errors only                                    */
     /*             2 = print (almost all) warnings and error messages       */
     /*             3 = print all warnings and error messages                */
     /*                                                                      */
     /* Return:                                                              */
     /*   old verbosity level (integer)                                      */
     /*----------------------------------------------------------------------*/
{
  int old_level;
  UNUR_ERROR_HANDLER *old_handler;

  /* Switch UNU.RAN error handler */
  switch(level) {
  case 0:
    old_handler = unur_set_error_handler( _Runuran_error_handler_suppress );
    break;
  case 1:
    old_handler = unur_set_error_handler( _Runuran_error_handler_error );
    break;
  case 2:
  default:
    old_handler = unur_set_error_handler( _Runuran_error_handler_warning );
    break;
  case 3:
    old_handler = unur_set_error_handler( _Runuran_error_handler_print );
    break;
  }

  /* Return level of old error handler */
  if (old_handler == _Runuran_error_handler_suppress) {
    old_level = 0L;
  }
  else if (old_handler == _Runuran_error_handler_error) {
    old_level = 1L;
  }
  else if (old_handler == _Runuran_error_handler_warning) {
    old_level = 2L;
  }
  else if (old_handler == _Runuran_error_handler_print) {
    old_level = 3L;
  }
  else {
    old_level = RUNURAN_DEFAULT_ERROR_HANDLER_LEVEL;
  }

  return old_level;
}

/*---------------------------------------------------------------------------*/

void
_Runuran_error_handler_print( const char *objid     ATTRIBUTE__UNUSED,
			      const char *file      ATTRIBUTE__UNUSED,
			      int line              ATTRIBUTE__UNUSED,
			      const char *errortype,
			      int errorcode,
			      const char *reason )   
     /*----------------------------------------------------------------------*/
     /* Error handler for UNU.RAN routines:                                  */
     /* Show all errors and warnings                                         */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... source file name (__FILE__)                          */
     /*   line      ... source line number (__LINE__)                        */
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNU.RAN error code                                   */
     /*   reason    ... short description of reason                          */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* print warning or error message */
  Rprintf("[UNU.RAN - %s] %s",errortype,unur_get_strerror(errorcode));
  if (reason && strlen(reason))
    Rprintf(": %s\n", reason);
  else
    Rprintf("\n");
  
#ifdef RUNURAN_DEBUG
  /* print file and line number */
  Rprintf("\tfile: %s, line: %d\n",file,line);
#endif

#ifdef UNUR_ENABLE_LOGGING
  /* print error message into log file (use UNU.RAN default error handler) */
  _unur_error_handler_default( objid, file, line, errortype, errorcode, reason );
#endif

} /* end of _Runuran_error_handler_print() */

/*---------------------------------------------------------------------------*/

void
_Runuran_error_handler_warning( const char *objid     ATTRIBUTE__UNUSED,
				const char *file      ATTRIBUTE__UNUSED,
				int line              ATTRIBUTE__UNUSED,
				const char *errortype,
				int errorcode,
				const char *reason )   
     /*----------------------------------------------------------------------*/
     /* Error handler for UNU.RAN routines:                                  */
     /* Show errors and most of the warnings                                 */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... source file name (__FILE__)                          */
     /*   line      ... source line number (__LINE__)                        */
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNU.RAN error code                                   */
     /*   reason    ... short description of reason                          */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* we suppress some warnings */
  if (errortype[0] == 'w') {
    switch (errorcode) {
      /* we do not print warnings for the following codes: */
    case UNUR_ERR_DISTR_REQUIRED:
      return;

    default:
      break;
    }
  }

  /* print error message or warning */ 
  _Runuran_error_handler_print(objid, file, line, errortype, errorcode, reason);

} /* end of _Runuran_error_handler_warning() */

/*---------------------------------------------------------------------------*/

void
_Runuran_error_handler_error( const char *objid     ATTRIBUTE__UNUSED,
			      const char *file      ATTRIBUTE__UNUSED,
			      int line              ATTRIBUTE__UNUSED,
			      const char *errortype,
			      int errorcode,
			      const char *reason )   
     /*----------------------------------------------------------------------*/
     /* Error handler for UNU.RAN routines:                                  */
     /* Show only errors                                                     */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... source file name (__FILE__)                          */
     /*   line      ... source line number (__LINE__)                        */
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNU.RAN error code                                   */
     /*   reason    ... short description of reason                          */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* we suppress some warnings */
  if (errortype[0] == 'w') {
      return;
  }

  /* print error message or warning */ 
  _Runuran_error_handler_print(objid, file, line, errortype, errorcode, reason);

} /* end of _Runuran_error_handler_error() */

/*---------------------------------------------------------------------------*/

void
_Runuran_error_handler_suppress( const char *objid     ATTRIBUTE__UNUSED,
				 const char *file      ATTRIBUTE__UNUSED,
				 int line              ATTRIBUTE__UNUSED,
				 const char *errortype ATTRIBUTE__UNUSED,
				 int errorcode         ATTRIBUTE__UNUSED,
				 const char *reason    ATTRIBUTE__UNUSED )
     /*----------------------------------------------------------------------*/
     /* Error handler that suppresses all warnings/errors.                   */
     /* Error handler for UNU.RAN routines                                   */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... source file name (__FILE__)                          */
     /*   line      ... source line number (__LINE__)                        */
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNU.RAN error code                                   */
     /*   reason    ... short description of reason                          */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  ; /* nothing to do */
} /* end of _Runuran_error_handler_suppress() */

/*---------------------------------------------------------------------------*/

double
_Runuran_R_unif_rand (void *unused  ATTRIBUTE__UNUSED)
     /*----------------------------------------------------------------------*/
     /* Wrapper for R built-in uniform random number generator               */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unused ... argument required for UNU.RAN API                       */
     /*                                                                      */
     /* Return:                                                              */
     /*   uniform random number                                              */ 
     /*----------------------------------------------------------------------*/
{
  /* TODO (?): Calling GetRNGstate() and PutRNGstate() is very slow! */
  /* So we have called in Runuran_sample(). (Good?) */
  /*   double x; */
  /*   GetRNGstate(); */
  /*   x = unif_rand(); */
  /*   PutRNGstate(); */
  /*   return x; */

  return unif_rand();
} /* end _Runuran_R_unif_rand() */

/*---------------------------------------------------------------------------*/

SEXP _Runuran_tag(void) 
     /*----------------------------------------------------------------------*/
     /* Make tag for R UNU.RAN generator object                              */
     /*                                                                      */
     /* Parameters: none                                                     */
     /*                                                                      */
     /* Return:                                                              */
     /*   tag (R object)                                                     */ 
     /*----------------------------------------------------------------------*/
{
  static SEXP tag = NULL;

  /* make tag for R object */
  if (!tag) tag = install("R_UNURAN_TAG");

  return tag;
} /* end _Runuran_tag() */

/*---------------------------------------------------------------------------*/

SEXP Runuran_use_aux_urng (SEXP sexp_unur, SEXP sexp_set)
     /*----------------------------------------------------------------------*/
     /* check, set or unset auxiliary URNG for given generator object        */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur ... 'Runuran' generator object                                */ 
     /*   set  ... whether we use an auxiliary URNG                          */
     /*                                                                      */
     /* Return:                                                              */
     /*   old value of 'set'                                                 */ 
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_gen;                   /* S4 class containing generator object */
  struct unur_gen *gen = NULL;     /* UNU.RAN generator object */
  int set;                         /* value which we have to set */
  SEXP sexp_old = R_NilValue;      /* old value of set */
  const char *class;               /* class name of 'unr' */ 

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN generator object");

  /* we need a generator object */ 
  class = translateChar(STRING_ELT(GET_CLASS(sexp_unur), 0));
  if (strcmp(class,"unuran")) 
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN generator object");

  /* extract pointer to UNU.RAN generator object */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  if (! isNull(sexp_gen)) {
    CHECK_UNUR_PTR(sexp_gen);
    gen = R_ExternalPtrAddr(sexp_gen);
    if (gen == NULL) {
      error("[UNU.RAN - error] broken UNU.RAN object");
    }
  }

  /* read old value of 'set' */
  PROTECT(sexp_old = NEW_LOGICAL(1));
  if (unur_get_urng_aux(gen) == NULL) {
    LOGICAL(sexp_old)[0] = NA_LOGICAL;
  }
  else {
    LOGICAL(sexp_old)[0] = (unur_get_urng(gen) == unur_get_urng_aux(gen)) ? FALSE : TRUE;
  }
  UNPROTECT(1);

  /* set new value */
  if (! isNull(sexp_set)) {
    set = LOGICAL(sexp_set)[0];
    if (unur_get_urng_aux(gen) == NULL) {
      error("[UNU.RAN - error] generator object does not support auxiliary URNG");
    }
    else {
      if (set)
	unur_chgto_urng_aux_default(gen);
      else
	unur_chg_urng_aux(gen,unur_get_urng(gen));
    }
  }

  /* return old value of 'set' */
  return (sexp_old);

} /* end of Runuran_use_aux_urng() */

/*---------------------------------------------------------------------------*/

SEXP Runuran_set_aux_seed (SEXP sexp_seed)
     /*----------------------------------------------------------------------*/
     /* set seed for auxiliary URNG                                          */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   seed  ... seed for auxiliary URNG                                  */ 
     /*                                                                      */
     /* Return:                                                              */
     /*   NULL                                                               */ 
     /*----------------------------------------------------------------------*/
{
  unsigned long seed = INTEGER(sexp_seed)[0];
  if (seed <= 0) 
    error("[UNU.RAN - error] seed is non-positive");
  unur_urng_seed (unur_get_default_urng_aux(), seed);
  return R_NilValue;
} /* end of Runuran_set_aux_seed() */

/*---------------------------------------------------------------------------*/
