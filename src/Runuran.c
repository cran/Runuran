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

#include "Runuran.h"
#include <unuran.h>

#include <unur_source.h>
#include <methods/unur_methods_source.h>

#include "Runuran_source.h"

/*---------------------------------------------------------------------------*/

static SEXP _Runuran_sample_unur (struct unur_gen *gen, int n);
/*---------------------------------------------------------------------------*/
/* Sample from generator object: use UNU.RAN object                          */
/*---------------------------------------------------------------------------*/

static SEXP _Runuran_sample_data (SEXP sexp_data, int n);
/*---------------------------------------------------------------------------*/
/* Sample from generator object: use R data list (packed object)             */
/*---------------------------------------------------------------------------*/

static SEXP _Runuran_quantile_unur (struct unur_gen *gen, SEXP sexp_U);
/*---------------------------------------------------------------------------*/
/* Evaluate approximate quantile function: use UNU.RAN object                */
/*---------------------------------------------------------------------------*/

static SEXP _Runuran_quantile_data (SEXP sexp_data, SEXP sexp_U, SEXP sexp_unur);
/*---------------------------------------------------------------------------*/
/* Evaluate approximate quantile function:  use R data list (packed object)  */
/*---------------------------------------------------------------------------*/

static void _Runuran_free(SEXP sexp_gen);
/*---------------------------------------------------------------------------*/
/* Free UNU.RAN generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _Runuran_error_handler( 
	const char *objid, const char *file, int line,
        const char *errortype, int errorcode, const char *reason );
/*---------------------------------------------------------------------------*/
/* Error handler for UNU.RAN routines.                                       */
/*---------------------------------------------------------------------------*/

static double _Runuran_R_unif_rand (void *unused);
/*---------------------------------------------------------------------------*/
/* Wrapper for R built-in uniform random number generator.                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* check pointer to generator object */
#define CHECK_UNUR_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != _Runuran_tag) \
        error("[UNU.RAN - error] invalid UNU.RAN object"); \
    } while (0)

/* Use an external reference to store the UNU.RAN generator objects */
static SEXP _Runuran_tag = NULL;

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

  /* make tag for R object */
  if (!_Runuran_tag) _Runuran_tag = install("R_UNURAN_TAG");

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

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_gen = R_MakeExternalPtr(gen, _Runuran_tag, sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_gen, _Runuran_free);

  /* return pointer to R */
  UNPROTECT(1);
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
#ifdef RUNURAN_DEBUG
  CHECK_UNUR_PTR(sexp_gen);
#endif
  gen = R_ExternalPtrAddr(sexp_gen);
  if (gen != NULL) {
    return _Runuran_sample_unur(gen,n);
  }

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (TYPEOF(sexp_data)!=NILSXP) {
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
      NUMERIC_POINTER(sexp_res)[i] = unur_sample_cont(gen); }
    break;

  case UNUR_DISTR_DISCR:  /* discrete univariate distribution */
    PROTECT(sexp_res = NEW_NUMERIC(n));
    for (i=0; i<n; i++) {
      NUMERIC_POINTER(sexp_res)[i] = (double) unur_sample_discr(gen); }
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
    sexp_res = _Runuran_sample_pinv(sexp_data,n);
    break;
  default:
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* update state for the R built-in URNG */
  PutRNGstate();

  /* return result to R */
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
  SEXP sexp_gen;
  SEXP sexp_data;
  struct unur_gen *gen;

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* check type of U */
  if (TYPEOF(sexp_U)!=REALSXP)
    error("[UNU.RAN - error] argument invalid: 'U' must be number or vector");

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
#ifdef RUNURAN_DEBUG
  CHECK_UNUR_PTR(sexp_gen);
#endif
  gen = R_ExternalPtrAddr(sexp_gen);
  if (gen != NULL) {
    return _Runuran_quantile_unur(gen,sexp_U);
  }

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (TYPEOF(sexp_data)!=NILSXP) {
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
  U = NUMERIC_POINTER(sexp_U);
  n = length(sexp_U);

  /* check whether UNU.RAN object implements inversion method */
  switch (gen->method) {
  case UNUR_METH_DGT:
  case UNUR_METH_HINV:
  case UNUR_METH_NINV:
  case UNUR_METH_PINV:
    break;
  default:
    error("[UNU.RAN - error] invalid UNU.RAN object: inversion method required!\n\tUse methods 'HINV', 'NINV', 'PINV'; or 'DGT'");
  }

  /* evaluate inverse CDF */
  PROTECT(sexp_res = NEW_NUMERIC(n));
  for (i=0; i<n; i++) {
    if (ISNAN(U[i]))
      /* if NA or NaN is given then we simply return the same value */
      NUMERIC_POINTER(sexp_res)[i] = U[i];
    else 
      NUMERIC_POINTER(sexp_res)[i] = unur_quantile(gen,U[i]); 
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
Runuran_print (SEXP sexp_unur, SEXP sexp_help)
     /*----------------------------------------------------------------------*/
     /* Print information about UNU.RAN generator object.                    */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*   help ... whether to print detailed information (integer / boolean) */
     /*                                                                      */
     /* Return:                                                              */
     /*   NilValue  [since (void) is not possible]                           */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_gen;
  SEXP sexp_data;
  struct unur_gen *gen = NULL;
  int help;
  const char *info;

  /* Extract data list */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (TYPEOF(sexp_data)!=NILSXP) {
    Rprintf("Object is PACKED !\n\n");
    return R_NilValue;
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
#ifdef RUNURAN_DEBUG
  CHECK_UNUR_PTR(sexp_gen);
#endif
  gen = R_ExternalPtrAddr(sexp_gen);
  if (gen == NULL) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* Extract help flag */
  help = *(INTEGER (AS_INTEGER (sexp_help)));

  /* get info string */
  info = unur_gen_info(gen,help);

  /* print info string */
  if (info) {
    Rprintf("%s",info);
  }
  /*   else { */
  /*     /\* no info string available *\/ */
  /*     Rprintf("nix\n");  */
  /*   } */
  
  /* nothing to return */
  return R_NilValue;
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
  if (TYPEOF(sexp_data)!=NILSXP) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] object already packed");
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
#ifdef RUNURAN_DEBUG
  CHECK_UNUR_PTR(sexp_gen);
#endif
  gen = R_ExternalPtrAddr(sexp_gen);
  if (gen == NULL) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* call packing subroutine for given generator object */
  switch (gen->method) {
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

void 
R_init_Runuran (DllInfo *info  ATTRIBUTE__UNUSED) 
     /*----------------------------------------------------------------------*/
     /* Initialization routine when loading the DLL                          */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   info ... passed by R and describes the DLL                         */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* Set new UNU.RAN error handler */
  unur_set_error_handler( _Runuran_error_handler );

  /* Set R built-in generator as default URNG */
  unur_set_default_urng( unur_urng_new( _Runuran_R_unif_rand, NULL) );
  unur_set_default_urng_aux( unur_urng_new( _Runuran_R_unif_rand, NULL) );

  /* Register native routines */ 
  /* Not implemented yet */ 
  /*   R_registerRoutines(info, NULL, Runuran_CallEntries, NULL, NULL); */
  /*   R_useDynamicSymbols(info, FALSE); */

} /* end of R_init_Runuran() */

/*---------------------------------------------------------------------------*/

void
R_unload_Runuran (DllInfo *info  ATTRIBUTE__UNUSED)
     /*----------------------------------------------------------------------*/
     /* Clear memory before unloading the DLL.                               */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   info ... passed by R and describes the DLL                         */
     /*                                                                      */
     /* Return:                                                              */
     /*   (void)                                                             */
     /*----------------------------------------------------------------------*/
{
  unur_urng_free(unur_get_default_urng());
  unur_urng_free(unur_get_default_urng_aux());
} /* end of R_unload_Runuran() */

/*---------------------------------------------------------------------------*/

void
_Runuran_error_handler( const char *objid     ATTRIBUTE__UNUSED,
			const char *file      ATTRIBUTE__UNUSED,
			int line              ATTRIBUTE__UNUSED,
			const char *errortype,
			int errorcode,
			const char *reason )   
     /*----------------------------------------------------------------------*/
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
#ifdef RUNURAN_DEBUG
#else
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
#endif

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

} /* end of _Runuran_error_handler() */

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
