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
 *   Copyright (c) 2007 Wolfgang Hoermann and Josef Leydold                  *
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

/*---------------------------------------------------------------------------*/

/* #define DEBUG 1 */

/*---------------------------------------------------------------------------*/

static void Runuran_free(SEXP sexp_gen);
/*---------------------------------------------------------------------------*/
/* Free UNU.RAN generator object.                                            */
/*---------------------------------------------------------------------------*/

static void Runuran_error_handler( 
	const char *objid, const char *file, int line,
        const char *errortype, int errorcode, const char *reason );
/*---------------------------------------------------------------------------*/
/* Error handler for UNU.RAN routines.                                       */
/*---------------------------------------------------------------------------*/

static double Runuran_R_unif_rand (void *unused);
/*---------------------------------------------------------------------------*/
/* Wrapper for R built-in uniform random number generator.                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* check pointer to generator object */
#define CHECK_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != Runuran_tag) \
        error("[UNU.RAN - error] invalid UNU.RAN object"); \
    } while (0)

/* Use an external reference to store the UNU.RAN generator objects */
static SEXP Runuran_tag = NULL;

/*---------------------------------------------------------------------------*/

SEXP
Runuran_init (SEXP sexp_distr, SEXP sexp_method)
     /*----------------------------------------------------------------------*/
     /* Create and initialize UNU.RAN generator object.                      */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_gen;
  struct unur_gen *gen;
  const char *method;

  /* make tag for R object */
  if (!Runuran_tag) Runuran_tag = install("R_UNURAN_TAG");

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
  PROTECT(sexp_gen = R_MakeExternalPtr(gen, Runuran_tag, R_NilValue));
  UNPROTECT(1);
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_gen, Runuran_free);

  /* return pointer to R */
  return (sexp_gen);

} /* end of Runuran_init() */

/*---------------------------------------------------------------------------*/

SEXP
Runuran_sample (SEXP sexp_gen, SEXP sexp_n)
     /*----------------------------------------------------------------------*/
     /* Sample from UNU.RAN generator object.                                */
     /*----------------------------------------------------------------------*/
{
  int n;
  struct unur_gen *gen;
  int i;
  SEXP res = NULL;

#ifdef DEBUG
  /* check pointer */
  CHECK_PTR(sexp_gen);
#endif

  /* Extract and check sample size */
  n = *(INTEGER (AS_INTEGER (sexp_n)));
  if (n<=0) {
    error("sample size 'n' must be positive integer");
  }

  /* Extract pointer to generator */
  gen = R_ExternalPtrAddr(sexp_gen);
  
#ifdef DEBUG
  /* this must not be a NULL pointer */
  if (gen == NULL)
    error("[UNU.RAN - error] bad UNU.RAN object");
#endif

  /* TODO: this must not be called when we do not use the R built-in URNG */
  GetRNGstate();

  /* generate random vector of length n */
  switch (unur_distr_get_type(unur_get_distr(gen))) {

  case UNUR_DISTR_CONT:   /* univariate continuous distribution */
  case UNUR_DISTR_CEMP:   /* empirical continuous univariate distribution */
    PROTECT(res = NEW_NUMERIC(n));
    for (i=0; i<n; i++) {
      NUMERIC_POINTER(res)[i] = unur_sample_cont(gen); }
    UNPROTECT(1);
    break;

  case UNUR_DISTR_DISCR:  /* discrete univariate distribution */
    PROTECT(res = NEW_NUMERIC(n));
    for (i=0; i<n; i++) {
      NUMERIC_POINTER(res)[i] = (double) unur_sample_discr(gen); }
    UNPROTECT(1);
    break;

  case UNUR_DISTR_CVEC:   /* continuous mulitvariate distribution */
  case UNUR_DISTR_CVEMP:  /* empirical continuous multivariate distribution */
  case UNUR_DISTR_MATR:   /* matrix distribution */
  default:
    error("[UNU.RAN - error] '%s': Distribution type yet not support",
	  unur_distr_get_name(unur_get_distr(gen)) );
  }

  /* TODO: this must not be called when we do not use the R built-in URNG */
  PutRNGstate();

  /* return result to R */
  return res;
 
} /* end of Runuran_sample() */

/*---------------------------------------------------------------------------*/

void
Runuran_free (SEXP sexp_gen)
     /*----------------------------------------------------------------------*/
     /* Free UNU.RAN generator object.                                       */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

#ifdef DEBUG
  /* check pointer */
  CHECK_PTR(sexp_gen);
  printf("Runuran_free called!\n");
#endif

  /* Extract pointer to generator */
  gen = R_ExternalPtrAddr(sexp_gen);

  /* free generator object */
  unur_free(gen);

  R_ClearExternalPtr(sexp_gen);

} /* end of Runuran_free() */

/*---------------------------------------------------------------------------*/

void 
R_init_Runuran (DllInfo *info) 
     /*----------------------------------------------------------------------*/
     /* Initialization routine when loading the DLL                          */
     /*----------------------------------------------------------------------*/
{
  /* Set new UNU.RAN error handler */
  unur_set_error_handler( Runuran_error_handler );

  /* Set R built-in generator as default URNG */
  unur_set_default_urng( unur_urng_new( Runuran_R_unif_rand, NULL) );
  unur_set_default_urng_aux( unur_urng_new( Runuran_R_unif_rand, NULL) );

  /* Register native routines */ 
  /* Not implemented yet */ 
  /*   R_registerRoutines(info, NULL, Runuran_CallEntries, NULL, NULL); */
  /*   R_useDynamicSymbols(info, FALSE); */

} /* end of R_init_Runuran() */

/*---------------------------------------------------------------------------*/

void
Runuran_error_handler( 
	const char *objid,     /* id/type of object              */
        const char *file,      /* source file name (__FILE__)    */
        int line,              /* source line number (__LINE__)  */ 
        const char *errortype, /* "warning" or "error"           */
        int errorcode,         /* UNU.RAN error code             */
        const char *reason     /* short description of reason    */
     )   
     /*----------------------------------------------------------------------*/
     /* Error handler for UNU.RAN routines                                   */
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

  /* print warning or error message */
  Rprintf("[UNU.RAN - %s] %s",errortype,unur_get_strerror(errorcode));
  if (reason && strlen(reason))
    Rprintf(": %s\n", reason);
  else
    Rprintf("\n");
  
#ifdef DEBUG
  /* print file and line number */
  Rprintf("\tfile: %s, line: %d\n",file,line);
#endif

} /* end of Runuran_error_handler() */

/*---------------------------------------------------------------------------*/

double
Runuran_R_unif_rand (void *unused)
     /*----------------------------------------------------------------------*/
     /* Wrapper for R built-in uniform random number generator               */
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
} /* end Runuran_R_unif_rand() */

/*---------------------------------------------------------------------------*/
