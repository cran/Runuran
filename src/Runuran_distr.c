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

static void Runuran_distr_free(SEXP sexp_distr);
/*---------------------------------------------------------------------------*/
/* Free UNU.RAN distribution object.                                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* check pointer to generator object */
#define CHECK_PTR(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != Runuran_distr_tag) \
        error("[UNU.RAN - error] invalid UNU.RAN distribution object"); \
    } while (0)

/* Use an external reference to store the UNU.RAN generator objects */
static SEXP Runuran_distr_tag = NULL;

/*---------------------------------------------------------------------------*/

SEXP
Runuran_discr_init (SEXP sexp_pv)
     /*----------------------------------------------------------------------*/
     /* Create and initialize UNU.RAN object for discrete distribution.      */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_distr;
  struct unur_distr *distr;
  const double *pv;
  int n_pv;

  /* make tag for R object */
  if (!Runuran_distr_tag) Runuran_distr_tag = install("R_UNURAN_DISTR_TAG");

  /* check argument */
  if (!sexp_pv || TYPEOF(sexp_pv) != REALSXP)
    errorcall_return(R_NilValue,"[UNU.RAN - error] invalid argument 'pv'");

  /* get probability vector */
  pv = REAL(sexp_pv);
  n_pv = length(sexp_pv);

  /* create distribution object */
  distr = unur_distr_discr_new();
  if (unur_distr_discr_set_pv(distr,pv,n_pv) != UNUR_SUCCESS) {
    unur_distr_free(distr);
    errorcall_return(R_NilValue,"[UNU.RAN - error] cannot create UNU.RAN distribution object");
  }

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_distr = R_MakeExternalPtr(distr, Runuran_distr_tag, R_NilValue));
  UNPROTECT(1);
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_distr, Runuran_distr_free);

  /* return pointer to R */
  return (sexp_distr);

} /* end of Runuran_discr_init() */

/*---------------------------------------------------------------------------*/

void
Runuran_distr_free (SEXP sexp_distr)
     /*----------------------------------------------------------------------*/
     /* Free UNU.RAN distribution object.                                    */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr;

#ifdef DEBUG
  /* check pointer */
  CHECK_PTR(sexp_distr);
  printf("Runuran_distr_free called!\n");
#endif

  /* Extract pointer to distribution object */
  distr = R_ExternalPtrAddr(sexp_distr);

  /* free distribution object */
  unur_distr_free(distr);

  R_ClearExternalPtr(sexp_distr);

} /* end of Runuran_distr_free() */

/*---------------------------------------------------------------------------*/
