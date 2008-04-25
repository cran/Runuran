/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: Runuran_special.c                                                 *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         R interface for special UNU.RAN routines                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
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

/*---------------------------------------------------------------------------*/

SEXP
Runuran_qhinv (SEXP sexp_unur, SEXP sexp_U)
     /*----------------------------------------------------------------------*/
     /* Evaluate approximate quantile function when a UNU.RAN object of      */
     /* type HINV is given.                                                  */
     /*----------------------------------------------------------------------*/
{
  double *U;
  int n = 1;
  struct unur_gen *gen;
  int i;
  SEXP sexp_res = NULL;
  SEXP sexp_gen;
  SEXP sexp_slotunur;

#ifdef RUNURAN_DEBUG
  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] invalid UNU.RAN object");
#endif

  /* Extract and check U */
#ifdef RUNURAN_DEBUG
  /* check type of U */
  if (TYPEOF(sexp_U)!=REALSXP)
    error("[UNU.RAN - error] invalid argument 'U'");
#endif

  /* Extract U */
  U = NUMERIC_POINTER(sexp_U);
  n = length(sexp_U);

  /* Extract pointer to UNU.RAN generator */
  sexp_slotunur = Rf_install("unur");
  sexp_gen = GET_SLOT(sexp_unur, sexp_slotunur);
#ifdef RUNURAN_DEBUG
  /*   CHECK_PTR(sexp_gen); */
  /** CHECK_PTR is defined in Runuran.c **/
#endif
  gen = R_ExternalPtrAddr(sexp_gen);
#ifdef RUNURAN_DEBUG
  if (gen == NULL)
    error("[UNU.RAN - error] bad UNU.RAN object");
#endif
  if ( gen->method != UNUR_METH_HINV ) {
    error("[UNU.RAN - error] invalid UNU.RAN object: method 'HINV' required");
  }

  /* evaluate inverse CDF */
  PROTECT(sexp_res = NEW_NUMERIC(n));
  for (i=0; i<n; i++) {
    if (ISNAN(U[i]))
      /* if NA or NaN is given then we simply return the same value */
      NUMERIC_POINTER(sexp_res)[i] = U[i];
    else 
      NUMERIC_POINTER(sexp_res)[i] = unur_hinv_eval_approxinvcdf(gen,U[i]); 
  }
  UNPROTECT(1);

  /* return result to R */
  return sexp_res;
 
} /* end of Runuran_qhinv() */

/*---------------------------------------------------------------------------*/

