/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: verify.c                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         Check UNU.RAN generator object.                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2010 Wolfgang Hoermann and Josef Leydold                  *
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

#include "Runuran.h"

/*---------------------------------------------------------------------------*/
/* run tests                                                                 */

static int run_verify_hat(struct unur_gen *gen, int n);


/*****************************************************************************/

SEXP
Runuran_verify_hat (SEXP sexp_unur, SEXP sexp_n)
/*---------------------------------------------------------------------------*/
/* Verify hat of UNU.RAN generator object that implements rejection method.  */
/*                                                                           */
/* Parameters:                                                               */
/*   unur ... 'Runuran' object (S4 class)                                    */ 
/*   n    ... sample size (positive integer)                                 */
/*                                                                           */
/* Return:                                                                   */
/*   ratio failed / samplesize   (numeric)                                   */
/*---------------------------------------------------------------------------*/
{
  int n;                     /* sample size */
  SEXP sexp_gen;             /* R object containing UNU.RAN generator object */
  struct unur_gen *gen;      /* UNU.RAN generator object */    
  SEXP sexp_failed = R_NilValue;  /* result (ratio) */

  /* function for toggling verify mode on/off */
  int (*chg_verify)(UNUR_GEN *generator,int verify);

  /* first argument must be S4 class */
  if (!IS_S4_OBJECT(sexp_unur))
    error("[UNU.RAN - error] argument invalid: 'unr' must be UNU.RAN object");

  /* UNU.RAN must not be packed */
  if (! isNull( GET_SLOT(sexp_unur, install("data")) ) ) {
    /* the generator object is packed */
    error("[UNU.RAN - error] cannot run this function on packed UNU.RAN objects");
  }

  /* Extract and check sample size */
  n = *(INTEGER (AS_INTEGER (sexp_n)));
  if (n<=0) {
    error("sample size 'n' must be positive integer");
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  CHECK_UNUR_PTR(sexp_gen);
  if (isNull(sexp_gen) || 
      ((gen=R_ExternalPtrAddr(sexp_gen)) == NULL) ) {
    errorcall_return(R_NilValue,"[UNU.RAN - error] broken UNU.RAN object");
  }

  /* get function for toggling verify mode */

#define METHOD(meth)  chg_verify = unur_##meth##_chg_verify; break;

  switch (unur_get_method(gen)) {

    /* continuous univariate distributions */
  case UNUR_METH_AROU:   METHOD(arou);
  case UNUR_METH_ARS:    METHOD(ars);
  case UNUR_METH_HRB:    METHOD(hrb);
  case UNUR_METH_HRD:    METHOD(hrd);
  case UNUR_METH_HRI:    METHOD(hri);
  case UNUR_METH_ITDR:   METHOD(itdr);
  case UNUR_METH_NROU:   METHOD(nrou);
  case UNUR_METH_SROU:   METHOD(srou);
  case UNUR_METH_SSR:    METHOD(ssr);
  case UNUR_METH_TABL:   METHOD(tabl);
  case UNUR_METH_TDR:    METHOD(tdr);
  case UNUR_METH_UTDR:   METHOD(utdr);

    /* discrete univariate distributions */
  case UNUR_METH_DARI:   METHOD(dari);
  case UNUR_METH_DSROU:  METHOD(dsrou);

    /* multivariate distributions */
  case UNUR_METH_MVTDR:  METHOD(mvtdr);
  case UNUR_METH_VNROU:  METHOD(vnrou);

  default:
    errorcall_return(R_NilValue,"[UNU.RAN - error] Method not supported");
  }

#undef METHOD

  /* create R object for storing result */
  PROTECT(sexp_failed = NEW_INTEGER(1));

  /* run generator in verify mode and store result */
  chg_verify(gen,TRUE);
  INTEGER_POINTER(sexp_failed)[0] = run_verify_hat(gen,n);
  chg_verify(gen,FALSE);

  /* return ratio 'failed' / 'sample size' to R */
  UNPROTECT(1);
  return sexp_failed;

} /* end of Runuran_verify_hat() */

/*---------------------------------------------------------------------------*/

int
run_verify_hat(struct unur_gen *gen, int n)
{
  int i;
  int dim;                   /* dimension of distribution object */
  double *x = NULL;
  int failed = 0;

  /* get state for the R built-in URNG */
  GetRNGstate();

  /* get dimension of distribution */
  dim = unur_get_dimension (gen);
  
  /* allocate working space */
  if (dim > 0) {
      x = (double*) R_alloc(dim, sizeof(double) );
  }

  /* switch off error messages */
  _Runuran_set_error_handler(0);

  /* run generator */
  for (i=0; i<n; i++) {

    /* reset errno */
    unur_reset_errno();
    switch(unur_distr_get_type(unur_get_distr(gen))) {
    case UNUR_DISTR_CONT:   /* univariate continuous distribution */
      unur_sample_cont(gen);
      break;
    
    case UNUR_DISTR_DISCR:  /* discrete univariate distribution */
      unur_sample_discr(gen);
      break;

    case UNUR_DISTR_CVEC:   /* continuous mulitvariate distribution */
      unur_sample_vec(gen,x);
      break;

    case UNUR_DISTR_CEMP:   /* empirical continuous univariate distribution */
    case UNUR_DISTR_CVEMP:  /* empirical continuous multivariate distribution */
    case UNUR_DISTR_MATR:   /* matrix distribution */
    default:
      _Runuran_set_error_handler(1);
      error("[UNU.RAN - error] '%s': Distribution type not support",
	    unur_distr_get_name(unur_get_distr(gen)) );
    }

    /* check for sampling error */
    if (unur_get_errno()) {
      /* == UNUR_ERR_GEN_CONDITION */
      failed++;
    }
  }

  /* switch on error messages */
  _Runuran_set_error_handler(1);
  
  /* update state for the R built-in URNG */
  PutRNGstate();

  /* portion of failed samples */
  return failed;

} /* end of run_verify_hat() */

/*---------------------------------------------------------------------------*/

