/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: Runuran_pinv.c                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         R interface for UNU.RAN -- PINV                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2009 Wolfgang Hoermann and Josef Leydold                  *
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

/* internal header files for UNU.RAN */
#include <unur_source.h>
#include <methods/pinv_struct.h>

/*---------------------------------------------------------------------------*/

/* number of entries (slots) in data list */
#define n_slots (5)

/* names of slots */
static const char *slot_name[n_slots] = {"mid","order","Umax","guide","iv"};

/* positions in data list */
enum {
  pmid = 0,      /* method ID [ This MUST be 0 ! ] */
  porder = 1,    /* order of polynomial */
  pUmax = 2,     /* Umax */
  pguide = 3,    /* guide table */
  piv = 4        /* coefficents of polynomials */
};

/*---------------------------------------------------------------------------*/

static double _pinv_eval (double U, double Umax, int order,
			  int guide_size, int *guide, double *iv);
/*---------------------------------------------------------------------------*/
/* Evaluate approximating polynomial.                                        */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

#define GEN    ((struct unur_pinv_gen*)gen->datap)
/* data for generator object */

#define DISTR  (gen->distr->data.cont)
/* data for distribution in generator object */

/*****************************************************************************/

void
_Runuran_pack_pinv (struct unur_gen *gen, SEXP sexp_unur)
     /*----------------------------------------------------------------------*/
     /* Pack Runuran generator object for method PINV into R list            */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   gen  ... pointer to UNU.RAN generator object                       */
     /*   unur ... 'Runuran' object (S4 class)                               */ 
     /*                                                                      */
     /* Return:                                                              */
     /*   data for generation method (R list)                                */
     /*----------------------------------------------------------------------*/
{
  int i,n,k;
  int iv_size, n_coeff;
  double *iv;

  /* names of list entries */
  SEXP sexp_data_names;

  /* data list and its entries */
  SEXP sexp_data, sexp_dom;
  SEXP sexp_mid, sexp_order, sexp_Umax, sexp_guide, sexp_iv;

  /* number of doubles stored for one interval */
  iv_size = 1 + 2*GEN->order;

  /* create entries for data list */

  /* method ID (int) */
  PROTECT(sexp_mid = Rf_allocVector(INTSXP, 1));
  INTEGER(sexp_mid)[0] = UNUR_METH_PINV;

  /* order (int) */
  PROTECT(sexp_order = Rf_allocVector(INTSXP, 1));
  INTEGER(sexp_order)[0] = GEN->order;

  /* Umax (double) */
  PROTECT(sexp_Umax = Rf_allocVector(REALSXP, 1));
  REAL(sexp_Umax)[0] = GEN->Umax;

  /* guide table (int[]) */
  PROTECT(sexp_guide = Rf_allocVector(INTSXP, GEN->guide_size));
  for (i=0; i<GEN->guide_size; i++) {
    INTEGER(sexp_guide)[i] = iv_size*GEN->guide[i];
  }

  /* total number of coefficients for polynomials */
  n_coeff = (GEN->n_ivs+1) * iv_size;

  /* table of coefficients for approximating polynomial */
  /* sequence for each interval: 
   *   cdfi, z[order-1], u[order-2], z[order-2], ..., u[0], z[0], xi  
   */
  PROTECT(sexp_iv = Rf_allocVector(REALSXP, n_coeff));
  iv = REAL(sexp_iv);
  for (i=0,n=-1; i<=GEN->n_ivs; i++) {
    iv[++n] = GEN->iv[i].cdfi;
    k = GEN->order - 1;
    iv[++n] = GEN->iv[i].zi[k];
    for (k--; k>=0; k--) {
      iv[++n] = GEN->iv[i].ui[k];
      iv[++n] = GEN->iv[i].zi[k];
    }
    iv[++n] = GEN->iv[i].xi;
  }

  /* list of "names" attribute of the objects in our list */
  PROTECT(sexp_data_names = Rf_allocVector(STRSXP, n_slots));
  for (i=0; i<n_slots; i++)
    SET_STRING_ELT(sexp_data_names, i, Rf_mkChar(slot_name[i]));

  /* create data list */
  PROTECT(sexp_data = Rf_allocVector(VECSXP, n_slots));
  SET_VECTOR_ELT(sexp_data, pmid,   sexp_mid);      /* attach 'mid' element   */
  SET_VECTOR_ELT(sexp_data, porder, sexp_order);    /* attach 'order' element */
  SET_VECTOR_ELT(sexp_data, pUmax,  sexp_Umax);     /* attach 'Umax' element  */
  SET_VECTOR_ELT(sexp_data, pguide, sexp_guide);    /* attach 'guide' element */
  SET_VECTOR_ELT(sexp_data, piv,    sexp_iv);       /* attach 'iv' element    */

  /* attach vector names */
  Rf_setAttrib(sexp_data, R_NamesSymbol, sexp_data_names);

  /* store in slot 'data' of S4 object 'unur' */
  R_do_slot_assign(sexp_unur, Rf_install("data"), sexp_data);

  /* set domain of distribution and store in slot 'dom' */
  PROTECT(sexp_dom = Rf_allocVector(REALSXP, 2));
  REAL(sexp_dom)[0] = DISTR.domain[0];
  REAL(sexp_dom)[1] = DISTR.domain[1];
  R_do_slot_assign(sexp_unur, Rf_install("dom"), sexp_dom);

  /* o.k. */
  UNPROTECT(8);
  return;

} /* end of _Runuran_pack_pinv() */

/*---------------------------------------------------------------------------*/

SEXP
_Runuran_sample_pinv (SEXP sexp_data, int n)
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
  int i;
  double U;
  SEXP sexp_res = R_NilValue;

  /* extract data */
  int order = INTEGER(VECTOR_ELT(sexp_data, porder))[0];
  double Umax = REAL(VECTOR_ELT(sexp_data, pUmax))[0];
  int *guide = INTEGER(VECTOR_ELT(sexp_data, pguide));
  int guide_size = Rf_length(VECTOR_ELT(sexp_data, pguide));
  double *iv = REAL(VECTOR_ELT(sexp_data, piv));
  
  /* generate sample */
  PROTECT(sexp_res = Rf_allocVector(REALSXP, n));
  for (i=0; i<n; i++) {
    U = unif_rand();   /* FIXME: R built-in URNG hard coded ! */
    REAL(sexp_res)[i] = 
      _pinv_eval (U, Umax, order, guide_size, guide, iv);
  }

  /* return result to R */
  UNPROTECT(1);
  return sexp_res;
 
} /* end of _Runuran_sample_pinv() */

/*---------------------------------------------------------------------------*/

SEXP
_Runuran_quantile_pinv (SEXP sexp_data, SEXP sexp_U, SEXP sexp_unur)
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
  int i,n;
  double *U;
  SEXP sexp_res = R_NilValue;

  /* domain of distribution */
  SEXP sexp_dom  = R_NilValue;

  /* extract data */
  int order = INTEGER(VECTOR_ELT(sexp_data, porder))[0];
  double Umax = REAL(VECTOR_ELT(sexp_data, pUmax))[0];
  int *guide = INTEGER(VECTOR_ELT(sexp_data, pguide));
  int guide_size = Rf_length(VECTOR_ELT(sexp_data, pguide));
  double *iv = REAL(VECTOR_ELT(sexp_data, piv));

  /* Extract U */
  U = REAL(sexp_U);
  n = Rf_length(sexp_U);

  /* domain of distribution */
  PROTECT(sexp_dom = R_do_slot(sexp_unur, Rf_install("dom")));
  
  /* evaluate inverse CDF */
  PROTECT(sexp_res = Rf_allocVector(REALSXP, n));
  for (i=0; i<n; i++) {
    if (ISNAN(U[i]))
      /* if NA or NaN is given then we simply return the same value */
      REAL(sexp_res)[i] = U[i];
    else {

      if (U[i] <= 0. ||  U[i] >= 1.) {
	/* same bahavior as in UNU.RAN */

	if (U[i] < 0. ||  U[i] > 1.)
	  Rf_warning("[UNU.RAN - warning] argument out of domain: U not in [0,1]");
	if (U[i] < 0.5 )
	  REAL(sexp_res)[i] = REAL(sexp_dom)[0];
	if (U[i] > 0.5 )
	  REAL(sexp_res)[i] = REAL(sexp_dom)[1];
      }

      else {
	REAL(sexp_res)[i] = 
	  _pinv_eval (U[i], Umax, order, guide_size, guide, iv);
      }
    }
  }

  /* return result to R */
  UNPROTECT(2);
  return sexp_res;
 
} /* end of _Runuran_quantile_pinv() */

/*---------------------------------------------------------------------------*/

double
_pinv_eval (double U, double Umax, int order, int guide_size, int *guide, double *iv)
     /*----------------------------------------------------------------------*/
     /* Evaluate approximating polynomial.                                   */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   U          ... u-value ~ U(0,1)                                    */
     /*   Umax       ... upper bound of computational  domain for U          */
     /*   order      ... order of Newton polynomial                          */
     /*   guide_size ... size of guide table                                 */
     /*   guide      ... guide table                                         */
     /*   iv         ... array of coefficients for all Newton polynomials    */
     /*                                                                      */
     /* Return:                                                              */
     /*   (approximate) quantiles for given 'U' values                       */
     /*----------------------------------------------------------------------*/
{
  int I;
  double V,X;
  const double *c;
  int k, width;

  /* number of entries per interval */
  width = 2*order + 1;
  
  /* transform U ~ U(0,1) to V ~ U(0,Umax) */
  V = Umax * U;

  /* find interval */
  I = guide[(int) (U * guide_size)];
  while (V > iv[I+width]) I+=width;

  /* compute interpolating polynomial for corresponding interval */
  V -= iv[I];
  c = iv+I+1;

  X = (c[0]*(V-c[1])+c[2])*(V-c[3])+c[4];
  for (k=3; k<order; k++)
    X = X*(V-c[2*k-1])+c[2*k];
  X = V*X+c[2*k-1];

  /* return result */
  return X;

} /* end of _pinv_eval() */

/*---------------------------------------------------------------------------*/
