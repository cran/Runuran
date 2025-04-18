/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mixture.c                                                         *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         Create UNU.RAN generator object for mixture of distribution.      *
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

/*****************************************************************************/

SEXP
Runuran_mixt (SEXP sexp_obj, SEXP sexp_prob, SEXP sexp_comp, SEXP sexp_inversion)
/*---------------------------------------------------------------------------*/
/* Create UNU.RAN generator object for mixture of distribution.              */
/*                                                                           */
/* Parameters:                                                               */
/*   obj       ... S4 class that contains 'Runuran' generator object         */ 
/*   prob      ... probability vector                                        */ 
/*   comp      ... vector of components: 'Runuran' objects (S4 classes)      */ 
/*   inversion ... whether inversion method should be used (boolean)         */
/*                                                                           */
/* Return:                                                                   */
/*   pointer to UNU.RAN generator object                                     */
/*---------------------------------------------------------------------------*/
{
  struct unur_gen **comp;      /* pointer to list of UNU.RAN objects */
  double *prob;                /* probality vector */
  int n_comp;                  /* number of components */
  int useinversion;            /* whether inversion method should be used */
  
  SEXP sexp_unur;              /* pointer to element in R list 'comp' */
  SEXP sexp_gen;               /* R pointer to generator object */
  SEXP sexp_is_inversion;      /* Slot 'inversion' */

  struct unur_par *par;
  struct unur_gen *gen = NULL; /* pointer to UNU.RAN object */
  
  int i;

  /* extract boolean */
  useinversion = *LOGICAL(Rf_coerceVector(sexp_inversion, LGLSXP));

  /* extract length of component vector */
  n_comp = Rf_length(sexp_comp);
  if(n_comp != Rf_length(sexp_prob)) {
    Rf_errorcall(R_NilValue,"[UNU.RAN - error] 'prob' and 'comp' must have same length");
  }

  /* extract components */
  if (! Rf_isVectorList(sexp_comp)) {
    Rf_errorcall(R_NilValue,"[UNU.RAN - error] invalid argument 'comp'");
  }
  comp = (struct unur_gen**) R_alloc(n_comp, sizeof(struct unur_gen *) );

  for (i=0; i<n_comp; i++) {
    sexp_unur = VECTOR_ELT(sexp_comp, i);
    if (! Rf_isS4(sexp_unur)) {
      Rf_error("[UNU.RAN - error] argument 'comp[%d]' does not contain UNU.RAN objects",i+1);
    }
    sexp_gen = R_do_slot(sexp_unur, Rf_install("unur"));
    CHECK_UNUR_PTR(sexp_gen);
    if (Rf_isNull(sexp_gen) || 
	((comp[i]=R_ExternalPtrAddr(sexp_gen)) == NULL) ) {
      Rf_error("[UNU.RAN - error] invalid argument 'comp[%d]'. maybe packed?",i+1);
    }
  }

  /* extract probability vector */
  PROTECT(sexp_prob = Rf_coerceVector(sexp_prob, REALSXP));
  prob = REAL(sexp_prob);

  /* create UNU.RAN generator for mixture */
  if (! ISNAN(prob[0])) {
    par = unur_mixt_new(n_comp, prob, comp);
    if (useinversion) {
      unur_mixt_set_useinversion(par,TRUE);
    }
    gen = unur_init(par);
  }
  /* we do not need 'sexp_prob' any more */
  UNPROTECT(1);

  /* check that 'prob' had contained useful data */
  if (ISNAN(prob[0])) {
    Rf_errorcall(R_NilValue,"[UNU.RAN - error] invalid argument 'prob'");
  }
  
  /* 'gen' must not be a NULL pointer */
  if (gen == NULL) {
    Rf_errorcall(R_NilValue,"[UNU.RAN - error] cannot create UNU.RAN object");
  }

  /* set slot 'inversion' to true when 'gen' implements an inversion method. */
  PROTECT(sexp_is_inversion = Rf_allocVector(LGLSXP, 1));
  LOGICAL(sexp_is_inversion)[0] = useinversion;
  R_do_slot_assign(sexp_obj, Rf_install("inversion"), sexp_is_inversion);
    
  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_gen = R_MakeExternalPtr(gen, _Runuran_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_gen, _Runuran_free);

  /* return pointer to R */
  UNPROTECT(2);
  return (sexp_gen);

} /* end of Runuran_mixture() */

/*---------------------------------------------------------------------------*/
