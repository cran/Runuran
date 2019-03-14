/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: performance.c                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         Get information about performance of UNU.RAN generator objects    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2011 Wolfgang Hoermann and Josef Leydold                  *
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
#include <distr/distr_source.h>

/* structures used by particular UNU.RAN methods */
#include <methods/arou_struct.h>
#include <methods/ars_struct.h>
#include <methods/cstd_struct.h>
#include <methods/dari_struct.h>
#include <methods/dsrou_struct.h>
#include <methods/dstd_struct.h>
#include <methods/hinv_struct.h>
#include <methods/itdr_struct.h>
#include <methods/nrou_struct.h>
#include <methods/pinv_struct.h>
#include <methods/srou_struct.h>
#include <methods/tabl_struct.h>
#include <methods/tdr_struct.h>


/*****************************************************************************/
/* array for storing list elements                                           */

#define MAX_LIST  (11)       /* maximum number of list entries */

struct Rlist {
  int len;                   /* length of list (depends on method) */
  char *names[MAX_LIST];     /* names of list elements */
  SEXP values;               /* pointer to R list */
};

/* functions for appending list elements */
static void add_string(struct Rlist *list, char *key, const char *string);
static void add_numeric(struct Rlist *list, char *key, double num);
static void add_numeric_vec(struct Rlist *list, char *key, double *num, int n_num);
static void add_integer(struct Rlist *list, char *key, int inum);
static void add_integer_vec(struct Rlist *list, char *key, int *inum, int n_num);

/*****************************************************************************/
/* add list elements                                                         */

void add_string(struct Rlist *list, char *key, const char *string)
{
  /* check length of list */
  if (list->len >= MAX_LIST)
    error("Runuran: Internal error! Please send bug report.");

  /* store name of list entry */
  list->names[list->len] = key;

  /* create R object for list entry */
  SET_VECTOR_ELT(list->values, list->len, mkString(string));
  
  /* update length of list */
  ++list->len;
} /* end of add_string() */

/* ------------------------------------------------------------------------- */

void add_numeric(struct Rlist *list, char *key, double num)
{
  if (list->len >= MAX_LIST)
    error("Runuran: Internal error! Please send bug report.");

  list->names[list->len] = key;
  SET_VECTOR_ELT(list->values, list->len, ScalarReal(num));
  ++list->len;
} /* end of add_numeric() */

/* ------------------------------------------------------------------------- */

void add_numeric_vec(struct Rlist *list, char *key, double *num, int n_num)
{
  int i;
  SEXP val;

  if (list->len >= MAX_LIST)
    error("Runuran: Internal error! Please send bug report.");

  list->names[list->len] = key;

  val = NEW_NUMERIC(n_num);
  for (i=0; i<n_num; i++)
    REAL(val)[i] = num[i];
  SET_VECTOR_ELT(list->values, list->len, val);

  ++list->len;
} /* end of add_numeric_list() */

/* ------------------------------------------------------------------------- */

void add_integer(struct Rlist *list, char *key, int inum)
{
  if (list->len >= MAX_LIST)
    error("Runuran: Internal error! Please send bug report.");

  list->names[list->len] = key;
  SET_VECTOR_ELT(list->values, list->len, ScalarInteger(inum));
  ++list->len;
} /* end of add_integer() */

/* ------------------------------------------------------------------------- */

void add_integer_vec(struct Rlist *list, char *key, int *inum, int n_num)
{
  int i;
  SEXP val;

  if (list->len >= MAX_LIST)
    error("Runuran: Internal error! Please send bug report.");

  list->names[list->len] = key;

  val = NEW_INTEGER(n_num);
  for (i=0; i<n_num; i++)
    INTEGER(val)[i] = inum[i];
  SET_VECTOR_ELT(list->values, list->len, val);

  ++list->len;
} /* end of add_integer_list() */

/* ------------------------------------------------------------------------- */

/*****************************************************************************/

SEXP
Runuran_performance (SEXP sexp_unur, SEXP sexp_debug)
     /*----------------------------------------------------------------------*/
     /* Get some informations about UNU.RAN generator object in an R list.   */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   unur  ... 'Runuran' object (S4 class)                              */ 
     /*   debug ... get more information if TRUE                             */ 
     /*                                                                      */
     /* Return:                                                              */
     /*   R list                                                             */
     /*----------------------------------------------------------------------*/
{
  SEXP sexp_list;              /* pointer to R list */
  SEXP sexp_names;             /* array of keywords */

  SEXP sexp_gen;               /* R pointer to generator object */
  struct unur_gen *gen = NULL; /* pointer to UNU.RAN object */
  SEXP sexp_data;              /* R pointer to data list in generator object 
				  (must be empty) */
  int debug;                   /* whether more information is provided */
  int i;                       /* aux loop variable */

  /* array of list elements */
  struct Rlist list;

  /* debug or not debug */
  debug = *(LOGICAL( AS_LOGICAL(sexp_debug) ));

  /* slot 'data' should not be pesent */
  sexp_data = GET_SLOT(sexp_unur, install("data"));
  if (! isNull(sexp_data)) {
    Rprintf("Object is PACKED !\n\n");
    return R_NilValue;
  }

  /* Extract pointer to UNU.RAN generator */
  sexp_gen = GET_SLOT(sexp_unur, install("unur"));
  CHECK_UNUR_PTR(sexp_gen);
  if (isNull(sexp_gen) || 
      ((gen=R_ExternalPtrAddr(sexp_gen)) == NULL) ) {
    warningcall_immediate(R_NilValue,"[UNU.RAN - warning] empty UNU.RAN object");
    return R_NilValue;
  }

  /* create temporary R list */
  PROTECT(list.values = allocVector(VECSXP, MAX_LIST));
  list.len = 0;

  /* we use macros to get an overview of used keywords and */
  /* to minimize the risk of typos.                        */

  /* set method */
#define METHOD(string)       add_string(&list,"method",(string))

  /* set kind (type) of generation method */
#define KIND_AR              add_string(&list,"type","ar")
#define KIND_INV             add_string(&list,"type","inv")
#define KIND_IAR             add_string(&list,"type","iar")
#define KIND_MCMC            add_string(&list,"type","mcmc")
#define KIND_OTHER           add_string(&list,"type","other")

  /* set class (type) of distribution */
  /* Remark we currently to not distinguish between "cont" and "cemp" */
#define CLASS_CEMP           add_string(&list,"distr.class","cont")
#define CLASS_CONT           add_string(&list,"distr.class","cont")
#define CLASS_CVEC           add_string(&list,"distr.class","cmv")
#define CLASS_DISCR          add_string(&list,"distr.class","discr")

  /* rejection method */
#define REJECTIONCONST(num)  add_numeric(&list,"rejection.constant",(num))
#define AREA_HAT(num)        add_numeric(&list,"area.hat",(num))
#define AREA_SQUEEZE(num)    add_numeric(&list,"area.squeeze",(num))
#define NINTS(inum)          add_integer(&list,"intervals",(inum))

  /* (approximate) inversion method */
#define TRUNC(left,right)    {				\
    double tmp[2];					\
    tmp[0]=(left); tmp[1]=(right);			\
    add_numeric_vec(&list,"truncated.domain",(tmp),2);	\
  } 
#define AREA_PDF(num)        add_numeric(&list,"area.pdf",(num))


  /* get data */
  switch (unur_get_method(gen)) {

    /************************/
    /* discrete, univariate */
    /************************/
    /* data for distribution in generator object */
#define DISTR  gen->distr->data.discr

    /* ..................................................................... */
  case UNUR_METH_DARI:
#define GEN ((struct unur_dari_gen*)gen->datap)
    METHOD("DARI"); KIND_AR; CLASS_DISCR;
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PMFSUM)
		   ? GEN->vt/DISTR.sum : NA_REAL);
    AREA_HAT (GEN->vt);
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_DAU:
    METHOD("DAU"); KIND_OTHER; CLASS_DISCR;
    break;
    /* ..................................................................... */
  case UNUR_METH_DGT:
    METHOD("DGT"); KIND_INV; CLASS_DISCR;
    break;
    /* ..................................................................... */
  case UNUR_METH_DSROU:
#define GEN ((struct unur_dsrou_gen*)gen->datap)
    METHOD("DSROU"); KIND_AR; CLASS_DISCR;
    REJECTIONCONST (2.*(-GEN->al+GEN->ar) / DISTR.sum);
    AREA_HAT (2.*(-GEN->al+GEN->ar));
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_DSS:
    METHOD("DSS"); KIND_INV; CLASS_DISCR;
    break;
    /* ..................................................................... */
  case UNUR_METH_DSTD:
#define GEN ((struct unur_dstd_gen*)gen->datap)
    METHOD("DSTD"); KIND_OTHER; CLASS_DISCR;
    if (debug) {
      add_numeric_vec(&list, "genparam", GEN->gen_param, GEN->n_gen_param);
      add_integer_vec(&list, "geniparam", GEN->gen_iparam, GEN->n_gen_iparam);
    }
#undef GEN
    break;
    /* ..................................................................... */
    
#undef DISTR


    /**************************/
    /* continuous, univariate */
    /**************************/

#define DISTR  gen->distr->data.cont 

    /* ..................................................................... */
  case UNUR_METH_AROU:
#define GEN ((struct unur_arou_gen*)gen->datap)
    METHOD("AROU"); KIND_AR; CLASS_CONT;
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PDFAREA) 
		   ? 2.*GEN->Atotal/DISTR.area : NA_REAL);
    AREA_HAT (2.*GEN->Atotal);
    AREA_SQUEEZE (2.*GEN->Asqueeze);
    NINTS (GEN->n_segs);
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_ARS:
#define GEN ((struct unur_ars_gen*)gen->datap)
    METHOD("ARS"); KIND_IAR; CLASS_CONT;
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PDFAREA) 
		   ? GEN->Atotal*exp(GEN->logAmax)/DISTR.area : NA_REAL);
    AREA_HAT (GEN->Atotal*exp(GEN->logAmax));
    NINTS (GEN->n_ivs);
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_CSTD:
#define GEN ((struct unur_cstd_gen*)gen->datap)
    METHOD("CSTD"); KIND_OTHER; CLASS_CONT;
    if (debug) {
      add_numeric_vec(&list, "genparam", GEN->gen_param, GEN->n_gen_param);
    }
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_HINV:
#define GEN ((struct unur_hinv_gen*)gen->datap)
    METHOD("HINV"); KIND_INV; CLASS_CONT;
    TRUNC(GEN->bleft,GEN->bright);
    NINTS (GEN->N-1);
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_HRB:
    METHOD("HRB"); KIND_OTHER; CLASS_CONT;
    break;
    /* ..................................................................... */
  case UNUR_METH_HRD:
    METHOD("HRD"); KIND_OTHER; CLASS_CONT;
    break;
    /* ..................................................................... */
  case UNUR_METH_HRI:
    METHOD("HRI"); KIND_OTHER; CLASS_CONT;
    break;
    /* ..................................................................... */
  case UNUR_METH_ITDR:
#define GEN ((struct unur_itdr_gen*)gen->datap)
    METHOD("ITDR"); KIND_AR; CLASS_CONT;
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PDFAREA) 
		    ? GEN->Atot/DISTR.area : NA_REAL);
    AREA_HAT (GEN->Atot);
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_NINV:
    METHOD("NINV"); KIND_INV; CLASS_CONT;
    break;
    /* ..................................................................... */
  case UNUR_METH_NROU:
#define GEN ((struct unur_nrou_gen*)gen->datap)
    METHOD("NROU"); KIND_AR; CLASS_CONT;
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PDFAREA) 
		   ? 2.*(GEN->umax - GEN->umin)*GEN->vmax / DISTR.area : NA_REAL);
    AREA_HAT (2.*(GEN->umax - GEN->umin) * GEN->vmax);
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_PINV:
#define GEN ((struct unur_pinv_gen*)gen->datap)
    METHOD("PINV"); KIND_INV; CLASS_CONT;
    TRUNC(GEN->bleft,GEN->bright);
    AREA_PDF(GEN->area);
    NINTS (GEN->n_ivs);

    if (debug) {
      int j,n;
      SEXP val;

      if (list.len+4 >= MAX_LIST)
	error("Runuran: Internal error! Please send bug report.");

      list.names[list.len] = "cdfi";
      val = NEW_NUMERIC(GEN->n_ivs + 1);
      for (n=0; n<=GEN->n_ivs; n++)
	REAL(val)[n] = GEN->iv[n].cdfi;
      SET_VECTOR_ELT(list.values, list.len, val);
      ++list.len;

      list.names[list.len] = "xi";
      val = NEW_NUMERIC(GEN->n_ivs + 1);
      for (n=0; n<=GEN->n_ivs; n++)
	REAL(val)[n] = GEN->iv[n].xi;
      SET_VECTOR_ELT(list.values, list.len, val);
      ++list.len;

      /* points for constructing Newton interpolation */
      list.names[list.len] = "ui";
      val = allocMatrix(REALSXP, GEN->n_ivs, GEN->order);
      for (n=0; n<GEN->n_ivs; n++) {
	for (j=0; j<GEN->order; j++)
	  REAL(val)[n + j*GEN->n_ivs] = GEN->iv[n].ui[j];
      }
      SET_VECTOR_ELT(list.values, list.len, val);
      ++list.len;

      /* coefficients for Newton interpolation */
      list.names[list.len] = "zi";
      val = allocMatrix(REALSXP, GEN->n_ivs, GEN->order);
      for (n=0; n<GEN->n_ivs; n++) {
	for (j=0; j<GEN->order; j++)
	  REAL(val)[n + j*GEN->n_ivs] = GEN->iv[n].zi[j];
      }
      SET_VECTOR_ELT(list.values, list.len, val);
      ++list.len;
    }
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_SROU:
#define GEN ((struct unur_srou_gen*)gen->datap)
#define SROU_VARFLAG_MIRROR   0x008u   /* use mirror principle */
    METHOD("SROU"); KIND_AR; CLASS_CONT;
    if (!_unur_isone(GEN->r)) {
      REJECTIONCONST (NA_REAL);
      AREA_HAT (NA_REAL);
    }
    else {
      REJECTIONCONST((GEN->Fmode >= 0.) ? 2. : ((gen->variant & SROU_VARFLAG_MIRROR) ? 2.829 : 4.));
      AREA_HAT(2.*(GEN->vr - GEN->vl) * GEN->um);
    }
#undef SROU_VARFLAG_MIRROR
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_SSR:
    METHOD("SSR"); KIND_AR; CLASS_CONT;
    REJECTIONCONST (NA_REAL);
    AREA_HAT (NA_REAL);
    break;
    /* ..................................................................... */
  case UNUR_METH_TABL:
#define GEN ((struct unur_tabl_gen*)gen->datap)
#define TABL_VARIANT_IA   0x0001u   /* use immediate acceptance */
    METHOD("TABL"); KIND_AR; CLASS_CONT;
    if (gen->variant&TABL_VARIANT_IA) {KIND_AR;} else {KIND_IAR;}
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PDFAREA) 
		    ? GEN->Atotal/DISTR.area : NA_REAL);
    AREA_HAT (GEN->Atotal);
    AREA_SQUEEZE (GEN->Asqueeze);
    NINTS (GEN->n_ivs);
#undef TABL_VARIANT_IA
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_TDR:
#define GEN ((struct unur_tdr_gen*)gen->datap)
#define TDR_VARIANT_IA       0x0030u   /* use immediate acceptance */
#define TDR_VARMASK_VARIANT  0x00f0u   /* indicates which variant  */
    METHOD("TDR");
    if ((gen->variant & TDR_VARMASK_VARIANT) == TDR_VARIANT_IA) {KIND_AR;} else {KIND_IAR;}
    CLASS_CONT;
    REJECTIONCONST ((gen->distr->set & UNUR_DISTR_SET_PDFAREA) 
		    ? GEN->Atotal/DISTR.area : NA_REAL);
    AREA_HAT (GEN->Atotal);
    AREA_SQUEEZE (GEN->Asqueeze);
    NINTS (GEN->n_ivs);
#undef TDR_VARMASK_VARIANT
#undef TDR_VARIANT_IA
#undef GEN
    break;
    /* ..................................................................... */
  case UNUR_METH_UTDR:
    METHOD("UTDR"); KIND_AR; CLASS_CONT;
    REJECTIONCONST (NA_REAL);
    AREA_HAT (NA_REAL);
    break;
    /* ..................................................................... */

  case UNUR_METH_MIXT:
#define GEN ((struct unur_mixt_gen*)gen->datap)
    METHOD("MIXT"); KIND_OTHER; CLASS_CONT;
#undef GEN
    break;
    /* ..................................................................... */

#undef DISTR


    /*************************/
    /* continuous, empirical */
    /*************************/

    /* ..................................................................... */
  case UNUR_METH_EMPK:
    METHOD("EMPK"); KIND_OTHER; CLASS_CEMP;
    break;
    /* ..................................................................... */
  case UNUR_METH_EMPL:
    METHOD("EMPL"); KIND_INV; CLASS_CEMP;
    break;
    /* ..................................................................... */
  case UNUR_METH_HIST:
    METHOD("HIST"); KIND_INV; CLASS_CEMP;
    break;
    /* ..................................................................... */
  case UNUR_METH_VEMPK:
    METHOD("VEMPK"); KIND_OTHER; CLASS_CEMP;
    break;
    /* ..................................................................... */


    /********************************************/
    /* continuous, multivariate (random vector) */
    /********************************************/

    /* ..................................................................... */
  case UNUR_METH_GIBBS:
    METHOD("GIBBS"); KIND_MCMC; CLASS_CVEC;
    break;
    /* ..................................................................... */
  case UNUR_METH_HITRO:
    METHOD("HITRO"); KIND_MCMC; CLASS_CVEC;
    break;
    /* ..................................................................... */
  case UNUR_METH_MVSTD:
    METHOD("MVSTD"); KIND_OTHER; CLASS_CVEC;
    break;
    /* ..................................................................... */
  case UNUR_METH_MVTDR:
    METHOD("MVTDR"); KIND_AR; CLASS_CVEC;
    break;
    /* ..................................................................... */
  case UNUR_METH_VNROU:
    METHOD("VNROU"); KIND_AR; CLASS_CVEC;
    break;
    /* ..................................................................... */

    /* case UNUR_METH_NORTA: */
    /* case UNUR_METH_MCORR: */


    /********/
    /* misc */
    /********/

    /* ..................................................................... */
  case UNUR_METH_UNIF:
    METHOD("UNIF"); KIND_INV; CLASS_CONT;
    break;
    /* ..................................................................... */

    /* case UNUR_METH_DEXT: */
    /* case UNUR_METH_CEXT: */

    /* automatic method --> meta method: selects one of the other methods automatically */
    /* case UNUR_METH_AUTO: */

    /* ..................................................................... */
  default: /* unknown ! */
    METHOD("NA"); 
  }

  /* create final list */
  PROTECT(sexp_list = allocVector(VECSXP, list.len)); 
  for(i = 0; i < list.len; i++)
    SET_VECTOR_ELT(sexp_list, i, VECTOR_ELT(list.values, i));
    
  /* an array of the "names" attribute of the objects in our list */
  PROTECT(sexp_names = allocVector(STRSXP, list.len));
  for(i = 0; i < list.len; i++)
    SET_STRING_ELT(sexp_names, i,  mkChar(list.names[i]));
 
  /* attach attribute names */
  setAttrib(sexp_list, R_NamesSymbol, sexp_names);

  /* return list */
  UNPROTECT(3);
  return sexp_list;

} /* end of Runuran_performance() */

/*---------------------------------------------------------------------------*/
