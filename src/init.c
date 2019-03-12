/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: init.c                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         Initialize package 'Runuran'                                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "Runuran.h"
#include "Runuran_ext.h"
#include <time.h>

/* internal header files for UNU.RAN */
#include <unur_source.h>

/*---------------------------------------------------------------------------*/

/* List of functions to be registered as native routines */
static const R_CallMethodDef CallEntries[] = {
    {"Runuran_CDF",            (DL_FUNC) &Runuran_CDF,            2},
    {"Runuran_PDF",            (DL_FUNC) &Runuran_PDF,            3},
    {"Runuran_cmv_init",       (DL_FUNC) &Runuran_cmv_init,       9},
    {"Runuran_cont_init",      (DL_FUNC) &Runuran_cont_init,     11},
    {"Runuran_discr_init",     (DL_FUNC) &Runuran_discr_init,     9},
    {"Runuran_init",           (DL_FUNC) &Runuran_init,           3},
    {"Runuran_mixt",           (DL_FUNC) &Runuran_mixt,           4},
    {"Runuran_pack",           (DL_FUNC) &Runuran_pack,           1},
    {"Runuran_performance",    (DL_FUNC) &Runuran_performance,    2},
    {"Runuran_print",          (DL_FUNC) &Runuran_print,          2},
    {"Runuran_quantile",       (DL_FUNC) &Runuran_quantile,       2},
    {"Runuran_sample",         (DL_FUNC) &Runuran_sample,         2},
    {"Runuran_set_aux_seed",   (DL_FUNC) &Runuran_set_aux_seed,   1},
    {"Runuran_std_cont",       (DL_FUNC) &Runuran_std_cont,       4},
    {"Runuran_std_discr",      (DL_FUNC) &Runuran_std_discr,      4},
    {"Runuran_use_aux_urng",   (DL_FUNC) &Runuran_use_aux_urng,   2},
    {"Runuran_verify_hat",     (DL_FUNC) &Runuran_verify_hat,     2},
    {"Runuran_set_error_level",(DL_FUNC) &Runuran_set_error_level,1},
    {NULL, NULL, 0}
};

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
  unur_set_error_handler( _Runuran_error_handler_warning );

  /* Set R built-in generator as default URNG */
  unur_set_default_urng( unur_urng_new( _Runuran_R_unif_rand, NULL) );

  /* We use a built-in generator from the UNU.RAN library for the auxiliary URNG */
  {
    UNUR_URNG *aux;
    /* create URNG object */
    aux = unur_urng_new( unur_urng_MRG31k3p, NULL );
    unur_urng_set_reset( aux, unur_urng_MRG31k3p_reset );
    unur_urng_set_seed( aux, unur_urng_MRG31k3p_seed);
    /* seed URNG object */
    unur_urng_seed (aux, (long) time(NULL));
    /* set as auxiliary generator */
    unur_set_default_urng_aux( aux );
  }

  /* Register native routines */ 
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE); 
  R_forceSymbols(info, TRUE);
 
  /* Declare some C routines to be callable from other packages */ 

  /* For project 'RunuranTEMPL': */
  R_RegisterCCallable("Runuran", "cont_init",   (DL_FUNC) Runuran_ext_cont_init);
  R_RegisterCCallable("Runuran", "cont_params", (DL_FUNC) unur_distr_cont_get_pdfparams);

  /* For project 'rvgtdist': */

#define RREGDEF(name)  R_RegisterCCallable("Runuran", #name, (DL_FUNC) name)

  RREGDEF(unur_init);
  RREGDEF(unur_free);
  RREGDEF(unur_sample_cont);

  RREGDEF(unur_distr_free);

  RREGDEF(unur_urng_new);
  RREGDEF(unur_urng_free);

  RREGDEF(unur_set_default_debug);
  RREGDEF(unur_set_default_urng);
  RREGDEF(unur_set_default_urng_aux);
  RREGDEF(unur_get_default_urng);

  RREGDEF(unur_get_strerror);
  RREGDEF(unur_set_error_handler);

  RREGDEF(unur_distr_gig);
  
  RREGDEF(unur_arou_new);
  RREGDEF(unur_arou_get_sqhratio);

  RREGDEF(unur_ars_new);

  RREGDEF(unur_tabl_new);
  RREGDEF(unur_tabl_get_sqhratio);
  RREGDEF(unur_tabl_set_max_sqhratio);
  RREGDEF(unur_tabl_set_max_intervals);
  RREGDEF(unur_tabl_set_boundary);

  RREGDEF(unur_tdr_new);
  RREGDEF(unur_tdr_set_variant_ia);
  RREGDEF(unur_tdr_get_sqhratio);

  RREGDEF(unur_pinv_new);
  
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

