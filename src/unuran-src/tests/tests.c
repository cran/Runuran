/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"
#define TEST_TIMING_LOG_SAMPLESIZE 5    
#define TEST_COUNTER_SAMPLESIZE 100000  
#define TEST_SAMPLE_ROWS    10     
#define TEST_SAMPLE_COLS    10
#define TEST_CHI2_VERBOSE   1     
#define TEST_CHI2_INTERVALS 100   
static int _unur_print_method( struct unur_par *par, FILE *out );
void 
unur_run_tests( struct unur_par *par, unsigned tests, FILE *out )
{
  struct unur_gen *gen = NULL;
  struct unur_par *par_clone = NULL;
  double time_setup, time_sample;
  _unur_check_NULL("Tests",par,RETURN_VOID);
  if (!out) out = stdout;
  if (_unur_print_method(par,out)!=UNUR_SUCCESS)
    return;  
  par_clone = _unur_par_clone(par);
  if (tests & UNUR_TEST_TIME)
    gen = unur_test_timing(par,TEST_TIMING_LOG_SAMPLESIZE, &time_setup, &time_sample, TRUE, out);
  else
    gen = _unur_init(par);
  if (!gen) { _unur_par_free(par_clone); return; }
  if (tests & UNUR_TEST_N_URNG )
    unur_test_count_urn(gen,TEST_COUNTER_SAMPLESIZE, TRUE, out);
  if (tests & UNUR_TEST_N_PDF )
    unur_test_par_count_pdf(par_clone,TEST_COUNTER_SAMPLESIZE, TRUE, out);
  if (tests & UNUR_TEST_SAMPLE )
    unur_test_printsample(gen,TEST_SAMPLE_ROWS,TEST_SAMPLE_COLS, out);
  if (tests & UNUR_TEST_CHI2)
    unur_test_chi2(gen,TEST_CHI2_INTERVALS,0,0,TEST_CHI2_VERBOSE, out);
  _unur_free(gen);
  _unur_par_free(par_clone); 
  return;
} 
int 
_unur_print_method( struct unur_par *par, FILE *out )
{
  switch (par->distr->type) {
  case UNUR_DISTR_DISCR:
    fprintf(out,"\nTYPE:\t\tdiscrete univariate distribution\n");
    break;
  case UNUR_DISTR_CONT:
    fprintf(out,"\nTYPE:\t\tcontinuous univariate distribution\n");
    break;
  case UNUR_DISTR_CEMP:
    fprintf(out,"\nTYPE:\t\tcontinuous univariate empirical distribution\n");
    break;
  case UNUR_DISTR_CVEC:
    fprintf(out,"\nTYPE:\t\tcontinuous multivariate distribution\n");
    break;
  default: 
    _unur_error("Tests",UNUR_ERR_GENERIC,"type of method unknown!");
    return UNUR_ERR_GENERIC;
  }
  switch (par->method) {
  case UNUR_METH_AUTO:
    COOKIE_CHECK(par,CK_AUTO_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tautomatic selection (AUTO)\n");
    break;
  case UNUR_METH_DAU:
    COOKIE_CHECK(par,CK_DAU_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\talias and alias-urn method (DAU)\n");
    break;
  case UNUR_METH_DGT:
    COOKIE_CHECK(par,CK_DGT_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tindexed search -- guide table (DGT)\n");
    break;
  case UNUR_METH_DSROU:
    COOKIE_CHECK(par,CK_DSROU_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tdiscrete simple universal ratio-of-uniforms search (DSROU)\n");
    break;
  case UNUR_METH_DSS:
    COOKIE_CHECK(par,CK_DSS_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tsequential search (DSS)\n");
    break;
  case UNUR_METH_DSTD:
    COOKIE_CHECK(par,CK_DSTD_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tspecial (DSTD)\n");
    break;
  case UNUR_METH_DEXT:
    COOKIE_CHECK(par,CK_DEXT_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\texternal generator (DEXT)\n");
    break;
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tautomatic ratio-of-uniforms method (NINV)\n");
    break;
  case UNUR_METH_HINV:
    COOKIE_CHECK(par,CK_HINV_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tnumerical inversion of CDF by Hermite Interpolation (HINV)\n");
    break;
  case UNUR_METH_ITDR:
    COOKIE_CHECK(par,CK_ITDR_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tinverse transformed density rejection (ITDR)\n");
    break;
  case UNUR_METH_NINV:
    COOKIE_CHECK(par,CK_NINV_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tnumerical inversion of CDF (NINV)\n");
    break;
  case UNUR_METH_SROU:
    COOKIE_CHECK(par,CK_SROU_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tsimple universal ratio-of-uniforms method (SROU)\n");
    break;
  case UNUR_METH_NROU:
    COOKIE_CHECK(par,CK_NROU_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tnaive universal ratio-of-uniforms method (NROU)\n");
    break;
  case UNUR_METH_SSR:
    COOKIE_CHECK(par,CK_SSR_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tsimple transformed density rejection with universal bounds (SSR)\n");
    break;
  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\trejection from piecewise constant hat (TABL)\n");
    break;
  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\ttransformed density rejection (TDR)\n");
    break;
  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\ttransformed density rejection, 3-point method (UTDR)\n");
    break;
  case UNUR_METH_CSTD:
    COOKIE_CHECK(par,CK_CSTD_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tspecial (CSTD)\n");
    break;
  case UNUR_METH_CEXT:
    COOKIE_CHECK(par,CK_CEXT_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\texternal generator (CEXT)\n");
    break;
  case UNUR_METH_EMPK:
    COOKIE_CHECK(par,CK_EMPK_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tempirical distribution with kernel smoothing (EMPK)\n");
    break;
  case UNUR_METH_GIBBS:
    COOKIE_CHECK(par,CK_GIBBS_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tMarkov Chain - GIBBS sampler (GIBBS)\n");
    break;
  case UNUR_METH_HITRO:
    COOKIE_CHECK(par,CK_HITRO_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\thit&run ratio-of-uniforms (HITRO)\n");
    break;
  case UNUR_METH_MVSTD:
    COOKIE_CHECK(par,CK_MVSTD_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tspecial (MVSTD)\n");
    break;
  case UNUR_METH_MVTDR:
    COOKIE_CHECK(par,CK_MVTDR_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tmultivariate transformed density rejection (MVTDR)\n");
    break;
  case UNUR_METH_NORTA:
    COOKIE_CHECK(par,CK_NORTA_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tnormal to anything (NORTA)\n");
    break;
  case UNUR_METH_VMT:
    COOKIE_CHECK(par,CK_VMT_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tvector matrix transformation (VMT)\n");
    break;
  case UNUR_METH_VNROU:
    COOKIE_CHECK(par,CK_VNROU_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\tvector naive ratio-of-uniforms (VNROU)\n");
    break;
  case UNUR_METH_UNIF:
    COOKIE_CHECK(par,CK_UNIF_PAR,UNUR_ERR_COOKIE);
    fprintf(out,"METHOD:\t\twrapper for uniform (UNIF)\n");
    break;
  default: 
    _unur_error("Tests",UNUR_ERR_GENERIC,"method unknown!");
    return UNUR_ERR_GENERIC;
  }
  return UNUR_SUCCESS;
} 
