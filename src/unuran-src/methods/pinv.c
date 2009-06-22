/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include <tests/unuran_tests.h>
#include <utils/lobatto_source.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "pinv.h"
#include "pinv_struct.h"
#define MAX_ORDER   (12)
#define PINV_UERROR_CORRECTION  (0.9)
#define PINV_DEFAULT_MAX_IVS  (10000)
#define PINV_MAX_LOBATTO_IVS  (20001)
#define PINV_PDFLLIM    (1.e-13)
#define PINV_UERROR_AREA_APPROX  (1.e-5)
#define PINV_TAILCUTOFF_FACTOR   (0.05)
#define PINV_TAILCUTOFF_MAX      (1.e-10) 
#define PINV_UTOL_CORRECTION  (0.05)
#define PINV_MAX_ITER_IVS    (10 * GEN->max_ivs)
#define PINV_GUIDE_FACTOR  (1)
#define PINV_VARIANT_PDF      0x0010u   
#define PINV_VARIANT_CDF      0x0020u   
#define PINV_DEBUG_REINIT    0x00000002u   
#define PINV_DEBUG_TABLE     0x00000010u   
#define PINV_DEBUG_SEARCHBD  0x00010000u   
#define PINV_DEBUG_ITABLE    0x00020000u   
#define PINV_SET_ORDER          0x001u  
#define PINV_SET_U_RESOLUTION   0x002u  
#define PINV_SET_BOUNDARY       0x008u  
#define PINV_SET_SEARCHBOUNDARY 0x010u  
#define PINV_SET_VARIANT        0x020u  
#define PINV_SET_MAX_IVS        0x040u  
#define GENTYPE "PINV"         
static struct unur_gen *_unur_pinv_init (struct unur_par *par);
static struct unur_gen *_unur_pinv_create (struct unur_par *par);
static int _unur_pinv_check_par (struct unur_gen *gen);
static struct unur_gen *_unur_pinv_clone (const struct unur_gen *gen);
static void _unur_pinv_free (struct unur_gen *gen);
static int _unur_pinv_make_guide_table (struct unur_gen *gen);
static double _unur_pinv_eval_PDF (double x, struct unur_gen *gen);
static double _unur_pinv_sample (struct unur_gen *gen);
static double _unur_pinv_eval_approxinvcdf (const struct unur_gen *gen, double u);
static int _unur_pinv_preprocessing (struct unur_gen *gen);
static int _unur_pinv_relevant_support (struct unur_gen *gen);
static double _unur_pinv_searchborder (struct unur_gen *gen, double x0, double bound,
				       double *dom, int *search);
static int _unur_pinv_approx_pdfarea (struct unur_gen *gen);
static int _unur_pinv_pdfarea (struct unur_gen *gen);
static int _unur_pinv_computational_domain (struct unur_gen *gen);
static double _unur_pinv_cut (struct unur_gen *gen, double dom, double w, double dw, double crit);
static int _unur_pinv_computational_domain_CDF (struct unur_gen *gen);
static double _unur_pinv_cut_CDF( struct unur_gen *gen, double dom, double x0, double ul, double uu );
static double _unur_pinv_Udiff (struct unur_gen *gen, double x, double h, double *fx);
static int _unur_pinv_create_table( struct unur_gen *gen );
static int _unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx );
static int _unur_pinv_lastinterval( struct unur_gen *gen );
static int _unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
				     double *xval);
static int _unur_pinv_linear_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
				     double *xval);
static int _unur_pinv_chebyshev_points (int order, double *pt);
static double _unur_pinv_newton_eval (double q, double ui[], double zi[], int order);
static double _unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv, 
					  double xval[]);
static int _unur_pinv_newton_testpoints (int g, double ui[], double utest[]);
static int _unur_pinv_linear_testpoints (int g, double ui[], double utest[]);
#ifdef UNUR_ENABLE_LOGGING
static void _unur_pinv_debug_init_start (const struct unur_gen *gen);
static void _unur_pinv_debug_init (const struct unur_gen *gen, int ok);
static void _unur_pinv_debug_relevant_support (const struct unur_gen *gen);
static void _unur_pinv_debug_pdfarea (const struct unur_gen *gen, int approx);
static void _unur_pinv_debug_computational_domain (const struct unur_gen *gen);
static void _unur_pinv_debug_intervals (const struct unur_gen *gen );
static void _unur_pinv_debug_create_table (const struct unur_gen *gen,
					   int iter, int n_incr_h, int n_decr_h,
					   int n_use_linear);
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_pinv_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_pinv_par*)par->datap) 
#define GEN       ((struct unur_pinv_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define SAMPLE    gen->sample.cont      
#define PDF(x)  (_unur_pinv_eval_PDF((x),(gen)))      
#define CDF(x)  (_unur_cont_CDF((x),(gen->distr)))    
#define _unur_pinv_getSAMPLE(gen)  (_unur_pinv_sample)
#include "pinv_newset.ch"
#include "pinv_init.ch"
#include "pinv_sample.ch"
#include "pinv_prep.ch"
#include "pinv_newton.ch"
#include "pinv_debug.ch"
#include "pinv_info.ch"
