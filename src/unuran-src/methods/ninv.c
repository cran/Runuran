/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "ninv.h"
#include "ninv_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define INTERVAL_COVERS  (0.5)
#define NINV_VARFLAG_NEWTON   0x1u   
#define NINV_VARFLAG_REGULA   0x2u   
#define NINV_VARFLAG_BISECT   0x4u   
#define NINV_DEBUG_REINIT    0x00000002u   
#define NINV_DEBUG_TABLE     0x00000010u   
#define NINV_DEBUG_CHG       0x00001000u   
#define NINV_DEBUG_SAMPLE    0x01000000u   
#define NINV_SET_MAX_ITER     0x001u   
#define NINV_SET_X_RESOLUTION 0x002u   
#define NINV_SET_U_RESOLUTION 0x004u   
#define NINV_SET_START        0x008u   
#define GENTYPE "NINV"         
static struct unur_gen *_unur_ninv_init( struct unur_par *par );
static int _unur_ninv_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_ninv_create( struct unur_par *par );
static int _unur_ninv_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_ninv_clone( const struct unur_gen *gen );
static void _unur_ninv_free( struct unur_gen *gen );
static int _unur_ninv_create_table( struct unur_gen *gen );
static int _unur_ninv_compute_start( struct unur_gen *gen );
static double _unur_ninv_sample_newton( struct unur_gen *gen );
static double _unur_ninv_sample_regula( struct unur_gen *gen );
static double _unur_ninv_sample_bisect( struct unur_gen *gen );
static double _unur_ninv_newton( const struct unur_gen *gen, double u);
static double _unur_ninv_regula( const struct unur_gen *gen, double u );
static double _unur_ninv_bisect( const struct unur_gen *gen, double u );
static int _unur_ninv_bracket( const struct unur_gen *gen, double u, 
			       double *xl, double *fl, double *xu, double *fu );
static int _unur_ninv_accuracy( const struct unur_gen *gen,
				double x_resol, double u_resol,
				double x0, double f0, double x1, double f1 );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_ninv_debug_init( const struct unur_gen *gen );
static void _unur_ninv_debug_start( const struct unur_gen *gen );
static void _unur_ninv_debug_sample( const struct unur_gen *gen, 
				     double u, double x, double fx, int iter );
static void _unur_ninv_debug_chg_truncated( const struct unur_gen *gen);
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_ninv_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_ninv_par*)par->datap) 
#define GEN       ((struct unur_ninv_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define SAMPLE    gen->sample.cont      
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    
#define CDF(x)    _unur_cont_CDF((x),(gen->distr))    
static UNUR_SAMPLING_ROUTINE_CONT *
_unur_ninv_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    return _unur_ninv_sample_newton;
  case NINV_VARFLAG_BISECT:
    return _unur_ninv_sample_bisect;
  case NINV_VARFLAG_REGULA:
  default:
    return _unur_ninv_sample_regula;
  }
} 
#include "ninv_newset.ch"
#include "ninv_init.ch"
#include "ninv_sample.ch"
#include "ninv_newton.ch"
#include "ninv_regula.ch"
#include "ninv_debug.ch"
#include "ninv_info.ch"
