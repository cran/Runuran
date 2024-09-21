/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "tabl.h"
#include "tabl_struct.h"
#define TABL_DEFAULT_COMPUTATION_LIMIT  1.e20
#define TABL_N_RETRY_DARS  5
#define TABL_N_RUN_ARS     10
#define TABL_VARMASK_VARIANT      0x000fu   
#define TABL_VARIANT_IA           0x0001u   
#define TABL_VARIANT_FAST         0x0002u   
#define TABL_VARMASK_SPLIT        0x00f0u  
#define TABL_VARFLAG_SPLIT_POINT  0x0010u  
#define TABL_VARFLAG_SPLIT_MEAN   0x0020u  
#define TABL_VARFLAG_SPLIT_ARC    0x0040u  
#define TABL_VARFLAG_USEEAR       0x0100u  
#define TABL_VARFLAG_USEDARS      0x0200u  
#define TABL_VARFLAG_PEDANTIC     0x0400u  
#define TABL_VARFLAG_VERIFY       0x0800u  
#define TABL_DEBUG_IV        0x00000100u 
#define TABL_DEBUG_IV_START  0x00000200u 
#define TABL_DEBUG_EAR       0x00000400u 
#define TABL_DEBUG_DARS      0x00000800u 
#define TABL_SET_GUIDEFACTOR      0x0001u
#define TABL_SET_SLOPES           0x0004u
#define TABL_SET_AREAFRACTION     0x0008u
#define TABL_SET_MAX_IVS          0x0010u
#define TABL_SET_MAX_SQHRATIO     0x0020u
#define TABL_SET_N_STP            0x0040u
#define TABL_SET_STP              0x0080u
#define TABL_SET_BOUNDARY         0x0100u
#define TABL_SET_USE_EAR          0x0200u
#define TABL_SET_USE_DARS         0x0400u
#define TABL_SET_DARS_FACTOR      0x0800u
#define GENTYPE "TABL"         
static struct unur_gen *_unur_tabl_init( struct unur_par *par );
static struct unur_gen *_unur_tabl_create( struct unur_par *par );
static struct unur_gen *_unur_tabl_clone( const struct unur_gen *gen );
static void _unur_tabl_free( struct unur_gen *gen);
static double _unur_tabl_rh_sample( struct unur_gen *gen );
static double _unur_tabl_rh_sample_check( struct unur_gen *gen );
static double _unur_tabl_ia_sample( struct unur_gen *gen );
static double _unur_tabl_ia_sample_check( struct unur_gen *gen );
static int _unur_tabl_get_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen );
static int _unur_tabl_get_intervals_from_cpoints( struct unur_par *par, struct unur_gen *gen );
static int _unur_tabl_compute_intervals( struct unur_par *par, struct unur_gen *gen );
static struct unur_tabl_interval *
_unur_tabl_run_equalarearule( struct unur_par *par, struct unur_gen *gen, struct unur_tabl_interval *iv_slope );
static int _unur_tabl_run_dars( struct unur_gen *gen );
static int
_unur_tabl_split_interval( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			   double x, double fx, unsigned split_mode );
static int
_unur_tabl_improve_hat( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			double x, double fx );
static int _unur_tabl_make_guide_table( struct unur_gen *gen );
static double _unur_tabl_eval_cdfhat( struct unur_gen *gen, double x );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_tabl_debug_init_start( const struct unur_par *par, const struct unur_gen *gen );
static void _unur_tabl_debug_init_finished( const struct unur_gen *gen );
static void _unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen );
static void _unur_tabl_debug_free( const struct unur_gen *gen );
static void _unur_tabl_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_tabl_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_tabl_par*)par->datap) 
#define GEN       ((struct unur_tabl_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    
static UNUR_SAMPLING_ROUTINE_CONT *
_unur_tabl_getSAMPLE( struct unur_gen *gen )
{
  if (gen->variant & TABL_VARIANT_IA)
    return (gen->variant & TABL_VARFLAG_VERIFY) 
      ? _unur_tabl_ia_sample_check 
      : _unur_tabl_ia_sample;
  else
    return (gen->variant & TABL_VARFLAG_VERIFY) 
      ? _unur_tabl_rh_sample_check 
      : _unur_tabl_rh_sample;
} 
#include "tabl_newset.ch"
#include "tabl_init.ch"
#include "tabl_sample.ch"
#include "tabl_debug.ch"
#include "tabl_info.ch"
