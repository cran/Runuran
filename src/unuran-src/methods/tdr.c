/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "tdr.h"
#include "tdr_struct.h"
#define TDR_VARMASK_T          0x000fu   
#define TDR_VAR_T_SQRT         0x0001u   
#define TDR_VAR_T_LOG          0x0002u   
#define TDR_VAR_T_POW          0x0003u   
#define TDR_VARMASK_VARIANT    0x00f0u   
#define TDR_VARIANT_GW         0x0010u   
#define TDR_VARIANT_PS         0x0020u   
#define TDR_VARIANT_IA         0x0030u   
#define TDR_VARFLAG_VERIFY     0x0100u   
#define TDR_VARFLAG_USECENTER  0x0200u   
#define TDR_VARFLAG_USEMODE    0x0400u   
#define TDR_VARFLAG_PEDANTIC   0x0800u   
#define TDR_VARFLAG_USEDARS    0x1000u   
#define TDR_DEBUG_IV           0x00000010u
#define TDR_DEBUG_SPLIT        0x00010000u
#define TDR_DEBUG_DARS         0x00020000u
#define TDR_DEBUG_SAMPLE       0x01000000u
#define TDR_DEBUG_REINIT       0x00000020u  
#define TDR_SET_CENTER         0x0002u
#define TDR_SET_STP            0x0001u
#define TDR_SET_N_STP          0x0002u
#define TDR_SET_PERCENTILES    0x0004u
#define TDR_SET_N_PERCENTILES  0x0008u
#define TDR_SET_RETRY_NCPOINTS 0x0010u
#define TDR_SET_GUIDEFACTOR    0x0020u
#define TDR_SET_C              0x0040u
#define TDR_SET_MAX_SQHRATIO   0x0080u
#define TDR_SET_MAX_IVS        0x0100u
#define TDR_SET_USE_DARS       0x0200u
#define TDR_SET_DARS_FACTOR    0x0400u
#define GENTYPE "TDR"          
static struct unur_gen *_unur_tdr_init( struct unur_par *par );
static int _unur_tdr_reinit( struct unur_gen *gen );
static int _unur_tdr_make_gen( struct unur_gen *gen );
static struct unur_gen *_unur_tdr_create( struct unur_par *par );
static double _unur_tdr_gw_sample( struct unur_gen *generator );
static double _unur_tdr_gw_sample_check( struct unur_gen *generator );
static double _unur_tdr_ps_sample( struct unur_gen *generator );
static double _unur_tdr_ps_sample_check( struct unur_gen *generator );
static double _unur_tdr_ia_sample( struct unur_gen *generator );
static double _unur_tdr_ia_sample_check( struct unur_gen *generator );
static double _unur_tdr_gw_eval_invcdfhat( const struct unur_gen *generator, double u,
					   double *hx, double *fx, double *sqx,
					   struct unur_tdr_interval **iv,
					   struct unur_tdr_interval **cpt );
static double _unur_tdr_ps_eval_invcdfhat( const struct unur_gen *generator, double u,
					   double *hx, double *fx, double *sqx,
					   struct unur_tdr_interval **iv );
static void _unur_tdr_free( struct unur_gen *gen);
static struct unur_gen *_unur_tdr_clone( const struct unur_gen *gen );
static int _unur_tdr_starting_cpoints( struct unur_gen *gen );
static int _unur_tdr_starting_intervals( struct unur_gen *gen );
static int _unur_tdr_gw_starting_intervals( struct unur_gen *gen );
static int _unur_tdr_ps_starting_intervals( struct unur_gen *gen );
static int _unur_tdr_run_dars( struct unur_gen *gen );
static int _unur_tdr_gw_dars( struct unur_gen *gen );
static int _unur_tdr_ps_dars( struct unur_gen *gen );
static int _unur_tdr_gw_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv );
static int _unur_tdr_ps_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv );
static struct unur_tdr_interval *_unur_tdr_interval_new( struct unur_gen *gen, 
							 double x, double fx, int is_mode );
static int _unur_tdr_tangent_intersection_point( struct unur_gen *gen,
						 struct unur_tdr_interval *iv, double *ipt );
static double _unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv,
				       double slope, double x );
static double _unur_tdr_interval_xxarea( struct unur_gen *gen, struct unur_tdr_interval *iv,
					 double slope, double x );
static double _unur_tdr_eval_intervalhat( struct unur_gen *gen,
					  struct unur_tdr_interval *iv, double x );
static double _unur_tdr_eval_cdfhat( struct unur_gen *gen, double x );
static int _unur_tdr_gw_interval_split( struct unur_gen *gen, 
					struct unur_tdr_interval *iv_old, double x, double fx );
static int _unur_tdr_ps_interval_split( struct unur_gen *gen, 
					struct unur_tdr_interval *iv_old, double x, double fx );
static int _unur_tdr_gw_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv, 
				     double x, double fx);
static int _unur_tdr_ps_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv, 
				     double x, double fx);
static int _unur_tdr_make_guide_table( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_tdr_debug_init_start( const struct unur_gen *gen );
static void _unur_tdr_debug_init_finished( const struct unur_gen *gen );
static void _unur_tdr_debug_dars_start( const struct unur_gen *gen );
static void _unur_tdr_debug_dars_finished( const struct unur_gen *gen );
static void _unur_tdr_debug_reinit_start( const struct unur_gen *gen );
static void _unur_tdr_debug_reinit_retry( const struct unur_gen *gen );
static void _unur_tdr_debug_reinit_finished( const struct unur_gen *gen );
static void _unur_tdr_debug_free( const struct unur_gen *gen );
static void _unur_tdr_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
static void _unur_tdr_gw_debug_intervals( const struct unur_gen *gen, int print_areas );
static void _unur_tdr_ps_debug_intervals( const struct unur_gen *gen, int print_areas );
static void _unur_tdr_gw_debug_sample( const struct unur_gen *gen,
				       const struct unur_tdr_interval *iv,
				       const struct unur_tdr_interval *pt, 
				       double x, double fx, double hx, double sqx );
static void _unur_tdr_ps_debug_sample( const struct unur_gen *gen, 
				       const struct unur_tdr_interval *iv,
				       double x, double fx, double hx, double sqx );
static void _unur_tdr_gw_debug_split_start( const struct unur_gen *gen, 
					    const struct unur_tdr_interval *iv,
					    double x, double fx );
static void _unur_tdr_gw_debug_split_stop( const struct unur_gen *gen, 
					   const struct unur_tdr_interval *iv_left,
					   const struct unur_tdr_interval *iv_right );
static void _unur_tdr_ps_debug_split_start( const struct unur_gen *gen, 
					    const struct unur_tdr_interval *iv_left,
					    const struct unur_tdr_interval *iv_right,
					    double x, double fx );
static void _unur_tdr_ps_debug_split_stop( const struct unur_gen *gen, 
					   const struct unur_tdr_interval *iv_left,
					   const struct unur_tdr_interval *iv_middle,
					   const struct unur_tdr_interval *iv_right );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_tdr_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cont      
#define PAR       ((struct unur_tdr_par*)par->datap) 
#define GEN       ((struct unur_tdr_gen*)gen->datap) 
#define DISTR     gen->distr->data.cont 
#define BD_LEFT   domain[0]             
#define BD_RIGHT  domain[1]             
#define SAMPLE    gen->sample.cont           
#define PDF(x)     _unur_cont_PDF((x),(gen->distr))      
#define dPDF(x)    _unur_cont_dPDF((x),(gen->distr))     
#define logPDF(x)  _unur_cont_logPDF((x),(gen->distr))   
#define dlogPDF(x) _unur_cont_dlogPDF((x),(gen->distr))  
static UNUR_SAMPLING_ROUTINE_CONT *
_unur_tdr_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    
    return (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_gw_sample_check : _unur_tdr_gw_sample;
  case TDR_VARIANT_IA:    
    return (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ia_sample_check : _unur_tdr_ia_sample;
  case TDR_VARIANT_PS:    
  default:
    return (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ps_sample_check : _unur_tdr_ps_sample;
  }
} 
#include "tdr_newset.ch"
#include "tdr_init.ch"
#include "tdr_sample.ch"
#include "tdr_debug.ch"
#include "tdr_info.ch"
