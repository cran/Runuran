/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <distr/cvec.h>
#include <distributions/unur_distributions.h>
#include <utils/fmax_source.h>
#include <utils/matrix_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "tdr.h"
#include "mvtdr.h"
#include "mvtdr_struct.h"
#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif
#define GUIDE_TABLE_SIZE    1
#define FIND_TP_TOL         0.001   
#define TOLERANCE           (1.e-8)
#define MVTDR_TDR_SQH_RATIO (0.95)
#define MVTDR_VARFLAG_VERIFY     0x01u   
#define MVTDR_DEBUG_VERTEX      0x00000010u   
#define MVTDR_DEBUG_CONE        0x00000020u   
#define MVTDR_SET_STEPSMIN        0x001u   
#define MVTDR_SET_MAXCONES        0x002u   
#define MVTDR_SET_BOUNDSPLITTING  0x004u   
#define GENTYPE "MVTDR"        
static struct unur_gen *_unur_mvtdr_init( struct unur_par *par );
static struct unur_gen *_unur_mvtdr_create( struct unur_par *par );
static int _unur_mvtdr_sample_cvec( struct unur_gen *gen, double *vec );
static void _unur_mvtdr_free( struct unur_gen *gen);
static struct unur_gen *_unur_mvtdr_clone( const struct unur_gen *gen );
static int _unur_mvtdr_simplex_sample( const struct unur_gen *gen, double *U );
static struct unur_gen *_unur_mvtdr_gammagen( struct unur_gen *gen, double alpha );
static int _unur_mvtdr_create_hat( struct unur_gen *gen );
static int _unur_mvtdr_initial_cones( struct unur_gen *gen );
static CONE *_unur_mvtdr_cone_new( struct unur_gen *gen );
static int _unur_mvtdr_cone_center( struct unur_gen *gen, CONE *c );
static int _unur_mvtdr_cone_params( struct unur_gen *gen, CONE *c );
static double _unur_mvtdr_cone_logH( struct unur_gen *gen, CONE *c );
static int _unur_mvtdr_cone_split( struct unur_gen *gen, CONE *c, int step );
static int _unur_mvtdr_triangulate( struct unur_gen *gen, int step, int all);
static int _unur_mvtdr_cone_height( struct unur_gen *gen, CONE *c );
static int _unur_mvtdr_max_gamma( struct unur_gen *gen );
static double _unur_mvtdr_tp_min_aux(double t, void *p);
static double _unur_mvtdr_tp_min( double t, void *p );
static int _unur_mvtdr_tp_find( struct unur_gen *gen, CONE *c );
static int _unur_mvtdr_tp_search( struct unur_gen *gen, TP_ARG *a );
static int _unur_mvtdr_tp_bracket( struct unur_gen *gen, TP_ARG *a );
static int _unur_mvtdr_initial_vertices( struct unur_gen *gen );
static VERTEX *_unur_mvtdr_vertex_new( struct unur_gen *gen );
static int _unur_mvtdr_number_vertices( struct unur_gen *gen, int level );
static VERTEX *_unur_mvtdr_vertex_on_edge( struct unur_gen *gen, VERTEX **vl );
static int _unur_mvtdr_etable_new( struct unur_gen *gen, int size );
static void _unur_mvtdr_etable_free( struct unur_gen *gen );
static VERTEX *_unur_mvtdr_etable_find_or_insert( struct unur_gen *gen, VERTEX **vidx );
static int _unur_mvtdr_make_guide_table( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_mvtdr_debug_init_start( const struct unur_gen *gen );
static void _unur_mvtdr_debug_init_finished( const struct unur_gen *gen, int successful );
static void _unur_mvtdr_debug_vertices( const struct unur_gen *gen );
static void _unur_mvtdr_debug_cones( const struct unur_gen *gen );
#endif
#ifdef UNUR_ENABLE_INFO
static void _unur_mvtdr_info( struct unur_gen *gen, int help );
#endif
#define DISTR_IN  distr->data.cvec      
#define PAR       ((struct unur_mvtdr_par*)par->datap) 
#define GEN       ((struct unur_mvtdr_gen*)gen->datap) 
#define DISTR     gen->distr->data.cvec 
#define SAMPLE    gen->sample.cvec           
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))         
#define logPDF(x) _unur_cvec_logPDF((x),(gen->distr))      
#define dPDF(r,x) _unur_cvec_dPDF((r),(x),(gen->distr))    
#define dlogPDF(r,x) _unur_cvec_dlogPDF((r),(x),(gen->distr))   
#define GEN_GAMMA  gen->gen_aux
#define T(x)       (log(x))       
#define T_deriv(x) (1./(x))       
#define T_inv(x)   (exp(x))       
#define TP_LEFT    1              
#define TP_MIDDLE  2              
#define TP_RIGHT   3              
#define TP_BRACKET 4              
#include "mvtdr_newset.ch"
#include "mvtdr_init.ch"
#include "mvtdr_sample.ch"
#include "mvtdr_debug.ch"
#include "mvtdr_info.ch"
