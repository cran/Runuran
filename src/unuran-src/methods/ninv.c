/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
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
#define MAX_STEPS (100)
#define STEPFAC  (0.4)
#define I_CHANGE_TO_BISEC (50)
#define NINV_VARFLAG_NEWTON   0x1u   
#define NINV_VARFLAG_REGULA   0x2u   
#define NINV_DEBUG_REINIT    0x00000002u   
#define NINV_DEBUG_TABLE     0x00000010u   
#define NINV_DEBUG_CHG       0x00001000u   
#define NINV_DEBUG_SAMPLE    0x01000000u   
#define NINV_SET_MAX_ITER     0x001u   
#define NINV_SET_X_RESOLUTION 0x002u   
#define NINV_SET_START        0x004u   
#define GENTYPE "NINV"         
static struct unur_gen *_unur_ninv_init( struct unur_par *par );
static int _unur_ninv_reinit( struct unur_gen *gen );
static struct unur_gen *_unur_ninv_create( struct unur_par *par );
static int _unur_ninv_check_par( struct unur_gen *gen );
static struct unur_gen *_unur_ninv_clone( const struct unur_gen *gen );
static void _unur_ninv_free( struct unur_gen *gen );
static double _unur_ninv_sample_regula( struct unur_gen *gen );
static double _unur_ninv_sample_newton( struct unur_gen *gen );
static double _unur_ninv_regula( struct unur_gen *gen, double u );
static double _unur_ninv_newton( struct unur_gen *gen, double u);
static int _unur_ninv_compute_start( struct unur_gen *gen );
static int _unur_ninv_create_table( struct unur_gen *gen );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_ninv_debug_init( const struct unur_gen *gen );
static void _unur_ninv_debug_start( const struct unur_gen *gen );
static void _unur_ninv_debug_sample_regula( const struct unur_gen *gen, 
					    double u, double x, double fx, int iter );
static void _unur_ninv_debug_sample_newton( const struct unur_gen *gen, 
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
  case NINV_VARFLAG_REGULA:
  default:
    return _unur_ninv_sample_regula;
  }
} 
struct unur_par *
unur_ninv_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_ninv_par) );
  COOKIE_SET(par,CK_NINV_PAR);
  par->distr   = distr;            
  PAR->max_iter  = 40;             
  PAR->rel_x_resolution = 1.0e-8;  
  PAR->s[0]      = 0.0;     
  PAR->s[1]      = 0.0;     
  PAR->table_on  = FALSE;   
  par->method   = UNUR_METH_NINV;          
  par->variant  = NINV_VARFLAG_REGULA;     
  par->set      = 0u;                      
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_ninv_init;
  return par;
} 
int
unur_ninv_set_usenewton( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   
    return UNUR_ERR_DISTR_REQUIRED;
 }
  par->variant = NINV_VARFLAG_NEWTON;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_useregula( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  par->variant = NINV_VARFLAG_REGULA;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_max_iter( struct unur_par *par, int max_iter )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximal iterations");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_iter = max_iter;
  par->set |= NINV_SET_MAX_ITER;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_max_iter( struct unur_gen *gen, int max_iter )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (max_iter < 1) {
    _unur_warning(gen->genid, UNUR_ERR_PAR_SET, "maximal iterations");
    return UNUR_ERR_PAR_SET;
  }
  GEN->max_iter = max_iter;
  gen->set |= NINV_SET_MAX_ITER;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_x_resolution( struct unur_par *par, double x_resolution )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (x_resolution < DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"x resolution");
    return UNUR_ERR_PAR_SET;
  }
  PAR->rel_x_resolution = x_resolution;
  par->set |= NINV_SET_X_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_x_resolution( struct unur_gen *gen, double x_resolution )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (x_resolution < DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"x resolution");
    return UNUR_ERR_PAR_SET;
  }
  GEN->rel_x_resolution = x_resolution;
  gen->set |= NINV_SET_X_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_start( struct unur_par *par, double s1, double s2 )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if ( s1 <= s2 ) {
     PAR->s[0] = s1;
     PAR->s[1] = s2;
  }
  else {
     PAR->s[0] = s2;
     PAR->s[1] = s1;
  }
  par->set |= NINV_SET_START;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_start( struct unur_gen *gen, double s1, double s2 )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if ( s1 <= s2 ) {
     GEN->s[0] = s1;
     GEN->s[1] = s2;
  }
  else {
     GEN->s[0] = s2;
     GEN->s[1] = s1;
  }
  GEN->table_on = FALSE;
  _unur_ninv_compute_start(gen);
  gen->set |= NINV_SET_START;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_table( struct unur_par *par, int tbl_pnts )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  PAR->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  PAR->table_on = TRUE;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_table( struct unur_gen *gen, int tbl_pnts )
{
  int result;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  GEN->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  result = _unur_ninv_create_table(gen); 
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_CHG) 
    if (result==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
#endif
  return result;
} 
int 
unur_ninv_chg_truncated( struct unur_gen *gen, double left, double right )
{
  double Umin, Umax;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  Umin = (left > -INFINITY) ? CDF(left)  : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;
  if (Umin > Umax) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  if (_unur_FP_equal(Umin,Umax)) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return UNUR_ERR_DISTR_SET;
    }
  }
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = Umin;
  GEN->Umax = Umax;
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_CHG) 
    _unur_ninv_debug_chg_truncated( gen );
#endif
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_ninv_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);
  if (par->variant == NINV_VARFLAG_NEWTON && ! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   
  }
  gen = _unur_ninv_create(par);
  _unur_par_free(par);
  if (!gen) { return NULL; }
  if (_unur_ninv_check_par(gen) != UNUR_SUCCESS) {
    _unur_ninv_free(gen); return NULL;
  }
  if (GEN->table_on) {
    if (_unur_ninv_create_table(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }
  else {
    if (_unur_ninv_compute_start(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_ninv_debug_init(gen);
#endif
  return gen;
} 
int
_unur_ninv_reinit( struct unur_gen *gen )
{
  int rcode;
  if ( (rcode = _unur_ninv_check_par(gen)) != UNUR_SUCCESS)
    return rcode;
  if (DISTR.upd_area != NULL)
    if ((DISTR.upd_area)(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cannot compute normalization constant");
      return UNUR_ERR_GEN_DATA;
    }
  if (GEN->table != NULL)
    rcode = _unur_ninv_create_table(gen);
  else 
    rcode = unur_ninv_chg_start( gen, 0., 0. );
  SAMPLE = _unur_ninv_getSAMPLE(gen);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_REINIT) {
    _unur_distr_cont_debug( gen->distr, gen->genid );
    if (rcode==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
  }
#endif
  return UNUR_SUCCESS;
} 
static struct unur_gen *
_unur_ninv_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NINV_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_ninv_gen) );
  COOKIE_SET(gen,CK_NINV_GEN);
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_ninv_getSAMPLE(gen);
  gen->destroy = _unur_ninv_free;
  gen->clone = _unur_ninv_clone;
  gen->reinit = _unur_ninv_reinit;
  GEN->max_iter = PAR->max_iter;      
  GEN->rel_x_resolution = PAR->rel_x_resolution; 
  GEN->table_on = PAR->table_on;      
  GEN->table_size = PAR->table_size;  
  GEN->s[0] = PAR->s[0];              
  GEN->s[1] = PAR->s[1];
  GEN->table = NULL;
  GEN->f_table = NULL;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_ninv_info;
#endif
  return gen;
} 
int
_unur_ninv_check_par( struct unur_gen *gen )
{
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];
  GEN->CDFmin = GEN->Umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN->CDFmax = GEN->Umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;
  if (_unur_FP_greater(GEN->CDFmin, GEN->CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    return UNUR_ERR_GEN_DATA;
  }
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_ninv_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_ninv_gen*)clone->datap)
  struct unur_gen *clone;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NINV_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  if (GEN->table) {
    CLONE->table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->table, GEN->table, GEN->table_size * sizeof(double) );
    CLONE->f_table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->f_table, GEN->f_table, GEN->table_size * sizeof(double) );
  }
  return clone;
#undef CLONE
} 
void
_unur_ninv_free( struct unur_gen *gen )
{ 
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  if (GEN->table)   free(GEN->table);
  if (GEN->f_table) free(GEN->f_table);
  _unur_generic_free(gen);
} 
double
_unur_ninv_sample_regula( struct unur_gen *gen )
{
  return _unur_ninv_regula( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
} 
double 
_unur_ninv_sample_newton( struct unur_gen *gen )
{
  return _unur_ninv_newton( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
}
double
unur_ninv_eval_approxinvcdf( struct unur_gen *gen, double u )
{ 
  double x;
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);
  if ( u<0. || u>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }
  if (u<=0.) return DISTR.domain[0];
  if (u>=1.) return DISTR.domain[1];
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    x = _unur_ninv_newton(gen,u);
    break;
  case NINV_VARFLAG_REGULA:
  default:
    x = _unur_ninv_regula(gen,u);
    break;
  }
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];
  return x;
} 
double 
_unur_ninv_regula( struct unur_gen *gen, double u )
{ 
  double x1, x2, a, xtmp;
  double x2abs;          
  double f1, f2,fa, ftmp;
  double length;         
  double lengthabs;      
  double lengthsgn;      
  double step;           
  double dx;             
  int count = 0;         
  int i;                 
  int step_count;        
  double rel_u_resolution;  
  CHECK_NULL(gen, INFINITY);  COOKIE_CHECK(gen, CK_NINV_GEN, INFINITY);
  rel_u_resolution = (GEN->Umax - GEN->Umin) * GEN->rel_x_resolution;
  if (GEN->table_on) {
    if ( _unur_FP_same(GEN->CDFmin, GEN->CDFmax) ) {
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (u - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }
    if ( ! _unur_FP_is_minus_infinity(GEN->table[i]) ){
      x1 = GEN->table[i];
      f1 = GEN->f_table[i]; 
    }
    else{
      x1 = GEN->table[i+1] + (GEN->table[i+1] - GEN->table[i+2]);
      f1 = CDF(x1);
    }
    if( ! _unur_FP_is_infinity(GEN->table[i+1]) ){
      x2 = GEN->table[i+1];
      f2 = GEN->f_table[i+1];
    }
    else{
      x2 = GEN->table[i] + (GEN->table[i] - GEN->table[i-1]);
      f2 = CDF(x2);
    }
  }
  else { 
   x1 =  GEN->s[0];      
   f1 =  GEN->CDFs[0];
   x2 =  GEN->s[1];         
   f2 =  GEN->CDFs[1];
  }   
  if ( x1 >= x2 ) { 
    xtmp = x1; ftmp = f1;
    x1   = x2; f1   = f2;
    x2 = xtmp + DBL_EPSILON;
    f2 = CDF(x2); 
  }
  if ( x1 < DISTR.trunc[0] || x1 >= DISTR.trunc[1] ){
    x1 = DISTR.trunc[0];
    f1 = GEN->Umin;    
  }
  if ( x2 > DISTR.trunc[1] || x2 <= DISTR.trunc[0] ){
    x2 = DISTR.trunc[1];
    f2 = GEN->Umax;    
  }
  f1 -= u;
  f2 -= u;
  step = (GEN->s[1]-GEN->s[0]) * STEPFAC;
  step_count = 0;
  while ( f1*f2 > 0. ) {
    if ( f1 > 0. ) {     
      x2  = x1;  
      f2  = f1;
      x1 -= step;   
      f1 = CDF(x1) - u;
    }
    else {         
      x1  = x2;
      f1  = f2;
      x2 += step;
      f2 = CDF(x2) - u;
    }
    if (step_count < MAX_STEPS) {
      ++step_count;
      step *= 2.;
      if( step_count > 20 && step < 1.) step = 1.; 
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "Regula Falsi can't find interval with sign change");
      x2 = 0.5*x1 + 0.5*x2;
      x2 = _unur_max( x2, DISTR.trunc[0]);
      x2 = _unur_min( x2, DISTR.trunc[1]);
      return x2;
    }
  }   
  a  = x1;       
  fa = f1;
  for (i=0; TRUE ; i++) {
    count++;
    if ( f1*f2<0. && fabs(f1) < fabs(f2) ) {   
      xtmp = x1; ftmp = f1;
      x1 = x2;   f1 = f2;
      x2 = xtmp; f2 = ftmp;
    }
    x2abs = fabs(x2);   
    if ( f1*f2 < 0.) {  
      count = 0;   
      a  = x1;     
      fa = f1;
    }
    length = x2 - a;  
    lengthabs = fabs(length);
    lengthsgn = (length < 0.) ? -1. : 1.;
    if ( _unur_iszero(f2)                              
	 || _unur_FP_same(fa, f2)                      
	 || lengthabs <= GEN->rel_x_resolution * x2abs 
	 || lengthabs <= GEN->rel_x_resolution * GEN->rel_x_resolution
       	 || fabs(f2) <= rel_u_resolution ) {           
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug & NINV_DEBUG_SAMPLE)
	_unur_ninv_debug_sample_regula( gen,u,x2,f2,i );
#endif
      return x2;
    }
    if (i >= GEN->max_iter)
      break;
    dx = (_unur_FP_same(f1,f2)) ? length/2. : f2*(x2-x1)/(f2-f1) ;  
    if ( fabs(dx) < GEN->rel_x_resolution * x2abs ){
      dx = lengthsgn * 0.99 * GEN->rel_x_resolution * x2abs;
      while ( x2 == x2 - dx ){ 
	if ( dx != 2.*dx)    
	  dx = 2.*dx;
        else
	  dx = length/2.;    
      }
    }
    if ( count > 1 || i > I_CHANGE_TO_BISEC ||
        (lengthabs-GEN->rel_x_resolution*x2abs)/(dx*lengthsgn) <= 1. )
      dx = length/2.; 
    x1 = x2;       f1 = f2;
    x2 = x2-dx;    f2 = CDF(x2) - u; 
  }  
  if (i >= GEN->max_iter) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded");
    x2 = _unur_max( x2, DISTR.trunc[0]);
    x2 = _unur_min( x2, DISTR.trunc[1]);
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample_regula( gen,u,x2,f2,i );
#endif
  return x2;
} 
double
_unur_ninv_newton( struct unur_gen *gen, double U )
{ 
  double x;           
  double fx;          
  double dfx;         
  double fxabs;       
  double xtmp, fxtmp; 
  double xold;        
  double fxtmpabs;    
  double damp;        
  double step;        
  int i;              
  int flat_count;     
  const int MAX_FLAT_COUNT = 40;
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);
  if (GEN->table_on) {
    if ( _unur_FP_same(GEN->CDFmin,GEN->CDFmax) ) {
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (U - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }
    if (_unur_FP_is_infinity(GEN->table[i+1])) {
      x  = GEN->table[i];
      fx = GEN->f_table[i];
    }
    else {
      x  = GEN->table[i+1];
      fx = GEN->f_table[i+1];
    }
  }
  else { 
    x  = GEN->s[0];
    fx = GEN->CDFs[0];
  }
  if ( x < DISTR.trunc[0] ){
    x  = DISTR.trunc[0];
    fx = GEN->Umin;    
  }
  else if ( x > DISTR.trunc[1] ){
    x  = DISTR.trunc[1];
    fx = GEN->Umax;    
  }
  fx   -= U;
  dfx   = PDF(x);
  fxabs = fabs(fx);
  xold  = x;    
  damp = 2.;          
  step = 1.;
  for (i=0; i < GEN->max_iter; i++) {
    flat_count = 0;
    while (_unur_iszero(dfx)) {   
      if (_unur_iszero(fx))  
	break; 
      if (fx > 0.) {         
        xtmp = x - step; 
	xtmp = _unur_max( xtmp, DISTR.domain[0] );
      }
      else {
        xtmp  = x + step;
	xtmp = _unur_min( xtmp, DISTR.domain[1] );
      }
      fxtmp    = CDF(xtmp) - U;
      fxtmpabs = fabs(fxtmp);
      if ( fxtmpabs < fxabs ) {       
        step = 1.;     
        x     = xtmp;
        fx    = fxtmp;
      }
      else if ( fxtmp*fx < 0. ) {     
        step /= 2.;                      
      } 
      else {                          
        step *= 2.;    
        x     = xtmp;
        fx    = fxtmp;
      }  
      dfx   = PDF(x);
      fxabs = fabs(fx);
      if (flat_count < MAX_FLAT_COUNT)
	flat_count++;
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
		    "Newton's method can't leave flat region");
	x = _unur_max( x, DISTR.trunc[0]);
	x = _unur_min( x, DISTR.trunc[1]);
	return x;
      }
    }   
    step = 1.;   
    if (_unur_iszero(fx))  
      break;
    if (_unur_isfinite(dfx)) {
      do {    
	damp /= 2.;
	xtmp = x - damp * fx/dfx;
	xtmp = _unur_min( xtmp, DISTR.trunc[1] );
	xtmp = _unur_max( xtmp, DISTR.trunc[0] );
	fxtmp = CDF(xtmp) - U;
      } while (fabs(fxtmp) > fxabs * (1.+GEN->rel_x_resolution));   
    }
    else {
      xtmp = 0.5*(x + xold);
      fxtmp = CDF(xtmp) - U;
    }
    damp  = 2.;       
    xold  = x;        
    x     = xtmp;     
    fx    = fxtmp;    
    dfx   = PDF(x);   
    fxabs = fabs(fx); 
    if ( fabs(x-xold) <= fabs(x) * GEN->rel_x_resolution 
         && fabs(fx) < GEN->rel_x_resolution ) {
      break;   
    }
  }  
  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded");
  x = _unur_max( x, DISTR.trunc[0]);
  x = _unur_min( x, DISTR.trunc[1]);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample_newton(gen, U, x, fx, i);
#endif
  return x;
} 
int
_unur_ninv_create_table( struct unur_gen *gen )
{
  int i;
  double x;
  int table_size = GEN->table_size;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);
  GEN->table    = _unur_xrealloc( GEN->table,   table_size * sizeof(double));
  GEN->f_table  = _unur_xrealloc( GEN->f_table, table_size * sizeof(double));
  GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
  GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
  GEN->CDFs[0]  = CDF(GEN->s[0]);
  GEN->CDFs[1]  = CDF(GEN->s[1]);
  GEN->table_on = FALSE;
  GEN->table[0]              = DISTR.domain[0];
  GEN->f_table[0]            = GEN->CDFmin;    
  GEN->table[table_size-1]   = DISTR.domain[1];
  GEN->f_table[table_size-1] = GEN->CDFmax;    
  for (i=1; i<table_size/2; i++){
    x = GEN->CDFmin + i * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[i]   = _unur_ninv_regula(gen,x);
    GEN->f_table[i] = CDF(GEN->table[i]);
    x = GEN->CDFmin + (table_size-i-1) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size-1-i] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size-1-i] = CDF(GEN->table[table_size-1-i]);
    if (GEN->table[i] > -INFINITY) {
      GEN->s[0] = GEN->table[i];
      GEN->CDFs[0] = GEN->f_table[i];
    }
    if (GEN->table[table_size-1-i] < INFINITY) {
      GEN->s[1] = GEN->table[table_size-1-i];
      GEN->CDFs[1] = GEN->f_table[table_size-1-i];
    }
  }  
  if (table_size & 1) { 
    x = GEN->CDFmin + (table_size/2) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size/2] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size/2] = CDF(GEN->table[table_size/2]);
  }  
  GEN->table_on = TRUE;
  return UNUR_SUCCESS;
}  
int
_unur_ninv_compute_start( struct unur_gen *gen )
{
  double u;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);
  if( GEN->table_on )
    return UNUR_SUCCESS;
  if( !_unur_FP_same(GEN->s[0],GEN->s[1]) ) {
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    return UNUR_SUCCESS;
  }
  switch (gen->variant) {
  case NINV_VARFLAG_REGULA:
    GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);    
    u = GEN->CDFmin + 0.5*(1.-INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[0] = _unur_ninv_regula(gen,u);
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    u = GEN->CDFmin + 0.5*(1.+INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[1] = _unur_ninv_regula(gen,u);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    break;    
  case NINV_VARFLAG_NEWTON:
    GEN->s[0] = _unur_max( DISTR.domain[0], -9.987655 );
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]); 
    GEN->CDFs[1] = CDF(GEN->s[1]);
    u = 0.5 * (GEN->CDFmin + GEN->CDFmax);
    GEN->s[0] = _unur_ninv_regula(gen, u); 
    GEN->CDFs[0] = CDF(GEN->s[0]);
    break;    
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }  
  return UNUR_SUCCESS;
}  
#ifdef UNUR_ENABLE_LOGGING
void
_unur_ninv_debug_init( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ninv (numerical inversion of CDF)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(LOG,"%s: sampling routine = _unur_ninv_sample",gen->genid);
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    fprintf(LOG,"_newton\n");
    break;
  case NINV_VARFLAG_REGULA: default:
    fprintf(LOG,"_regula\n");
    break;
  }
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_ninv_debug_start(gen);
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_ninv_debug_start( const struct unur_gen *gen )
{
  FILE *LOG;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  if (GEN->table_on) {
    fprintf(LOG,"%s: use table (size = %d)\n",gen->genid,GEN->table_size);
    if (gen->debug & NINV_DEBUG_TABLE)
      for (i=0; i<GEN->table_size; i++)
	fprintf(LOG,"%s:\tx = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->table[i],GEN->f_table[i]);
  }
  else { 
    fprintf(LOG,"%s: starting points:\n",gen->genid);
    fprintf(LOG,"%s:\ts[0] = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->s[0],GEN->CDFs[0]);
    if (gen->variant & NINV_VARFLAG_REGULA)
      fprintf(LOG,"%s:\ts[1] = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->s[1],GEN->CDFs[1]);
  }
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_ninv_debug_sample_regula( const struct unur_gen *gen, double u, double x, double fx, int iter )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d iterations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN->max_iter);
} 
void
_unur_ninv_debug_sample_newton( const struct unur_gen *gen, double u, double x, double fx, int iter )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d iterations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN->max_iter);
} 
void 
_unur_ninv_debug_chg_truncated( const struct unur_gen *gen )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(LOG,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);
} 
#endif   
#ifdef UNUR_ENABLE_INFO
void
_unur_ninv_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  double n_iter;
  int samplesize = 10000;
  int use_newton = (gen->variant==NINV_VARFLAG_NEWTON) ? TRUE : FALSE;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = CDF");
  if (use_newton) 
    _unur_string_append(info," PDF");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n\n");
  _unur_string_append(info,"method: NINV (Numerical INVersion)\n");
  if (use_newton) 
    _unur_string_append(info,"   Newton method\n");
  else
    _unur_string_append(info,"   Regula falsi\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  n_iter = unur_test_count_pdf(gen,samplesize,FALSE,NULL)/(2.*samplesize);
  if (!use_newton) n_iter *= 2.;
  _unur_string_append(info,"   average number of iterations = %.2f  [approx.]\n", n_iter);
  if (GEN->table_on) {
    _unur_string_append(info,"   starting points = table of size %d\n", GEN->table_size);
  }
  else {
    _unur_string_append(info,"   starting points = ");
    if (use_newton) {
      _unur_string_append(info,"%g (CDF = %g)  %s\n", GEN->s[0], GEN->CDFs[0],
			  (gen->set & NINV_SET_START) ? "" : "[default]");
    }
    else {
      _unur_string_append(info,"%g, %g  (CDF = %g, %g)   %s\n",
			  GEN->s[0],GEN->s[1], GEN->CDFs[0],GEN->CDFs[1],
			  (gen->set & NINV_SET_START) ? "" : "[default]");
    }
  }
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    if (use_newton) 
      _unur_string_append(info,"   usenewton\n");
    else
      _unur_string_append(info,"   useregula  [default]\n");
    _unur_string_append(info,"   x_resolution = %g  %s\n", GEN->rel_x_resolution,
			(gen->set & NINV_SET_X_RESOLUTION) ? "" : "[default]");
    _unur_string_append(info,"   max_iter = %d  %s\n", GEN->max_iter,
			(gen->set & NINV_SET_MAX_ITER) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if (! (gen->set & NINV_SET_X_RESOLUTION) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase accuracy by decreasing \"x_resolution\".");
    if (! (gen->set & NINV_SET_MAX_ITER) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase \"max_iter\" if you encounter problems with accuracy.");
    _unur_string_append(info,"\n");
  }
} 
#endif   
