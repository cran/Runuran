/* Copyright (c) 2000-2017 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "x_gen.h"
#include "auto.h"
#include "auto_struct.h"
#include "cstd.h"
#include "dari.h"
#include "dgt.h"
#include "dstd.h"
#include "empk.h"
#include "hist.h"
#include "mvstd.h"
#include "tdr.h"
#include "vempk.h"
#define AUTO_SET_LOGSS          0x001u
#define GENTYPE "AUTO"         
static struct unur_gen *_unur_auto_init( struct unur_par *par );
static struct unur_gen *_unur_init_cont( struct unur_par *par );
static struct unur_gen *_unur_init_cvec( struct unur_par *par );
static struct unur_gen *_unur_init_discr( struct unur_par *par );
static struct unur_gen *_unur_init_cemp( struct unur_par *par );
static struct unur_gen *_unur_init_cvemp( struct unur_par *par );
#define PAR       ((struct unur_auto_par*)par->datap) 
struct unur_par *
unur_auto_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL(GENTYPE,distr,NULL);
  par = _unur_par_new( sizeof(struct unur_auto_par) );
  COOKIE_SET(par,CK_AUTO_PAR);
  par->distr    = distr;           
  par->method   = UNUR_METH_AUTO;  
  par->variant  = 0u;              
  par->set      = 0u;                  
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = par->urng;               
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_auto_init;
  return par;
} 
int 
unur_auto_set_logss( UNUR_PAR *par, int logss )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, AUTO );
  if (logss < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"log < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->logss = logss;
  par->set |= AUTO_SET_LOGSS;
  return UNUR_SUCCESS;
} 
struct unur_gen *
_unur_auto_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);
  if ( par->method != UNUR_METH_AUTO ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_AUTO_PAR,NULL);
  switch (par->distr->type) {
  case UNUR_DISTR_CONT:
    gen = _unur_init_cont( par );
    break;
  case UNUR_DISTR_CVEC:
    gen = _unur_init_cvec( par );
    break;
  case UNUR_DISTR_DISCR:
    gen = _unur_init_discr( par );
    break;
  case UNUR_DISTR_CEMP:
    gen = _unur_init_cemp( par );
    break;
  case UNUR_DISTR_CVEMP:
    gen = _unur_init_cvemp( par );
    break;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    gen = NULL;
    break;
  }
  if (gen) {
    gen->urng = par->urng;
    gen->urng_aux = par->urng_aux;
    gen->debug = par->debug;
  }
  _unur_par_free(par);
  return gen;
} 
struct unur_gen *
_unur_init_cont( struct unur_par *par_auto )
{
  struct unur_par *par;
  struct unur_gen *gen;
  do {
    par = unur_tdr_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;
    par = unur_cstd_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;
  } while (0);
  return gen;
} 
struct unur_gen *
_unur_init_cvec( struct unur_par *par_auto )
{
  struct unur_par *par;
  struct unur_gen *gen;
  par = unur_mvstd_new(par_auto->distr);
  gen = unur_init(par);
  return gen;
} 
struct unur_gen *
_unur_init_discr( struct unur_par *par_auto )
{
  struct unur_par *par;
  struct unur_gen *gen;
  do {
    if (par_auto->distr->data.discr.pv != NULL) {
      par = unur_dgt_new(par_auto->distr);
      gen = unur_init(par);
      if (gen) break;
    }
    if (par_auto->distr->data.discr.pmf != NULL) {
      par = unur_dari_new(par_auto->distr);
      gen = unur_init(par);
      if (gen) break;
      par = unur_dgt_new(par_auto->distr);
      gen = unur_init(par);
      if (gen) break;
    }
    par = unur_dstd_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;
  } while (0);
  return gen;
} 
struct unur_gen *
_unur_init_cemp( struct unur_par *par_auto )
{
  struct unur_par *par;
  struct unur_gen *gen;
  do {
    par = unur_empk_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;
    par = unur_hist_new(par_auto->distr);
    gen = unur_init(par);
    if (gen) break;
  } while(0);
  return gen;
} 
struct unur_gen *
_unur_init_cvemp( struct unur_par *par_auto )
{
  struct unur_par *par;
  struct unur_gen *gen;
  par = unur_vempk_new(par_auto->distr);
  gen = unur_init(par);
  return gen;
} 
#ifdef UNUR_ENABLE_LOGGING
#endif   
