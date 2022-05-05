/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include <parser/functparser_source.h>
#include "distr_source.h"
#include "distr.h"
#include "discr.h"
#define DISTR distr->data.discr
#define MAX_PMF_DOMAIN_FOR_UPD_PMFSUM   (1000)
static double _unur_distr_discr_eval_pmf_tree( int k, const struct unur_distr *distr );
static double _unur_distr_discr_eval_cdf_tree( int k, const struct unur_distr *distr );
static void _unur_distr_discr_free( struct unur_distr *distr );
static int _unur_distr_discr_find_mode( struct unur_distr *distr );
struct unur_distr *
unur_distr_discr_new( void )
{
  register struct unur_distr *distr;
  register int i;
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;
  COOKIE_SET(distr,CK_DISTR_DISCR);
  distr->type = UNUR_DISTR_DISCR;
  distr->id = UNUR_DISTR_GENERIC;
  distr->dim = 1;   
  distr->destroy = _unur_distr_discr_free;
  distr->clone = _unur_distr_discr_clone;
  DISTR.pv        = NULL;          
  DISTR.n_pv      = 0;             
  DISTR.pmf       = NULL;          
  DISTR.cdf       = NULL;          
  DISTR.invcdf    = NULL;          
  DISTR.init      = NULL;          
  DISTR.set_params= NULL;          
  DISTR.n_params  = 0;             
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;
  DISTR.norm_constant = 1.;        
  DISTR.trunc[0] = DISTR.domain[0] = 0;         
  DISTR.trunc[1] = DISTR.domain[1] = INT_MAX;   
  DISTR.mode     = 0;              
  DISTR.upd_mode = _unur_distr_discr_find_mode;  
  DISTR.sum     = 1.;              
  DISTR.upd_sum = NULL;            
  DISTR.pmftree    = NULL;         
  DISTR.cdftree    = NULL;         
  return distr;
} 
struct unur_distr *
_unur_distr_discr_clone( const struct unur_distr *distr )
{
#define CLONE clone->data.discr
  struct unur_distr *clone;
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  memcpy( clone, distr, sizeof( struct unur_distr ) );
  CLONE.pmftree  = (DISTR.pmftree) ? _unur_fstr_dup_tree(DISTR.pmftree) : NULL;
  CLONE.cdftree  = (DISTR.cdftree) ? _unur_fstr_dup_tree(DISTR.cdftree) : NULL;
  if (DISTR.pv) {
    CLONE.pv = _unur_xmalloc( DISTR.n_pv * sizeof(double) );
    memcpy( CLONE.pv, DISTR.pv, DISTR.n_pv * sizeof(double) );
  }
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }
  return clone;
#undef CLONE
} 
void
_unur_distr_discr_free( struct unur_distr *distr )
{
  if( distr == NULL ) 
    return;
  _unur_check_distr_object( distr, DISCR, RETURN_VOID );
  if (DISTR.pmftree)  _unur_fstr_free(DISTR.pmftree);
  if (DISTR.cdftree)  _unur_fstr_free(DISTR.cdftree);
  if (DISTR.pv) free( DISTR.pv );
  if (distr->name_str) free(distr->name_str);
  COOKIE_CLEAR(distr);
  free( distr );
} 
int
unur_distr_discr_set_pv( struct unur_distr *distr, const double *pv, int n_pv )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pmf != NULL || DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PMF/CDF given, cannot set PV");
    return UNUR_ERR_DISTR_SET;
  }
  if (n_pv < 0) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"length of PV");
    return UNUR_ERR_DISTR_SET;
  }
  if ( (DISTR.domain[0] > 0) && ((unsigned)DISTR.domain[0] + (unsigned)n_pv > INT_MAX) ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"length of PV too large, overflow");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.domain[1] = DISTR.domain[0] + n_pv - 1;
  DISTR.pv = _unur_xrealloc( DISTR.pv, n_pv * sizeof(double) );
  if (!DISTR.pv) return UNUR_ERR_MALLOC;
  memcpy( DISTR.pv, pv, n_pv * sizeof(double) );
  DISTR.n_pv = n_pv;
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_make_pv( struct unur_distr *distr )
{
  double *pv;          
  int n_pv;            
  double cdf, cdf_old; 
  double thresh_cdf;   
  int valid;           
  int i;
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );
  if ( DISTR.pmf == NULL && DISTR.cdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"PMF or CDF");
    return 0;
  }
  if (DISTR.pv != NULL) {
    free(DISTR.pv); DISTR.n_pv = 0;
  }
  if ((unsigned)DISTR.domain[1] - (unsigned)DISTR.domain[0] < UNUR_MAX_AUTO_PV ) {
    n_pv = DISTR.domain[1] - DISTR.domain[0] + 1;
    pv = _unur_xmalloc( n_pv * sizeof(double) );
    if (DISTR.pmf) {
      for (i=0; i<n_pv; i++)
	pv[i] = _unur_discr_PMF(DISTR.domain[0]+i,distr);
    }
    else if (DISTR.cdf) {
      cdf_old = 0.;
      for (i=0; i<n_pv; i++) {
	cdf = _unur_discr_CDF(DISTR.domain[0]+i,distr);
	pv[i] = cdf - cdf_old;
	cdf_old = cdf;
      }
    }
    valid = TRUE;
  }
  else {
#define MALLOC_SIZE 1000 
    int n_alloc;         
    int max_alloc;       
    int size_alloc;      
    if ( (DISTR.domain[0] <= 0) || (INT_MAX - DISTR.domain[0] >= UNUR_MAX_AUTO_PV - 1) ) {
      size_alloc = MALLOC_SIZE;
      max_alloc = UNUR_MAX_AUTO_PV;
    }
    else { 
      size_alloc = max_alloc = INT_MAX - DISTR.domain[0];
    }
    n_pv = 0;
    pv = NULL;
    valid = FALSE;  
    cdf = 0.;       
    cdf_old = 0.;   
    thresh_cdf = (distr->set & UNUR_DISTR_SET_PMFSUM) ? (1.-1.e-8)*DISTR.sum : UNUR_INFINITY;
    for (n_alloc = size_alloc; n_alloc <= max_alloc; n_alloc += size_alloc) {
      pv = _unur_xrealloc( pv, n_alloc * sizeof(double) );
      if (DISTR.pmf) {
	for (i=0; i<size_alloc; i++) {
	  cdf += pv[n_pv] = _unur_discr_PMF(DISTR.domain[0]+n_pv,distr);
	  n_pv++;
	  if (cdf > thresh_cdf) { valid = TRUE; break; }
	}
      }
      else if (DISTR.cdf) {
	for (i=0; i<size_alloc; i++) {
	  cdf = _unur_discr_CDF(DISTR.domain[0]+n_pv,distr);
	  pv[n_pv] = cdf - cdf_old;
	  cdf_old = cdf;
	  n_pv++;
	  if (cdf > thresh_cdf) { valid = TRUE; break; }
	}
      }	  
      if (cdf > thresh_cdf) break;
    }
    if (distr->set & UNUR_DISTR_SET_PMFSUM) {
      if (valid != TRUE)
	_unur_warning(distr->name,UNUR_ERR_DISTR_GET,"PV truncated");
    }
    else { 
      valid = TRUE;
      DISTR.sum = cdf;
      distr->set |= UNUR_DISTR_SET_PMFSUM;
    }
#undef MALLOC_SIZE
  }
  DISTR.pv = pv;
  DISTR.n_pv = n_pv;
  DISTR.domain[1] = DISTR.domain[0] + n_pv - 1;
  return (valid) ? n_pv : -n_pv;
} 
int 
unur_distr_discr_get_pv( const struct unur_distr *distr, const double **pv )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );
  *pv = (DISTR.pv) ? DISTR.pv : NULL;
  return DISTR.n_pv;
} 
double
unur_distr_discr_eval_pv( int k, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_INFINITY );
  _unur_check_distr_object( distr, DISCR, UNUR_INFINITY );
  if (DISTR.pv != NULL) {
    if (k < DISTR.domain[0] || k > DISTR.domain[1])
      return 0.;
    else
      return (DISTR.pv[k-DISTR.domain[0]]);
  }
  if (DISTR.pmf != NULL) {
    double px = _unur_discr_PMF(k,distr);
    if (_unur_isnan(px)) {
      _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"PMF returns NaN");
      return 0.;
    }
    else
      return px;
  }
  _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
  return UNUR_INFINITY;
} 
int
unur_distr_discr_set_pmf( struct unur_distr *distr, UNUR_FUNCT_DISCR *pmf )
{
  _unur_check_NULL( NULL, distr,  UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pmf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pv != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"delete exisiting PV");
    free(DISTR.pv); DISTR.n_pv = 0;
  }
  if (DISTR.pmf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PMF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.pmf = pmf;
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_set_cdf( struct unur_distr *distr, UNUR_FUNCT_DISCR *cdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, cdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pv != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"delete exisiting PV");
    free(DISTR.pv); DISTR.n_pv = 0;
  }
  if (DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.cdf = cdf;
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_set_invcdf( struct unur_distr *distr, UNUR_IFUNCT_DISCR *invcdf )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, invcdf,UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (DISTR.invcdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of inverse CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_INVALID;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  DISTR.invcdf = invcdf;
  return UNUR_SUCCESS;
} 
UNUR_FUNCT_DISCR *
unur_distr_discr_get_pmf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  return DISTR.pmf;
} 
UNUR_FUNCT_DISCR *
unur_distr_discr_get_cdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  return DISTR.cdf;
} 
UNUR_IFUNCT_DISCR *
unur_distr_discr_get_invcdf( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  return DISTR.invcdf;
} 
double
unur_distr_discr_eval_pmf( int k, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_INFINITY );
  _unur_check_distr_object( distr, DISCR, UNUR_INFINITY );
  if (DISTR.pmf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_INFINITY;
  }
  return _unur_discr_PMF(k,distr);
} 
double
unur_distr_discr_eval_cdf( int k, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_INFINITY );
  _unur_check_distr_object( distr, DISCR, UNUR_INFINITY );
  if (DISTR.cdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_INFINITY;
  }
  return _unur_discr_CDF(k,distr);
} 
int
unur_distr_discr_eval_invcdf( double u, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INT_MAX );
  _unur_check_distr_object( distr, DISCR, INT_MAX );
  if (DISTR.invcdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INT_MAX;
  }
  if (u<=0.)
    return DISTR.domain[0];
  if (u>=1.)
    return DISTR.domain[1];
  else
    return _unur_discr_invCDF(u,distr);
} 
int
unur_distr_discr_set_pmfstr( struct unur_distr *distr, const char *pmfstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, pmfstr, UNUR_ERR_NULL );
  if (DISTR.pv != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PV given, cannot set PMF");
    return UNUR_ERR_DISTR_SET;
  }
  if (DISTR.pmf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PMF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_DATA;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.pmftree = _unur_fstr2tree(pmfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.pmf  = _unur_distr_discr_eval_pmf_tree;
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_set_cdfstr( struct unur_distr *distr, const char *cdfstr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, cdfstr, UNUR_ERR_NULL );
  if (DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }
  if (distr->base) return UNUR_ERR_DISTR_DATA;
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if ( (DISTR.cdftree = _unur_fstr2tree(cdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.cdf  = _unur_distr_discr_eval_cdf_tree;
  return UNUR_SUCCESS;
} 
double
_unur_distr_discr_eval_pmf_tree( int k, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_INFINITY );
  _unur_check_distr_object( distr, DISCR, UNUR_INFINITY );
  return ((DISTR.pmftree) ? _unur_fstr_eval_tree(DISTR.pmftree,(double)k) : 0.);
} 
double
_unur_distr_discr_eval_cdf_tree( int k, const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_INFINITY );
  _unur_check_distr_object( distr, DISCR, UNUR_INFINITY );
  return ((DISTR.cdftree) ? _unur_fstr_eval_tree(DISTR.cdftree,(double)k) : 0.);
} 
char *
unur_distr_discr_get_pmfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  _unur_check_NULL( NULL, DISTR.pmftree, NULL );
  return _unur_fstr_tree2string(DISTR.pmftree,"x","PMF",TRUE);
} 
char *
unur_distr_discr_get_cdfstr( const struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  _unur_check_NULL( NULL, DISTR.cdftree, NULL );
  return _unur_fstr_tree2string(DISTR.cdftree,"x","CDF",TRUE);
} 
int
unur_distr_discr_set_pmfparams( struct unur_distr *distr, const double *params, int n_params )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name, params, UNUR_ERR_NULL);
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  if (DISTR.set_params)
    return (DISTR.set_params(distr,params,n_params));
  DISTR.n_params = n_params;
  if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_get_pmfparams( const struct unur_distr *distr, const double **params )
{
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );
  *params = (DISTR.n_params) ? DISTR.params : NULL;
  return DISTR.n_params;
} 
int
unur_distr_discr_set_domain( struct unur_distr *distr, int left, int right )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (left >= right) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.trunc[0] = DISTR.domain[0] = left;
  DISTR.trunc[1] = DISTR.domain[1] = (DISTR.pv == NULL) ? right : left+DISTR.n_pv-1;
  distr->set |= UNUR_DISTR_SET_DOMAIN;
  distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_TRUNCATED | 
		  UNUR_DISTR_SET_MASK_DERIVED );
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_get_domain( const struct unur_distr *distr, int *left, int *right )
{
  *left = INT_MIN;
  *right = INT_MAX;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];
  return UNUR_SUCCESS;
} 
int
unur_distr_discr_set_mode( struct unur_distr *distr, int mode )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  DISTR.mode = mode;
  distr->set |= UNUR_DISTR_SET_MODE;
  return UNUR_SUCCESS;
} 
int 
unur_distr_discr_upd_mode( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (DISTR.upd_mode == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
  if ((DISTR.upd_mode)(distr)==UNUR_SUCCESS) {
    distr->set |= UNUR_DISTR_SET_MODE;
    return UNUR_SUCCESS;
  }
  else {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }
} 
int
unur_distr_discr_get_mode( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, INT_MAX );
  _unur_check_distr_object( distr, DISCR, INT_MAX );
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    if (DISTR.upd_mode == NULL) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INT_MAX;
    }
    else {
      if (unur_distr_discr_upd_mode(distr)!=UNUR_SUCCESS) {
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return INT_MAX;
      }
    }
  }
  return DISTR.mode;
} 
int
unur_distr_discr_set_pmfsum( struct unur_distr *distr, double sum )
{
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (sum <= 0.) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"pmf sum <= 0");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.sum = sum;
  distr->set |= UNUR_DISTR_SET_PMFSUM;
  return UNUR_SUCCESS;
} 
int 
unur_distr_discr_upd_pmfsum( struct unur_distr *distr )
{
  double sum = 0.;
  int k, left, right, length;
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_SET );
  distr->set |= UNUR_DISTR_SET_PMFSUM;
  if (DISTR.upd_sum != NULL) {
    if ((DISTR.upd_sum)(distr)==UNUR_SUCCESS)
      return UNUR_SUCCESS;
  }
  left  = DISTR.domain[0];
  right = DISTR.domain[1];
  length = right - left;
  if (DISTR.cdf != NULL) {
    if (left > INT_MIN) left -= 1;
    DISTR.sum = _unur_discr_CDF(right,distr) - _unur_discr_CDF(left,distr);
    return UNUR_SUCCESS;
  }
  if (DISTR.pv != NULL) {
    for (k = 0; k<= length; k++)
      sum += DISTR.pv[k];
    DISTR.sum = sum;
    return UNUR_SUCCESS;
  }
  if (DISTR.pmf != NULL && length > 0 && length <= MAX_PMF_DOMAIN_FOR_UPD_PMFSUM) {
    for (k = left; k<= right; k++)
      sum += _unur_discr_PMF(k,distr);
    DISTR.sum = sum;
    return UNUR_SUCCESS;
  }
  distr->set &= ~UNUR_DISTR_SET_PMFSUM;
  _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"Cannot compute sum");
  return UNUR_ERR_DISTR_DATA;
} 
double
unur_distr_discr_get_pmfsum( struct unur_distr *distr )
{
  _unur_check_NULL( NULL, distr, UNUR_INFINITY );
  _unur_check_distr_object( distr, DISCR, UNUR_INFINITY );
  if ( !(distr->set & UNUR_DISTR_SET_PMFSUM) ) {
    if ( unur_distr_discr_upd_pmfsum(distr) != UNUR_SUCCESS ) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"sum");
      return UNUR_INFINITY;
    }
  }
  return DISTR.sum;
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_discr_debug( const struct unur_distr *distr, const char *genid, unsigned printvector )
{
  FILE *LOG;
  int i;
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_DISCR,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = discrete univariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);
  if ( DISTR.pmf ) {
    fprintf(LOG,"%s:\tPMF with %d argument(s)\n",genid,DISTR.n_params);
    for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);
  }
  if (DISTR.n_pv>0) {
    fprintf(LOG,"%s:\tprobability vector of length %d",genid,DISTR.n_pv);
    if (printvector) {
      for (i=0; i<DISTR.n_pv; i++) {
	if (i%10 == 0)
	  fprintf(LOG,"\n%s:\t",genid);
	fprintf(LOG,"  %.5f",DISTR.pv[i]);
      }
    }
    fprintf(LOG,"\n%s:\n",genid);
  }
  if ( DISTR.pmf ) {
    fprintf(LOG,"%s:\tdomain for pmf = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[1]);
    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
    fprintf(LOG,"\n%s:\n",genid);
  }
  if (DISTR.n_pv>0) {
    fprintf(LOG,"%s:\tdomain for pv = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[0]-1+DISTR.n_pv);
    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
    fprintf(LOG,"\n%s:\n",genid);
  }
  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %d\n",genid,DISTR.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);
  if (distr->set & UNUR_DISTR_SET_PMFSUM)
    fprintf(LOG,"%s:\tsum over PMF = %g\n",genid,DISTR.sum);
  else
    fprintf(LOG,"%s:\tsum over PMF unknown\n",genid);
  fprintf(LOG,"%s:\n",genid);
} 
#endif    
int 
_unur_distr_discr_find_mode(struct unur_distr *distr )
{
#define N_TRIALS  (100)
  int x[3];                       
  double fx[3];                   
  int xnew;                       
  double fxnew;                   
  int step;                       
  int this, other;                
  int cutthis;                    
  const double r = (sqrt(5.)-1.)/2.;     
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  x[0] = DISTR.domain[0];
  x[2] = DISTR.domain[1];
  fx[0] = unur_distr_discr_eval_pv(x[0], distr);
  fx[2] = unur_distr_discr_eval_pv(x[2], distr);
  if (x[2] <= x[0] + 1) {
    DISTR.mode = (fx[0] <= fx[2]) ? x[2] : x[0];
    distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
    return UNUR_SUCCESS;
  }
  x[1]  = (x[0]/2) + (x[2]/2);
  if (x[1]<=x[0]) x[1]++;
  if (x[1]>=x[2]) x[1]--;
  fx[1] = unur_distr_discr_eval_pv(x[1], distr); 
  if ( !(fx[1]>0.)) {
    xnew = (DISTR.domain[0]!=INT_MIN) ? DISTR.domain[0] : 0;
    for (step = 1; step < N_TRIALS; step++) {
      xnew += step;
      if (xnew >= DISTR.domain[1]) break;
      if ((fxnew = unur_distr_discr_eval_pv(xnew,distr)) > 0.) {
	x[1] = xnew; fx[1] = fxnew; break;
      }
    }
  }
  if ( !(fx[1]>0.) && DISTR.domain[0]!=0) {
    xnew = 0;
    for (step = 1; step < N_TRIALS; step++) {
      xnew += step;
      if (xnew >= DISTR.domain[1]) break;
      if ((fxnew = unur_distr_discr_eval_pv(xnew,distr)) > 0.) {
	x[1] = xnew; fx[1] = fxnew; break;
      }
    }
  }
  if ( !(fx[1]>0.) && DISTR.domain[1]!=INT_MAX) {
    xnew = DISTR.domain[1];
    for (step = 1; step < N_TRIALS; step++) {
      xnew -= step;
      if (xnew <= DISTR.domain[0]) break;
      if ((fxnew = unur_distr_discr_eval_pv(xnew,distr)) > 0.) {
	x[1] = xnew; fx[1] = fxnew; break;
      }
    }
  }
  if ( !(fx[1]>0.)) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,
		"find_mode(): no positive entry in PV found");
    return UNUR_ERR_DISTR_DATA;
  }
  if (fx[1]<fx[0] && fx[1]<fx[2]) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA, "find_mode(): PV not unimodal");
    return UNUR_ERR_DISTR_DATA;
  }
  while (1) {
    if (x[0]+1 >= x[1] && x[1] >= x[2]-1) {
      DISTR.mode = (fx[0]>fx[2]) ? x[0] : x[2];
      if (fx[1]>DISTR.mode) DISTR.mode = x[1];
      distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
      return UNUR_SUCCESS;
    } 
    xnew  = (int) (r*x[0] + (1.-r)*x[2]);
    if (xnew == x[0])  ++xnew;
    if (xnew == x[2])  --xnew;
    if (xnew == x[1])  xnew += (x[1]-1==x[0]) ? 1 : -1;
    if (xnew < x[1]) {
      this = 0; other = 2; } 
    else {
      this = 2; other = 0; } 
    fxnew = unur_distr_discr_eval_pv(xnew,distr);
    if ( fxnew < fx[0] && fxnew < fx[2] ) {
      _unur_error(distr->name,UNUR_ERR_DISTR_DATA, "find_mode(): PV not unimodal");
      return UNUR_ERR_DISTR_DATA;
    }
    do {
      if (!_unur_FP_same(fxnew,fx[1])) {
	cutthis = (fxnew > fx[1]) ? FALSE : TRUE;
	break;
      }
      if (fx[this]  > fx[1]) { cutthis = FALSE; break; }
      if (fx[other] > fx[1]) { cutthis = TRUE;  break; }
      for (step = 1; step < N_TRIALS && xnew >= x[0]  && xnew <= x[2]; step++) {
	xnew += (this==0) ? -1 : 1;
	fxnew = unur_distr_discr_eval_pv(xnew,distr);
	if (_unur_FP_less(fxnew,fx[1])) {
	  DISTR.mode = x[1];
	  distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
	  return UNUR_SUCCESS;
	}
      }
      _unur_error(distr->name,UNUR_ERR_DISTR_DATA, "find_mode(): PV not unimodal");
      return UNUR_ERR_DISTR_DATA;
    } while (0);
    if (cutthis) {
      x[this] = xnew; fx[this] = fxnew;
    }
    else {
      x[other] = x[1]; fx[other] = fx[1];
      x[1] = xnew; fx[1] = fxnew;
    }
  }   
} 
#undef DISTR
