/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h>
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#define MROU_HOOKE_RHO     (0.5)
#define MROU_HOOKE_EPSILON (1.e-7)
#define MROU_HOOKE_MAXITER (1000L)
#define MROU_RECT_SCALING (1.e-4)
static double _unur_mrou_rectangle_aux_vmax(double *x, void *p );
static double _unur_mrou_rectangle_aux_umin(double *x, void *p );
static double _unur_mrou_rectangle_aux_umax(double *x, void *p );
#define PDF(x)    _unur_cvec_PDF((x),(distr))    
double
_unur_mrou_rectangle_aux_vmax(double *x, void *p )
{
  struct MROU_RECTANGLE *rr;
  rr = p; 
  return -pow( _unur_cvec_PDF((x),(rr->distr)) ,
	       1./(1.+ rr->r * rr->dim) );
}
double
_unur_mrou_rectangle_aux_umin(double *x, void *p)
{
  struct MROU_RECTANGLE *rr;
  rr = p; 
  return ( (x[rr->aux_dim] - rr->center[rr->aux_dim])
	   * pow( _unur_cvec_PDF((x),(rr->distr)),
		  rr->r / (1.+ rr->r * rr->dim) ) );
}
double
_unur_mrou_rectangle_aux_umax(double *x, void *p)
{
  return (- _unur_mrou_rectangle_aux_umin(x,p)) ;
}
struct MROU_RECTANGLE *
_unur_mrou_rectangle_new( void )
{
  struct MROU_RECTANGLE *rr;
  rr = _unur_xmalloc(sizeof(struct MROU_RECTANGLE ));
  rr->distr  = NULL;
  rr->dim    = 0;
  rr->umin   = NULL;
  rr->umax   = NULL;
  rr->r      = 1;
  rr->bounding_rectangle = 1;
  rr->center = NULL;
  rr->genid  = "";
  return rr;
} 
int
_unur_mrou_rectangle_compute( struct MROU_RECTANGLE *rr )
{
  struct unur_funct_vgeneric faux; 
  double *xstart, *xend, *xumin, *xumax; 
  int d, dim;             
  int hooke_iters_vmax;   
  int hooke_iters_umin;   
  int hooke_iters_umax;   
  double scaled_epsilon;  
  int flag_finite = TRUE; 
  dim = rr->dim;
  xstart = _unur_xmalloc(dim * sizeof(double));
  xend   = _unur_xmalloc(dim * sizeof(double));
  xumin  = _unur_xmalloc(dim * sizeof(double));
  xumax  = _unur_xmalloc(dim * sizeof(double));
  if ( (rr->distr->set & UNUR_DISTR_SET_MODE) && (rr->distr->data.cvec.mode != NULL)) { 
    faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_vmax;
    faux.params = rr;
    rr->vmax = -faux.f(rr->distr->data.cvec.mode, faux.params);
  }
  else {
    faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_vmax;
    faux.params = rr;
    memcpy(xstart, rr->center, dim * sizeof(double));
    hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend,
				    MROU_HOOKE_RHO, MROU_HOOKE_EPSILON, MROU_HOOKE_MAXITER);
    rr->vmax = -faux.f(xend, faux.params);
    if (hooke_iters_vmax >= MROU_HOOKE_MAXITER) {
      scaled_epsilon = MROU_HOOKE_EPSILON * rr->vmax;
      if (scaled_epsilon>MROU_HOOKE_EPSILON) scaled_epsilon=MROU_HOOKE_EPSILON;
      memcpy(xstart, xend, dim * sizeof(double));
      hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend,
                                      MROU_HOOKE_RHO, scaled_epsilon , MROU_HOOKE_MAXITER);
      rr->vmax = -faux.f(xend, faux.params);
      if (hooke_iters_vmax >= MROU_HOOKE_MAXITER) {
        _unur_warning(rr->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (vmax)");
      }
    }
    rr->vmax = rr->vmax * ( 1+ MROU_RECT_SCALING);
  }
  flag_finite = _unur_isfinite(rr->vmax);
  if (rr->bounding_rectangle) {
    if (rr->umin == NULL || rr->umax == NULL) {
      free(xstart); free(xend); free(xumin); free(xumax);
      _unur_error(rr->genid,UNUR_ERR_NULL,"");
      return UNUR_ERR_NULL;
    }
    for (d=0; d<dim; d++) {
      rr->aux_dim  = d;
      memcpy(xstart, rr->center, dim * sizeof(double));
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umin;
      faux.params = rr;
      hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend,
				      MROU_HOOKE_RHO, MROU_HOOKE_EPSILON, MROU_HOOKE_MAXITER);
      rr->umin[d] = faux.f(xend, faux.params);
      memcpy(xumin, xend, dim * sizeof(double));
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umax;
      faux.params = rr;
      hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend,
				      MROU_HOOKE_RHO, MROU_HOOKE_EPSILON, MROU_HOOKE_MAXITER);
      rr->umax[d] = -faux.f(xend, faux.params);
      memcpy(xumax, xend, dim * sizeof(double));
      if (hooke_iters_umin >= MROU_HOOKE_MAXITER) {      
	scaled_epsilon = MROU_HOOKE_EPSILON * (rr->umax[d]-rr->umin[d]);
	if (scaled_epsilon>MROU_HOOKE_EPSILON) scaled_epsilon=MROU_HOOKE_EPSILON;
	faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umin;
	faux.params = rr;
	memcpy(xstart, xumin, dim * sizeof(double));
	hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend,
					MROU_HOOKE_RHO, scaled_epsilon , MROU_HOOKE_MAXITER);
	rr->umin[d] = faux.f(xend, faux.params);
	if (hooke_iters_umin >= MROU_HOOKE_MAXITER) {
	  _unur_warning(rr->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umin)");
	}
      }
      if (hooke_iters_umax >= MROU_HOOKE_MAXITER) {
	scaled_epsilon = MROU_HOOKE_EPSILON * (rr->umax[d]-rr->umin[d]);
	if (scaled_epsilon>MROU_HOOKE_EPSILON) scaled_epsilon=MROU_HOOKE_EPSILON;
	faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umax;
	faux.params = rr;
	memcpy(xstart, xumax, dim * sizeof(double));
	hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend,
					MROU_HOOKE_RHO, scaled_epsilon , MROU_HOOKE_MAXITER);
	rr->umin[d] = faux.f(xend, faux.params);
	if (hooke_iters_umax >= MROU_HOOKE_MAXITER) {
	  _unur_warning(rr->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umax)");
	}
      }
      rr->umin[d] = rr->umin[d] - (rr->umax[d]-rr->umin[d])*MROU_RECT_SCALING/2.;
      rr->umax[d] = rr->umax[d] + (rr->umax[d]-rr->umin[d])*MROU_RECT_SCALING/2.;
      flag_finite = flag_finite && _unur_isfinite(rr->umin[d]) && _unur_isfinite(rr->umax[d]);
    }
  }
  free(xstart); free(xend); free(xumin); free(xumax);
  if (rr->vmax <= 0.) {
    _unur_error("RoU",UNUR_ERR_DISTR_DATA,"cannot find bounding rectangle");
    return UNUR_ERR_DISTR_DATA;
  }
  return (flag_finite ? UNUR_SUCCESS : UNUR_ERR_INF);
} 
#undef PDF
#undef MROU_HOOKE_RHO
#undef MROU_HOOKE_EPSILON
#undef MROU_HOOKE_MAXITER
#undef MROU_RECT_SCALING
