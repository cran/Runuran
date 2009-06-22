/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "lobatto_source.h"
#include "lobatto_struct.h"
static double 
_unur_lobatto5_simple (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		       double x, double h, double *fx);
static double
_unur_lobatto5_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			 double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			 struct unur_lobatto_table *Itable);
static double 
_unur_lobatto5_recursion (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			  double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			  double int1, double fl, double fr, double fc,
			  int *W_accuracy, struct unur_lobatto_table *Itable);
static int
_unur_lobatto_table_append (struct unur_lobatto_table *Itable, double x, double u);
static void
_unur_lobatto_table_resize (struct unur_lobatto_table *Itable);
#define FKT(x)  (funct((x),gen))      
#define W1 (0.17267316464601146)   
#define W2 (1.-W1)
double
_unur_lobatto5_simple (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		       double x, double h, double *fx)
{ 
  double fl, fr;
  if (fx==NULL) {
    fl = FKT(x);
    fr = FKT(x+h);
  }
  else {
    fl = (*fx>=0.) ? *fx : FKT(x);
    fr = *fx = FKT(x+h);
  }
  return (9*(fl+fr)+49.*(FKT(x+h*W1)+FKT(x+h*W2))+64*FKT(x+h/2.))*h/180.;
} 
double
_unur_lobatto_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			double x, double h, double tol, UNUR_LOBATTO_ERROR uerror)
{
  return _unur_lobatto5_adaptive(funct,gen,x,h,tol,uerror,NULL); 
} 
double
_unur_lobatto5_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen, 
			 double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			 struct unur_lobatto_table *Itable)
{
  double fl, fc, fr;  
  double int1, int2;  
  int W_accuracy = FALSE; 
  if (_unur_iszero(h))
    return 0.;
  if (!_unur_isfinite(x+h)) {
    _unur_error(gen->genid,UNUR_ERR_INF,"boundaries of integration domain not finite");
    return INFINITY;
  }
  fl = FKT(x);
  fc = FKT(x+h/2.);
  fr = FKT(x+h);
  int1 = (9*(fl+fr)+49.*(FKT(x+h*W1)+FKT(x+h*W2))+64*fc)*h/180.;
  int2 = _unur_lobatto5_recursion(funct,gen,x,h,tol,uerror,int1,fl,fc,fr,&W_accuracy,Itable);
  if (W_accuracy)
    _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
		  "numeric integration did not reach full accuracy");
  return int2;
} 
double
_unur_lobatto5_recursion (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			  double x, double h, double tol, UNUR_LOBATTO_ERROR uerror,
			  double int1, double fl, double fc, double fr,
			  int *W_accuracy, struct unur_lobatto_table *Itable)
{
  double flc, frc;    
  double int2;        
  double intl, intr;  
  double ierror;      
  flc = FKT(x+h/4);
  frc = FKT(x+3*h/4);
  intl = (9*(fl+fc)+49.*(FKT(x+h*W1*0.5)+FKT(x+h*W2*0.5))+64*flc)*h/360.;
  intr = (9*(fc+fr)+49.*(FKT(x+h*(0.5+W1*0.5))+FKT(x+h*(0.5+W2*0.5)))+64*frc)*h/360.;
  int2 = intl + intr;
  if (uerror!=NULL)
    ierror = uerror(gen, fabs(int1-int2), x+h/2.);
  else 
    ierror = fabs(int1-int2);
  if (ierror >= tol) {
    if (_unur_FP_equal(x+h/2.,x)) {
      *W_accuracy = TRUE;
    }
    else {
      return ( _unur_lobatto5_recursion(funct,gen,x,h/2,tol/1.,uerror,
					intl,fl,flc,fc, W_accuracy,Itable) +
	       _unur_lobatto5_recursion(funct,gen,x+h/2,h/2,tol/1.,uerror,
					intr,fc,frc,fr, W_accuracy,Itable) );
    }
  }
  if (Itable) {
    _unur_lobatto_table_append(Itable, x+h/2., intl);
    _unur_lobatto_table_append(Itable, x+h, intr);
  }
  return int2;
} 
double
_unur_lobatto_eval_diff (struct unur_lobatto_table *Itable, double x, double h, double *fx)
{
  int cur;                    
  double x1;                  
  double Q;                   
  struct unur_lobatto_nodes *values;
  int n_values;
#define clear_fx() if(fx!=NULL){*fx=-1.;}
  CHECK_NULL(Itable,INFINITY);
  values = Itable->values;
  n_values = Itable->n_values;
  if (!_unur_isfinite(x+h)) {
    clear_fx();
    return INFINITY;
  }
  if (x < Itable->bleft || x+h > Itable->bright) {
    clear_fx();
    return _unur_lobatto5_adaptive(Itable->funct, Itable->gen, x, h, 
				   Itable->tol, Itable->uerror, NULL);
  }
  cur = Itable->cur_iv;
  while (cur < n_values && values[cur].x < x)
    ++cur;
  if (cur >= n_values) {
    clear_fx();
    return _unur_lobatto5_adaptive(Itable->funct, Itable->gen, x, h, 
				   Itable->tol, Itable->uerror, NULL);
  }
  x1 = values[cur].x;
  ++cur;
  if (cur >= n_values ||
      values[cur].x > x+h) {
    return _unur_lobatto5_simple(Itable->funct, Itable->gen, x, h, fx);
  }
  Q = _unur_lobatto5_simple(Itable->funct, Itable->gen, x, x1-x, fx);
  do {
    Q += values[cur].u;
    x1 = values[cur].x;
    ++cur;
  } while (cur < n_values && values[cur].x <= x+h);
  clear_fx();
  if (cur >= n_values) {
    Q += _unur_lobatto5_adaptive(Itable->funct, Itable->gen, x1, x+h-x1,
				 Itable->tol, Itable->uerror, NULL);
  }
  else {
    Q += _unur_lobatto5_simple(Itable->funct, Itable->gen, x1, x+h-x1, fx);
  }
  return Q;
#undef clear_fx
} 
double
_unur_lobatto_integral (struct unur_lobatto_table *Itable)
{
  CHECK_NULL(Itable,INFINITY);
  return Itable->integral;
} 
struct unur_lobatto_table *
_unur_lobatto_init (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		    double left, double center, double right,
		    double tol, UNUR_LOBATTO_ERROR uerror, int size)
{
  struct unur_lobatto_table *Itable;
  if (size<2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"size<2");
    return NULL;
  }
  Itable = _unur_xmalloc( sizeof(struct unur_lobatto_table) );
  Itable->values = _unur_xmalloc(size * sizeof(struct unur_lobatto_nodes) );
  Itable->size = size;
  Itable->n_values = 0;
  Itable->cur_iv = 0;
  Itable->funct = funct;
  Itable->gen = gen;
  Itable->bleft = left;
  Itable->bright = right;
  Itable->tol = tol;
  Itable->uerror = uerror;
  _unur_lobatto_table_append(Itable,left,0.);
  Itable->integral = 
    _unur_lobatto5_adaptive(funct, gen, left, center-left, tol, uerror, Itable ) +
    _unur_lobatto5_adaptive(funct, gen, center, right-center, tol, uerror, Itable );
  _unur_lobatto_table_resize(Itable);
  return Itable;
} 
int _unur_lobatto_find_linear (struct unur_lobatto_table *Itable, double x)
{
  if (Itable != NULL) {
    while (Itable->cur_iv < Itable->n_values &&
	   Itable->values[Itable->cur_iv].x < x) 
      ++(Itable->cur_iv);
    return UNUR_SUCCESS;
  }
  else {
    return UNUR_ERR_SILENT;
  }
} 
int
_unur_lobatto_table_append (struct unur_lobatto_table *Itable, double x, double u)
{
  if (Itable==NULL)
    return UNUR_ERR_NULL;
  if (Itable->n_values >= Itable->size - 1)
    return UNUR_ERR_GENERIC;
  Itable->values[Itable->n_values].x = x;
  Itable->values[Itable->n_values].u = u;
  ++(Itable->n_values);
  return UNUR_SUCCESS;
} 
void
_unur_lobatto_table_resize (struct unur_lobatto_table *Itable)
{
  if (Itable) {
    Itable->size = Itable->n_values;
    Itable->values = _unur_xrealloc(Itable->values,
				    Itable->size * sizeof(struct unur_lobatto_nodes));
  }
} 
void
_unur_lobatto_free (struct unur_lobatto_table **Itable)
{
  if (*Itable) {
    free ((*Itable)->values);
    free (*Itable);
    *Itable = NULL;
  }
} 
void
_unur_lobatto_debug_table (struct unur_lobatto_table *Itable, const struct unur_gen *gen,
			   int print_Itable )
{
  FILE *LOG;
  int n;
  CHECK_NULL(Itable,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: subintervals for Lobatto integration: %d\n",gen->genid,
	  Itable->n_values - 1);
  for (n=0; print_Itable && n < Itable->n_values; n++) {
    fprintf(LOG,"%s:  [%3d] x = %g, u = %g\n",gen->genid,
	    n, Itable->values[n].x, Itable->values[n].u );
  }
} 
