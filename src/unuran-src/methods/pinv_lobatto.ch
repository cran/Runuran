/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#define W1 (0.17267316464601146)   
#define W2 (1.-W1)
double
_unur_pinv_lobatto5 (struct unur_gen *gen, double x, double h)
{ 
  return (9*(PDF(x)+PDF(x+h))+49.*(PDF(x+h*W1)+PDF(x+h*W2))+64*PDF(x+h/2.))*h/180.;
} 
double
_unur_pinv_adaptivelobatto5 (struct unur_gen *gen, double x, double h, double tol,
			     struct unur_pinv_CDFtable *CDFtable)
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
  fl = PDF(x);
  fc = PDF(x+h/2.);
  fr = PDF(x+h);
  int1 = (9*(fl+fr)+49.*(PDF(x+h*W1)+PDF(x+h*W2))+64*fc)*h/180.;
  int2 = _unur_pinv_adaptivelobatto5_rec(gen,x,h,tol,int1,fl,fc,fr,&W_accuracy,CDFtable);
  if (W_accuracy)
    _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
		  "numeric integration did not reach full accuracy");
  return int2;
} 
double
_unur_pinv_adaptivelobatto5_rec (struct unur_gen *gen, double x, double h, double tol,
				 double int1, double fl, double fc, double fr,
				 int *W_accuracy,
				 struct unur_pinv_CDFtable *CDFtable)
{
  double flc, frc;    
  double int2;        
  double intl, intr;  
  flc = PDF(x+h/4);
  frc = PDF(x+3*h/4);
  intl = (9*(fl+fc)+49.*(PDF(x+h*W1*0.5)+PDF(x+h*W2*0.5))+64*flc)*h/360.;
  intr = (9*(fc+fr)+49.*(PDF(x+h*(0.5+W1*0.5))+PDF(x+h*(0.5+W2*0.5)))+64*frc)*h/360.;
  int2 = intl + intr;
  if (fabs(int1-int2) < tol)
    ;
  else {
    if (_unur_FP_equal(x+h/2.,x)) {
      *W_accuracy = TRUE;
    }
    else {
      return ( _unur_pinv_adaptivelobatto5_rec(gen,x,    h/2,tol/1.,intl,fl,flc,fc,
					       W_accuracy,CDFtable) +
	       _unur_pinv_adaptivelobatto5_rec(gen,x+h/2,h/2,tol/1.,intr,fc,frc,fr,
					       W_accuracy,CDFtable) );
    }
  }
#ifdef PINV_USE_CDFTABLE
  if (CDFtable) {
    _unur_pinv_CDFtable_append(CDFtable, x+h/2., intl);
    _unur_pinv_CDFtable_append(CDFtable, x+h, intr);
  }
#endif
  return int2;
} 
#undef W1
#undef W2
double
_unur_pinv_Udiff_lobatto (struct unur_gen *gen, double x, double h, double utol)
{
#ifdef PINV_USE_CDFTABLE
  struct unur_pinv_CDFtable *CDFtable = GEN->CDFtable; 
  int cur;                    
  double x1;                  
  double Udiff;               
#endif
  if (!_unur_isfinite(x+h)) {
    return INFINITY;
  }
#ifdef PINV_USE_CDFTABLE
  if (CDFtable == NULL) {
    if (x < GEN->bleft || x+h > GEN->bright)
      return _unur_pinv_lobatto5(gen, x, h);
    else
      return _unur_pinv_adaptivelobatto5(gen, x, h, utol, NULL);
  }
  cur = CDFtable->cur_iv;
  while (cur < CDFtable->n_values &&
	 CDFtable->values[cur].x < x)
    ++cur;
  if (cur >= CDFtable->n_values) {
    return _unur_pinv_adaptivelobatto5(gen, x, h, utol, NULL);
  }
  x1 = CDFtable->values[cur].x;
  ++cur;
  if (cur >= CDFtable->n_values ||
      CDFtable->values[cur].x > x+h) {
    return _unur_pinv_lobatto5(gen, x, h);
  }
  Udiff = _unur_pinv_lobatto5(gen, x, x1 - x);
  do {
    Udiff += CDFtable->values[cur].u;
    x1 = CDFtable->values[cur].x;
    ++cur;
  } while (cur < CDFtable->n_values
	   && CDFtable->values[cur].x <= x+h);
  if (x+h < GEN->bright && cur >= CDFtable->n_values) {
    Udiff += _unur_pinv_adaptivelobatto5(gen, x1, x+h, utol, NULL);
  }
  else {
    Udiff += _unur_pinv_lobatto5(gen, x1, x+h - x1);
  }
  return Udiff;
#else
#ifdef PINV_USE_SIMPLE_LOBATTO
  return _unur_pinv_lobatto5(gen, x, h);
#else
  if (x < GEN->bleft || x+h > GEN->bright)
    return _unur_pinv_lobatto5(gen, x, h);
  else
    return _unur_pinv_adaptivelobatto5(gen, x, h, utol, NULL);
#endif
#endif 
} 
#ifdef PINV_USE_CDFTABLE
struct unur_pinv_CDFtable *
_unur_pinv_CDFtable_create (int size)
{
  struct unur_pinv_CDFtable *table;
  if (size<=0)
    return NULL;
  table = _unur_xmalloc( sizeof(struct unur_pinv_CDFtable) );
  table->values = _unur_xmalloc(size * sizeof(struct unur_pinv_CDFvalues) ); 
  table->size = size;
  table->n_values = 0;
  table->cur_iv = 0;
  return table;
} 
int
_unur_pinv_CDFtable_append (struct unur_pinv_CDFtable *table, double x, double u)
{
  if (table==NULL) 
    return UNUR_ERR_NULL;
  if (table->n_values >= table->size - 1)
    return UNUR_ERR_GENERIC;
  table->values[table->n_values].x = x;
  table->values[table->n_values].u = u;
  ++(table->n_values);
  return UNUR_SUCCESS;
} 
void 
_unur_pinv_CDFtable_resize (struct unur_pinv_CDFtable **table)
{
  if (*table) {
    *table = _unur_xrealloc(*table, (*table)->n_values * sizeof(double));
    (*table)->size = (*table)->n_values;
  }
} 
void
_unur_pinv_CDFtable_free (struct unur_pinv_CDFtable **table)
{
  if (*table) {
    free ((*table)->values);
    free (*table);
    *table = NULL;
  }
} 
#endif 
