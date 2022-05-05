/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
double *
_unur_vector_new(int dim)
{
  int i;
  double *v;
  v = _unur_xmalloc(dim*sizeof(double));
  for (i=0; i<dim; i++) v[i] = 0.;
  return v;
} 
void 
_unur_vector_free(double *v)
{
  if (v) free(v);
} 
double 
_unur_vector_norm(int dim, double *v)
{
  int i;
  double vsum;
  double vmax;
  double p;
  if (v==NULL) return 0.; 
  vmax = 0.;
  for (i=0; i<dim; i++) {
    if (vmax < fabs(v[i])) vmax = fabs(v[i]); 
  }
  if (vmax <= 0) return 0.;
  vsum = 0.;
  for (i=0; i<dim; i++) {
    p=v[i]/vmax;
    vsum += p*p;
  }
  return vmax * sqrt(vsum);
} 
void 
_unur_vector_normalize(int dim, double *v)
{
  int i;
  double norm;
  if (v==NULL) return; 
  norm = _unur_vector_norm(dim, v);
  for (i=0; i<dim; i++)  v[i] /= norm;
} 
double 
_unur_vector_scalar_product(int dim, double *v1, double *v2)
{
  int i;
  double scalar_product;
  if (v1==NULL || v2==NULL) return 0.; 
  scalar_product = 0.;
  for (i=0; i<dim; i++) {
    scalar_product += v1[i]*v2[i];
  }
  return scalar_product;
} 
