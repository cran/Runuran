/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double *_unur_vector_new(int dim);
void _unur_vector_free(double *v);
double _unur_vector_norm(int dim, double *v);
double _unur_vector_scalar_product(int dim, double *v1, double *v2);
void _unur_vector_normalize(int dim, double *v);
