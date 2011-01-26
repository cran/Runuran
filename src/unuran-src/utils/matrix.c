/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "matrix_source.h"
static int _unur_matrix_swap_rows (int dim, double *A, int i, int j);
static int _unur_matrix_permutation_swap (int dim, int *p, int i, int j);
static int _unur_matrix_LU_decomp (int dim, double *A, int *P, int *signum);
static int _unur_matrix_backsubstitution_dtrsv(int dim, double *LU, double *X);
static int _unur_matrix_forwardsubstitution_dtrsv(int dim, double *LU, double *X);
static int _unur_matrix_LU_invert (int dim, double *LU, int *p, double *inverse);
int
_unur_matrix_transform_diagonal (int dim, const double *M, const double *D, double *res)
{
#define idx(a,b) ((a)*dim+(b))
  int i,j,k;
  double sum;
  CHECK_NULL(M,UNUR_ERR_NULL);
  CHECK_NULL(D,UNUR_ERR_NULL);
  CHECK_NULL(res,UNUR_ERR_NULL);
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++) {
      for (sum=0., k=0; k<dim; k++)
	sum += D[k] * M[idx(k,i)] * M[idx(k,j)];
      res[idx(i,j)] = sum;
    }
  return UNUR_SUCCESS;
#undef idx
} 
int 
_unur_matrix_swap_rows (int dim, double *A, int i, int j)
{
  CHECK_NULL(A,UNUR_ERR_NULL);
  if (i != j)
  {
    double *row1 = A + i * dim;
    double *row2 = A + j * dim;
    int k;
    for (k = 0; k < dim; k++)
    {
       double tmp = row1[k] ;
       row1[k] = row2[k] ;
       row2[k] = tmp ;
    }
  }
  return UNUR_SUCCESS;
} 
int 
_unur_matrix_permutation_swap (int dim ATTRIBUTE__UNUSED, int *p, int i, int j)
{
  if (i != j)
    {
      int tmp = p[i];
      p[i] = p[j];
      p[j] = tmp;
    }
  return UNUR_SUCCESS;
} 
int 
_unur_matrix_LU_decomp (int dim, double *A, int *p, int *signum)
{
#define idx(a,b) ((a)*dim+(b))
  int i, j, k;
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(p,UNUR_ERR_NULL);
  CHECK_NULL(signum,UNUR_ERR_NULL);
  *signum = 1;
  for(i=0;i<dim;i++) p[i]=i;
  for (j = 0; j < dim - 1; j++){
    double ajj, max = fabs (A[idx(j,j)]);
    int i_pivot = j;
    for (i = j + 1; i < dim; i++){
      double aij = fabs (A[idx(i,j)]);
      if (aij > max){
	max = aij;
	i_pivot = i;
      }
    }
    if (i_pivot != j){
      _unur_matrix_swap_rows (dim, A, j, i_pivot);
      _unur_matrix_permutation_swap (dim, p, j, i_pivot);
      *signum = -(*signum);
    }
    ajj = A[idx(j,j)];
    if ( !_unur_iszero(ajj) ){
      for (i = j + 1; i < dim; i++){
	double aij = A[idx(i,j)] / ajj;
	A[idx(i,j)] = aij;
	for (k = j + 1; k < dim; k++){
	  double aik = A[idx(i,k)];
	  double ajk = A[idx(j,k)];
	  A[idx(i,k)]= aik - aij * ajk;
	}
      }
    }
  }
  return UNUR_SUCCESS;
#undef idx
} 
int 
_unur_matrix_backsubstitution_dtrsv(int dim, double *LU, double *X)
{
#define idx(a,b) ((a)*dim+(b))
  int ix,jx,i,j;
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(X,UNUR_ERR_NULL);
  ix = (dim - 1);
  X[ix] = X[ix] / LU[idx(ix,ix)];
  ix--;
  for (i = dim - 1; i > 0 && i--;) {
    double tmp = X[ix];
    jx = ix + 1;
    for (j = i + 1; j < dim; j++) {
      tmp -= LU[idx(i,j)] * X[jx];
      jx ++;
    }
    X[ix] = tmp / LU[idx(i,i)];
    ix --;
  }
  return UNUR_SUCCESS;
#undef idx
} 
int 
_unur_matrix_forwardsubstitution_dtrsv(int dim, double *LU, double *X)
{ 
#define idx(a,b) ((a)*dim+(b))
  int ix,jx,i,j;
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(X,UNUR_ERR_NULL);
  ix = 0;
  ix++;
  for (i = 1; i < dim; i++) {
    double tmp = X[ix];
    jx = 0;
    for (j = 0; j < i; j++) {
      tmp -= LU[idx(i,j)] * X[jx];
      jx += 1;
    }
    X[ix] = tmp;
    ix ++;
  }
  return UNUR_SUCCESS;
#undef idx
} 
int 
_unur_matrix_LU_invert (int dim, double *LU, int *p, double *inverse)
{ 
#define idx(a,b) ((a)*dim+(b))
  double *vector;
  int i,j;
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(p,UNUR_ERR_NULL);
  CHECK_NULL(inverse,UNUR_ERR_NULL);
  vector = _unur_xmalloc(dim*sizeof(double));
  for (i = 0; i < dim; i++){
    for(j=0;j<dim;j++) vector[j] = 0.;
    vector[i] = 1.;
    _unur_matrix_forwardsubstitution_dtrsv (dim, LU, vector);
    _unur_matrix_backsubstitution_dtrsv (dim, LU, vector);
    for(j=0;j<dim;j++){
      inverse[idx(j,p[i])] = vector[j]; 
    }
  }
  free(vector);
  return UNUR_SUCCESS;
#undef idx
} 
int 
_unur_matrix_invert_matrix(int dim, const double *A, double *Ainv, double *det)
{ 
#define idx(a,b) ((a)*dim+(b))
  int *p, s, i;
  double *LU;             
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(Ainv,UNUR_ERR_NULL);
  CHECK_NULL(det,UNUR_ERR_NULL);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }
  p = _unur_xmalloc(dim*sizeof(int));
  LU = _unur_xmalloc(dim*dim*sizeof(double));
  memcpy(LU, A, dim*dim*sizeof(double));
  _unur_matrix_LU_decomp(dim, LU, p, &s);
  *det = s;
  for(i=0;i<dim;i++)
    *det *= LU[idx(i,i)];
  _unur_matrix_LU_invert(dim, LU, p, Ainv);   
  free(LU);
  free(p);
  return UNUR_SUCCESS;
#undef idx
} 
int 
_unur_matrix_multiplication(int dim, const double *A, const double *B, double *AB)
{ 
#define idx(a,b) ((a)*dim+(b))
  int i, j, k;
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(B,UNUR_ERR_NULL);
  CHECK_NULL(AB,UNUR_ERR_NULL);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }
  for(i=0;i<dim;i++) 
  for(j=0;j<dim;j++) {
    AB[idx(i,j)]=0.;
    for(k=0;k<dim;k++) {
      AB[idx(i,j)] += A[idx(i,k)]*B[idx(k,j)];
    }
  }
  return UNUR_SUCCESS;
#undef idx
} 
double
_unur_matrix_determinant ( int dim, const double *A )
{
#define idx(a,b) ((a)*dim+(b))
  int *p, s, i;
  double *LU;     
  double det;
  CHECK_NULL(A,  INFINITY);
  if (dim==1) return A[0];
  p = _unur_xmalloc(dim*sizeof(int));
  LU = _unur_xmalloc(dim*dim*sizeof(double));
  memcpy(LU, A, dim*dim*sizeof(double));
  _unur_matrix_LU_decomp(dim, LU, p, &s);
  det = s;
  for(i=0;i<dim;i++)
    det *= LU[idx(i,i)];
  free(LU);
  free(p);
  return det;
#undef idx
} 
double 
_unur_matrix_qf (int dim, double *x, double *A)
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;
  double sum,outersum;
  CHECK_NULL(x,INFINITY);
  CHECK_NULL(A,INFINITY);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return INFINITY;
  }
  outersum=0.;
  for(i=0;i<dim;i++){
    sum=0.;
    for(j=0;j<dim;j++)
      sum+=A[idx(i,j)]*x[j];
    outersum+=sum*x[i];
  }
  return outersum;
#undef idx
} 
int
_unur_matrix_cholesky_decomposition (int dim, const double *S, double *L )
{ 
#define idx(a,b) ((a)*dim+(b))
  int i,j,k;
  double sum1,sum2;
  CHECK_NULL(S,UNUR_ERR_NULL);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }
  L[idx(0,0)] = sqrt( S[idx(0,0)] );
  for(j=1; j<dim; j++) {
    L[idx(j,0)] = S[idx(j,0)] / L[idx(0,0)];
    sum1 = L[idx(j,0)] * L[idx(j,0)];
    for(k=1; k<j; k++) {
      sum2 = 0.;
      for(i=0; i<k; i++)
	sum2 += L[idx(j,i)] * L[idx(k,i)];
      L[idx(j,k)] = (S[idx(j,k)] - sum2) / L[idx(k,k)];
      sum1 += L[idx(j,k)] * L[idx(j,k)];
    }
    if (!(S[idx(j,j)] > sum1)) {
      return UNUR_FAILURE;
    }
    L[idx(j,j)] = sqrt( S[idx(j,j)] - sum1 );
  }
  for(j=0; j<dim; j++)
    for(k=j+1; k<dim; k++)
      L[idx(j,k)]=0.;
  return UNUR_SUCCESS;
#undef idx
} 
void
_unur_matrix_print_vector ( int dim, const double *vec, const char *info,
			    FILE *LOG, const char *genid, const char *indent )
{
  int i;
  if (vec) {
    fprintf(LOG,"%s: %s\n", genid, info );
    fprintf(LOG,"%s: %s( %g", genid, indent, vec[0]);
    for (i=1; i<dim; i++) 
      fprintf(LOG,", %g", vec[i]);
    fprintf(LOG," )\n");
  }
  else {
    fprintf(LOG,"%s: %s [unknown]\n", genid, info );
  }
  fprintf(LOG,"%s:\n",genid);
} 
void
_unur_matrix_print_matrix ( int dim, const double *mat, const char *info,
			   FILE *LOG, const char *genid, const char *indent )
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;
  if (mat) {
    fprintf(LOG,"%s: %s\n", genid, info); 
    for (i=0; i<dim; i++) {
      fprintf(LOG,"%s: %s(% e", genid, indent, mat[idx(i,0)]);
      for (j=1; j<dim; j++)
	fprintf(LOG,",% e",mat[idx(i,j)]);
      fprintf(LOG," )\n");
    }
  }
  else {
    fprintf(LOG,"%s: %s [unknown]\n", genid, info );
  }
  fprintf(LOG,"%s:\n",genid);
#undef idx
} 
