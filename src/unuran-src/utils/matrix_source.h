/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

int _unur_matrix_transform_diagonal (int dim, const double *M, const double *D, double *res);
int _unur_matrix_cholesky_decomposition (int dim, const double *S, double *L );
int _unur_matrix_invert_matrix (int dim, const double *A, double *Ainv, double *det );
double _unur_matrix_determinant ( int dim, const double *A );
int _unur_matrix_multiplication(int dim, const double *A, const double *B, double *AB);
double _unur_matrix_qf (int dim, double *x, double *A );
int _unur_matrix_eigensystem (int dim, const double *M, double *values, double *vectors );
void _unur_matrix_print_vector ( int dim, const double *vec, const char *info,
				 FILE *LOG, const char *genid, const char *indent );
void _unur_matrix_print_matrix ( int dim, const double *mat, const char *info,
				 FILE *LOG, const char *genid, const char *indent );
