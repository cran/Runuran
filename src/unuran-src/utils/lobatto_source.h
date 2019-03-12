/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

typedef double UNUR_LOBATTO_FUNCT(double x, struct unur_gen *gen);
typedef double UNUR_LOBATTO_ERROR(struct unur_gen *gen, double delta, double x);
struct unur_lobatto_table;
double _unur_lobatto_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			       double x, double h, double tol, UNUR_LOBATTO_ERROR uerror);
struct unur_lobatto_table *
_unur_lobatto_init (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		    double left, double center, double right,
		    double tol, UNUR_LOBATTO_ERROR uerror, int size);
int _unur_lobatto_find_linear (struct unur_lobatto_table *Itable, double x);
double _unur_lobatto_eval_diff (struct unur_lobatto_table *Itable, double x, double h, double *fx);
double _unur_lobatto_eval_CDF (struct unur_lobatto_table *Itable, double x);
double _unur_lobatto_integral (struct unur_lobatto_table *Itable );
void _unur_lobatto_free (struct unur_lobatto_table **Itable);
void _unur_lobatto_debug_table (struct unur_lobatto_table *Itable,
				const struct unur_gen *gen, int print_Itable );
int _unur_lobatto_size_table (struct unur_lobatto_table *Itable);
