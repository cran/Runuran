/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double _unur_util_find_max( struct unur_funct_generic fs,
                            double interval_min, 
		            double interval_max,
		            double max_guess );
double _unur_util_brent(struct unur_funct_generic fs,
                        double a, double b, double c, double tol);
