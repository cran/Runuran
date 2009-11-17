/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNURAN_DISTRIBUTIONS_H_SEEN
#define UNURAN_DISTRIBUTIONS_H_SEEN
#include <distributions/unur_stddistr.h>
UNUR_DISTR *unur_distr_beta(const double *params, int n_params);
UNUR_DISTR *unur_distr_burr(const double *params, int n_params);
UNUR_DISTR *unur_distr_cauchy(const double *params, int n_params);
UNUR_DISTR *unur_distr_chi(const double *params, int n_params);
UNUR_DISTR *unur_distr_chisquare(const double *params, int n_params);
UNUR_DISTR *unur_distr_exponential(const double *params, int n_params);
UNUR_DISTR *unur_distr_extremeI(const double *params, int n_params);
UNUR_DISTR *unur_distr_extremeII(const double *params, int n_params);
UNUR_DISTR *unur_distr_F(const double *params, int n_params);
UNUR_DISTR *unur_distr_gamma(const double *params, int n_params);
UNUR_DISTR *unur_distr_ghyp(const double *params, int n_params);
UNUR_DISTR *unur_distr_gig(const double *params, int n_params);
UNUR_DISTR *unur_distr_gig2(const double *params, int n_params);
UNUR_DISTR *unur_distr_hyperbolic(const double *params, int n_params);
UNUR_DISTR *unur_distr_ig(const double *params, int n_params);
UNUR_DISTR *unur_distr_laplace(const double *params, int n_params);
UNUR_DISTR *unur_distr_logistic(const double *params, int n_params);
UNUR_DISTR *unur_distr_lognormal(const double *params, int n_params);
UNUR_DISTR *unur_distr_lomax(const double *params, int n_params);
UNUR_DISTR *unur_distr_normal( const double *params, int n_params );
UNUR_DISTR *unur_distr_pareto( const double *params, int n_params );
UNUR_DISTR *unur_distr_powerexponential(const double *params, int n_params);
UNUR_DISTR *unur_distr_rayleigh(const double *params, int n_params);
UNUR_DISTR *unur_distr_slash(const double *params, int n_params);
UNUR_DISTR *unur_distr_student(const double *params, int n_params);
UNUR_DISTR *unur_distr_triangular(const double *params, int n_params);
UNUR_DISTR *unur_distr_uniform(const double *params, int n_params);
UNUR_DISTR *unur_distr_weibull(const double *params, int n_params);
UNUR_DISTR *unur_distr_multinormal(int dim, const double *mean, const double *covar);
UNUR_DISTR *unur_distr_multicauchy(int dim, const double *mean, const double *covar);
UNUR_DISTR *unur_distr_multistudent(int dim, double nu, const double *mean, const double *covar);
UNUR_DISTR *unur_distr_multiexponential(int dim, const double *sigma, const double *theta);
UNUR_DISTR *unur_distr_copula(int dim, const double *rankcorr);
UNUR_DISTR *unur_distr_correlation( int n );
UNUR_DISTR *unur_distr_binomial(const double *params, int n_params);
UNUR_DISTR *unur_distr_geometric(const double *params, int n_params);
UNUR_DISTR *unur_distr_hypergeometric(const double *params, int n_params);
UNUR_DISTR *unur_distr_logarithmic(const double *params, int n_params);
UNUR_DISTR *unur_distr_negativebinomial(const double *params, int n_params);
UNUR_DISTR *unur_distr_poisson(const double *params, int n_params);
UNUR_DISTR *unur_distr_zipf(const double *params, int n_params);
#endif  
