/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_DISTRIBUTIONS_SOURCE_H_SEEN
#define UNUR_DISTRIBUTIONS_SOURCE_H_SEEN
int _unur_stdgen_beta_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_beta_bb( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_bc( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b00( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b01( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b1prs( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_binv( UNUR_GEN *generator );
int _unur_stdgen_burr_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_cauchy_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_chi_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_chi_chru( UNUR_GEN *generator );
int _unur_stdgen_exponential_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_exponential_inv( UNUR_GEN *generator );
int _unur_stdgen_extremeI_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_extremeII_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_gamma_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_gamma_gll( UNUR_GEN *generator );
double _unur_stdgen_sample_gamma_gs( UNUR_GEN *generator );
double _unur_stdgen_sample_gamma_gd( UNUR_GEN *generator );
int _unur_stdgen_gig_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_gig_gigru( UNUR_GEN *generator );
int _unur_stdgen_laplace_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_logistic_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_lomax_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_normal_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_normal_bm( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_pol( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_quo( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_nquo( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_leva( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_kr( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_acr( UNUR_GEN *generator );
double _unur_stdgen_sample_normal_sum( UNUR_GEN *generator );
int _unur_stdgen_pareto_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_powerexponential_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_powerexponential_epd( UNUR_GEN *generator );
int _unur_stdgen_student_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_student_tpol( UNUR_GEN *generator );
double _unur_stdgen_sample_student_trouo( UNUR_GEN *generator );
int _unur_stdgen_slash_init( UNUR_PAR *parameters, UNUR_GEN *generator );
double _unur_stdgen_sample_slash_slash( UNUR_GEN *generator );
int _unur_stdgen_triangular_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_uniform_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_weibull_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_multinormal_init( UNUR_GEN *generator );
int _unur_stdgen_sample_multinormal_cholesky( UNUR_GEN *generator, double *X );
int _unur_stdgen_binomial_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_geometric_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_sample_geometric_inv( UNUR_GEN *generator );
int _unur_stdgen_hypergeometric_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_logarithmic_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_sample_logarithmic_lsk( UNUR_GEN *generator );
int _unur_stdgen_poisson_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pdtabl( UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pdac( UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pprsc( UNUR_GEN *generator );
int _unur_stdgen_zipf_init( UNUR_PAR *parameters, UNUR_GEN *generator );
int _unur_stdgen_sample_zipf_zet( UNUR_GEN *generator );
#define _unur_cstd_set_sampling_routine(gen,routine) \
   do { \
     if ((gen)==NULL) return UNUR_SUCCESS;            \
     (gen)->sample.cont = (routine);                  \
     ((struct unur_cstd_gen*)gen->datap)->sample_routine_name = #routine;   \
   } while (0)
#define _unur_dstd_set_sampling_routine(gen,routine) \
   do { \
     if ((gen)==NULL) return UNUR_SUCCESS;            \
     (gen)->sample.discr = (routine);                 \
     ((struct unur_dstd_gen*)gen->datap)->sample_routine_name = #routine;   \
   } while (0)
#endif  
