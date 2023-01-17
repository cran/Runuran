/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNURAN_TESTS_H_SEEN
#define UNURAN_TESTS_H_SEEN
#define UNUR_TEST_ALL      (~0u)     
#define UNUR_TEST_TIME     0x001u    
#define UNUR_TEST_N_URNG   0x002u    
#define UNUR_TEST_N_PDF    0x004u    
#define UNUR_TEST_CHI2     0x008u    
#define UNUR_TEST_SAMPLE   0x010u    
void unur_run_tests( UNUR_PAR *parameters, unsigned tests, FILE *out );
void unur_test_printsample( UNUR_GEN *generator, int n_rows, int n_cols, FILE *out );
UNUR_GEN *unur_test_timing( UNUR_PAR *parameters, int log10_samplesize, 
			    double *time_setup, double *time_sample,
			    int verbosity, FILE *out );
double unur_test_timing_R( UNUR_PAR *parameters, const char *distrstr, const char *methodstr,
			   double log10_samplesize, double *time_setup, double *time_marginal );
double unur_test_timing_uniform( const UNUR_PAR *parameters, int log10_samplesize );
double unur_test_timing_exponential( const UNUR_PAR *parameters, int log10_samplesize );
double unur_test_timing_total( const UNUR_PAR *parameters, int samplesize, double avg_duration );
int unur_test_count_urn( UNUR_GEN *generator, int samplesize, int verbosity, FILE *out );
int unur_test_count_pdf( UNUR_GEN *generator, int samplesize, int verbosity, FILE *out );
int unur_test_par_count_pdf( UNUR_PAR *parameters, int samplesize, int verbosity, FILE *out );
double unur_test_chi2( UNUR_GEN *generator, int intervals, int samplesize, int classmin,
		       int verbosity, FILE *out );
int unur_test_moments( UNUR_GEN *generator, double *moments, int n_moments, int samplesize,
		       int verbosity, FILE *out );
double unur_test_correlation( UNUR_GEN *generator1, UNUR_GEN *generator2,
			      int samplesize, int verbosity, FILE *out );
int unur_test_quartiles( UNUR_GEN *generator,
			 double *q0, double *q1, double *q2, double *q3, double *q4, 
			 int samplesize, int verbosity, FILE *out );
double unur_test_u_error( const UNUR_GEN *generator, 
			  double *max_error, double *MAE, double threshold,
			  int samplesize, int randomized, int testtails,
			  int verbosity, FILE *out );
int unur_test_cvec_rankcorr( double *rc, UNUR_GEN *gen, int samplesize, int verbose, FILE *out );
#endif  
