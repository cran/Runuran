/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_discr_new( void );
int unur_distr_discr_set_pv( UNUR_DISTR *distribution, const double *pv, int n_pv );
int unur_distr_discr_make_pv( UNUR_DISTR *distribution );
int unur_distr_discr_get_pv( const UNUR_DISTR *distribution, const double **pv );
int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
int unur_distr_discr_set_invcdf( UNUR_DISTR *distribution, UNUR_IFUNCT_DISCR *invcdf );
UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( const UNUR_DISTR *distribution );
UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( const UNUR_DISTR *distribution );
UNUR_IFUNCT_DISCR *unur_distr_discr_get_invcdf( const UNUR_DISTR *distribution );
double unur_distr_discr_eval_pv(int k, const UNUR_DISTR *distribution );
double unur_distr_discr_eval_pmf( int k, const UNUR_DISTR *distribution );
double unur_distr_discr_eval_cdf( int k, const UNUR_DISTR *distribution );
int unur_distr_discr_eval_invcdf( double u, const UNUR_DISTR *distribution );
int unur_distr_discr_set_pmfstr( UNUR_DISTR *distribution, const char *pmfstr );
int unur_distr_discr_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
char *unur_distr_discr_get_pmfstr( const UNUR_DISTR *distribution );
char *unur_distr_discr_get_cdfstr( const UNUR_DISTR *distribution );
int unur_distr_discr_set_pmfparams( UNUR_DISTR *distribution, const double *params, int n_params );
int unur_distr_discr_get_pmfparams( const UNUR_DISTR *distribution, const double **params );
int unur_distr_discr_set_domain( UNUR_DISTR *distribution, int left, int right );
int unur_distr_discr_get_domain( const UNUR_DISTR *distribution, int *left, int *right );
int unur_distr_discr_set_mode( UNUR_DISTR *distribution, int mode );
int unur_distr_discr_upd_mode( UNUR_DISTR *distribution );
int unur_distr_discr_get_mode( UNUR_DISTR *distribution );
int unur_distr_discr_set_pmfsum( UNUR_DISTR *distribution, double sum );
int unur_distr_discr_upd_pmfsum( UNUR_DISTR *distribution );
double unur_distr_discr_get_pmfsum( UNUR_DISTR *distribution );
