/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_cvec_new( int dim );
int unur_distr_cvec_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *pdf );
int unur_distr_cvec_set_dpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dpdf );
int unur_distr_cvec_set_pdpdf( UNUR_DISTR *distribution, UNUR_FUNCTD_CVEC *pdpdf );
UNUR_FUNCT_CVEC *unur_distr_cvec_get_pdf( const UNUR_DISTR *distribution );
UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dpdf( const UNUR_DISTR *distribution );
UNUR_FUNCTD_CVEC *unur_distr_cvec_get_pdpdf( const UNUR_DISTR *distribution );
double unur_distr_cvec_eval_pdf( const double *x, UNUR_DISTR *distribution );
int unur_distr_cvec_eval_dpdf( double *result, const double *x, UNUR_DISTR *distribution );
double unur_distr_cvec_eval_pdpdf( const double *x, int coord, UNUR_DISTR *distribution );
int unur_distr_cvec_set_logpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *logpdf );
int unur_distr_cvec_set_dlogpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dlogpdf );
int unur_distr_cvec_set_pdlogpdf( UNUR_DISTR *distribution, UNUR_FUNCTD_CVEC *pdlogpdf );
UNUR_FUNCT_CVEC *unur_distr_cvec_get_logpdf( const UNUR_DISTR *distribution );
UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dlogpdf( const UNUR_DISTR *distribution );
UNUR_FUNCTD_CVEC *unur_distr_cvec_get_pdlogpdf( const UNUR_DISTR *distribution );
double unur_distr_cvec_eval_logpdf( const double *x, UNUR_DISTR *distribution );
int unur_distr_cvec_eval_dlogpdf( double *result, const double *x, UNUR_DISTR *distribution );
double unur_distr_cvec_eval_pdlogpdf( const double *x, int coord, UNUR_DISTR *distribution );
int unur_distr_cvec_set_mean( UNUR_DISTR *distribution, const double *mean );
const double *unur_distr_cvec_get_mean( const UNUR_DISTR *distribution );
int unur_distr_cvec_set_covar( UNUR_DISTR *distribution, const double *covar );
int unur_distr_cvec_set_covar_inv( UNUR_DISTR *distribution, const double *covar_inv );
const double *unur_distr_cvec_get_covar( const UNUR_DISTR *distribution );
const double *unur_distr_cvec_get_cholesky( const UNUR_DISTR *distribution );
const double *unur_distr_cvec_get_covar_inv( UNUR_DISTR *distribution );
int unur_distr_cvec_set_rankcorr( UNUR_DISTR *distribution, const double *rankcorr );
const double *unur_distr_cvec_get_rankcorr( const UNUR_DISTR *distribution );
const double *unur_distr_cvec_get_rk_cholesky( const UNUR_DISTR *distribution );
int unur_distr_cvec_set_marginals( UNUR_DISTR *distribution, UNUR_DISTR *marginal );
int unur_distr_cvec_set_marginal_array( UNUR_DISTR *distribution, UNUR_DISTR **marginals );
int unur_distr_cvec_set_marginal_list( UNUR_DISTR *distribution, ... );
const UNUR_DISTR *unur_distr_cvec_get_marginal( const UNUR_DISTR *distribution, int n );
int unur_distr_cvec_set_pdfparams( UNUR_DISTR *distribution, const double *params, int n_params );
int unur_distr_cvec_get_pdfparams( const UNUR_DISTR *distribution, const double **params );
int unur_distr_cvec_set_pdfparams_vec( UNUR_DISTR *distribution, int par, const double *param_vec, int n_params );
int unur_distr_cvec_get_pdfparams_vec( const UNUR_DISTR *distribution, int par, const double **param_vecs );
int unur_distr_cvec_set_domain_rect( UNUR_DISTR *distribution, const double *lowerleft, const double *upperright );
int unur_distr_cvec_is_indomain( const double *x, const UNUR_DISTR *distribution );
int unur_distr_cvec_set_mode( UNUR_DISTR *distribution, const double *mode );
int unur_distr_cvec_upd_mode( UNUR_DISTR *distribution );
const double *unur_distr_cvec_get_mode( UNUR_DISTR *distribution );
int unur_distr_cvec_set_center( UNUR_DISTR *distribution, const double *center );
const double *unur_distr_cvec_get_center( UNUR_DISTR *distribution );
int unur_distr_cvec_set_pdfvol( UNUR_DISTR *distribution, double volume );
int unur_distr_cvec_upd_pdfvol( UNUR_DISTR *distribution );
double unur_distr_cvec_get_pdfvol( UNUR_DISTR *distribution );
