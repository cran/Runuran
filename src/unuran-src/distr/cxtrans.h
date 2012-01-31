/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_cxtrans_new( const UNUR_DISTR *distribution );
const UNUR_DISTR *unur_distr_cxtrans_get_distribution( const UNUR_DISTR *distribution );
int unur_distr_cxtrans_set_alpha( UNUR_DISTR *distribution, double alpha );
int unur_distr_cxtrans_set_rescale( UNUR_DISTR *distribution, double mu, double sigma );
double unur_distr_cxtrans_get_alpha( const UNUR_DISTR *distribution );
double unur_distr_cxtrans_get_mu( const UNUR_DISTR *distribution );
double unur_distr_cxtrans_get_sigma( const UNUR_DISTR *distribution );
int unur_distr_cxtrans_set_logpdfpole( UNUR_DISTR *distribution, double logpdfpole, double dlogpdfpole );
#define unur_distr_cxtrans_get_pdf(distr)   unur_distr_cont_get_pdf((distr))
#define unur_distr_cxtrans_get_dpdf(distr)  unur_distr_cont_get_dpdf((distr))
#define unur_distr_cxtrans_get_cdf(distr)   unur_distr_cont_get_cdf((distr))
#define unur_distr_cxtrans_eval_pdf(x,distr)  unur_distr_cont_eval_pdf((x),(distr))
#define unur_distr_cxtrans_eval_dpdf(x,distr) unur_distr_cont_eval_dpdf((x),(distr))
#define unur_distr_cxtrans_eval_cdf(x,distr)  unur_distr_cont_eval_cdf((x),(distr))
int unur_distr_cxtrans_set_domain( UNUR_DISTR *distribution, double left, double right );
#define unur_distr_cxtrans_get_domain(distr,left,right)  unur_distr_cont_get_domain((distr),(left),(right))
#define unur_distr_cxtrans_get_truncated(distr,left,right)  unur_distr_cont_get_truncated((distr),(left),(right))
#define unur_distr_cxtrans_set_mode(distr,mode)   unur_distr_cont_set_mode((distr),(mode))
#define unur_distr_cxtrans_upd_mode(distr)   unur_distr_cont_upd_mode((distr))
#define unur_distr_cxtrans_get_mode(distr)   unur_distr_cont_get_mode((distr))
#define unur_distr_cxtrans_set_pdfarea(distr,area)   unur_distr_cont_set_pdfarea((distr),(area))
#define unur_distr_cxtrans_upd_pdfarea(distr)   unur_distr_cont_upd_pdfarea((distr))
#define unur_distr_cxtrans_get_pdfarea(distr)   unur_distr_cont_get_pdfarea((distr))
