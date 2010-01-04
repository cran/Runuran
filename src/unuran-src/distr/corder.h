/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_DISTR *unur_distr_corder_new( const UNUR_DISTR *distribution, int n, int k );
const UNUR_DISTR *unur_distr_corder_get_distribution( const UNUR_DISTR *distribution );
int unur_distr_corder_set_rank( UNUR_DISTR *distribution, int n, int k );
int unur_distr_corder_get_rank( const UNUR_DISTR *distribution, int *n, int *k );
#define unur_distr_corder_get_pdf(distr)   unur_distr_cont_get_pdf((distr))
#define unur_distr_corder_get_dpdf(distr)  unur_distr_cont_get_dpdf((distr))
#define unur_distr_corder_get_cdf(distr)   unur_distr_cont_get_cdf((distr))
#define unur_distr_corder_eval_pdf(x,distr)  unur_distr_cont_eval_pdf((x),(distr))
#define unur_distr_corder_eval_dpdf(x,distr) unur_distr_cont_eval_dpdf((x),(distr))
#define unur_distr_corder_eval_cdf(x,distr)  unur_distr_cont_eval_cdf((x),(distr))
#define unur_distr_corder_set_pdfparams(distr,params,n)  unur_distr_cont_set_pdfparams((distr),(params),(n))
#define unur_distr_corder_get_pdfparams(distr,params)  unur_distr_cont_get_pdfparams((distr),(params))
#define unur_distr_corder_set_domain(distr,left,right)  unur_distr_cont_set_domain((distr),(left),(right))
#define unur_distr_corder_get_domain(distr,left,right)  unur_distr_cont_get_domain((distr),(left),(right))
#define unur_distr_corder_get_truncated(distr,left,right)  unur_distr_cont_get_truncated((distr),(left),(right))
#define unur_distr_corder_set_mode(distr,mode)   unur_distr_cont_set_mode((distr),(mode))
#define unur_distr_corder_upd_mode(distr)   unur_distr_cont_upd_mode((distr))
#define unur_distr_corder_get_mode(distr)   unur_distr_cont_get_mode((distr))
#define unur_distr_corder_set_pdfarea(distr,area)   unur_distr_cont_set_pdfarea((distr),(area))
#define unur_distr_corder_upd_pdfarea(distr)   unur_distr_cont_upd_pdfarea((distr))
#define unur_distr_corder_get_pdfarea(distr)   unur_distr_cont_get_pdfarea((distr))
