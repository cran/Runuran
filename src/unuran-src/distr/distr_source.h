/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_DISTR_SOURCE_H_SEEN
#define UNUR_DISTR_SOURCE_H_SEEN
#define UNUR_DISTR_SET_MASK_ESSENTIAL 0xffff0000u
#define UNUR_DISTR_SET_DOMAIN         0x00010000u
#define UNUR_DISTR_SET_DOMAINBOUNDED  0x00020000u 
#define UNUR_DISTR_SET_STDDOMAIN      0x00040000u 
#define UNUR_DISTR_SET_TRUNCATED      0x00080000u 
#define UNUR_DISTR_SET_MEAN           0x01000000u 
#define UNUR_DISTR_SET_COVAR          0x02000000u 
#define UNUR_DISTR_SET_COVAR_INV      0x04000000u 
#define UNUR_DISTR_SET_COVAR_IDENT    0x40000000u 
#define UNUR_DISTR_SET_CHOLESKY       0x08000000u 
#define UNUR_DISTR_SET_RANKCORR       0x10000000u 
#define UNUR_DISTR_SET_RK_CHOLESKY    0x20000000u 
#define UNUR_DISTR_SET_STDMARGINAL    0x00100000u 
#define UNUR_DISTR_SET_MARGINAL       0x00200000u 
#define UNUR_DISTR_SET_GENERIC        0x00080000u 
#define UNUR_DISTR_SET_MASK_DERIVED   0x0000ffffu
#define UNUR_DISTR_SET_MODE           0x00000001u
#define UNUR_DISTR_SET_MODE_APPROX    0x00000020u 
#define UNUR_DISTR_SET_CENTER         0x00000002u
#define UNUR_DISTR_SET_CENTER_APPROX  0x00000040u 
#define UNUR_DISTR_SET_PDFAREA        0x00000004u
#define UNUR_DISTR_SET_PMFSUM         0x00000008u
#define UNUR_DISTR_SET_PDFVOLUME      0x00000010u
#define _unur_cont_PDF(x,distr)     ((*((distr)->data.cont.pdf)) ((x),(distr)))
#define _unur_cont_dPDF(x,distr)    ((*((distr)->data.cont.dpdf))((x),(distr)))
#define _unur_cont_logPDF(x,distr)  ((*((distr)->data.cont.logpdf)) ((x),(distr)))
#define _unur_cont_dlogPDF(x,distr) ((*((distr)->data.cont.dlogpdf))((x),(distr)))
#define _unur_cont_CDF(x,distr)     ((*((distr)->data.cont.cdf)) ((x),(distr)))
#define _unur_cont_logCDF(x,distr)  ((*((distr)->data.cont.logcdf)) ((x),(distr)))
#define _unur_cont_invCDF(u,distr)  ((*((distr)->data.cont.invcdf)) ((u),(distr)))
#define _unur_cont_HR(x,distr)      ((*((distr)->data.cont.hr))  ((x),(distr)))
#define _unur_discr_PMF(x,distr)    ((*((distr)->data.discr.pmf))((x),(distr)))
#define _unur_discr_CDF(x,distr)    ((*((distr)->data.discr.cdf))((x),(distr)))
#define _unur_discr_invCDF(u,distr) ((int) (*((distr)->data.discr.invcdf)) ((u),(distr)))
double _unur_cvec_PDF(const double *x, struct unur_distr *distr);
int _unur_cvec_dPDF(double *result, const double *x, struct unur_distr *distr);
double _unur_cvec_pdPDF(const double *x, int coord, struct unur_distr *distr);
double _unur_cvec_logPDF(const double *x, struct unur_distr *distr);
int _unur_cvec_dlogPDF(double *result, const double *x, struct unur_distr *distr);
double _unur_cvec_pdlogPDF(const double *x, int coord, struct unur_distr *distr);
#define _unur_cont_have_logPDF(distr)  (((distr)->data.cont.logpdf==NULL)?FALSE:TRUE)
#define _unur_cont_have_dlogPDF(distr) (((distr)->data.cont.dlogpdf==NULL)?FALSE:TRUE)
double _unur_distr_cont_eval_pdf_from_logpdf( double x, const struct unur_distr *distr );
double _unur_distr_cont_eval_dpdf_from_dlogpdf( double x, const struct unur_distr *distr );
double _unur_distr_cont_eval_cdf_from_logcdf( double x, const struct unur_distr *distr );
double _unur_distr_cvec_eval_pdf_from_logpdf( const double *x, struct unur_distr *distr );
int _unur_distr_cvec_eval_dpdf_from_dlogpdf( double *result, const double *x, struct unur_distr *distr );
double _unur_distr_cvec_eval_pdpdf_from_pdlogpdf( const double *x, int coord, struct unur_distr *distr );
struct unur_distr *_unur_distr_generic_new( void );
struct unur_distr *_unur_distr_cemp_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_cont_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_matr_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_cvec_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_cvemp_clone( const struct unur_distr *distr );
struct unur_distr *_unur_distr_discr_clone( const struct unur_distr *distr );
#define _unur_distr_clone(distr)    ((distr)->clone(distr))
#define _unur_distr_free(distr)    do {if (distr) (distr)->destroy(distr);} while(0)
#ifdef UNUR_ENABLE_LOGGING
void _unur_distr_cont_debug( const UNUR_DISTR *distribution, const char *genid );
void _unur_distr_corder_debug( const UNUR_DISTR *order_statistics, const char *genid );
void _unur_distr_cxtrans_debug( const UNUR_DISTR *cxtrans, const char *genid );
void _unur_distr_cemp_debug( const UNUR_DISTR *distribution, const char *genid, unsigned printvector );
void _unur_distr_matr_debug( const UNUR_DISTR *distribution, const char *genid );
void _unur_distr_cvec_debug( const UNUR_DISTR *distribution, const char *genid );
void _unur_distr_condi_debug( const UNUR_DISTR *distribution, const char *genid );
void _unur_distr_cvemp_debug( const UNUR_DISTR *distribution, const char *genid, unsigned printvector );
void _unur_distr_discr_debug( const UNUR_DISTR *distribution, const char *genid, unsigned printvector );
#endif
#ifdef UNUR_ENABLE_INFO
void _unur_distr_info_typename( struct unur_gen *gen );
void _unur_distr_info_vector( struct unur_gen *gen, const double *vec, int n );
void _unur_distr_cvec_info_domain( struct unur_gen *gen );
#endif
int _unur_distr_cont_find_center( struct unur_distr *distr );
int _unur_distr_cvec_marginals_are_equal( struct unur_distr **marginals, int dim );
int _unur_distr_cvec_duplicate_firstmarginal( struct unur_distr *distribution );
int _unur_distr_cvec_is_indomain( const double *x, const struct unur_distr *distribution);
int _unur_distr_cvec_has_boundeddomain( const struct unur_distr *distribution );
#define _unur_check_distr_object( distr,distrtype, rcode ) \
  do { \
    if ((distr)->type != UNUR_DISTR_##distrtype) { \
      _unur_warning((distr)->name,UNUR_ERR_DISTR_INVALID,""); \
      return rcode; } \
    COOKIE_CHECK(distr,CK_DISTR_##distrtype,rcode); } while (0)
#endif   
