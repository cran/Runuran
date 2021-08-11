/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr.h"
#include "cxtrans.h"
#include "cont.h"
#include "distr_source.h"
static const char distr_name[] = "transformed RV";
#define DISTR distr->data.cont    
#define CXT   cxt->data.cont      
#define BD_LEFT     domain[0]     
#define BD_RIGHT    domain[1]     
#define ALPHA       params[0]     
#define MU          params[1]     
#define SIGMA       params[2]     
#define logPDFPOLE  params[3]     
#define dlogPDFPOLE params[4]     
static int _unur_distr_cxtrans_compute_domain( struct unur_distr *cxt );
static double _unur_cdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_pdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_logpdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_dpdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_dlogpdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_pdf_at_pole( const struct unur_distr *cxt );
static double _unur_dpdf_at_pole( const struct unur_distr *cxt );
struct unur_distr *
unur_distr_cxtrans_new( const struct unur_distr *distr )
{
  struct unur_distr *cxt;
  _unur_check_NULL( distr_name,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  cxt = unur_distr_cont_new();
  if (!cxt) return NULL;
  cxt->id = UNUR_DISTR_CXTRANS;
  cxt->name = distr_name;
  cxt->base = _unur_distr_cont_clone( distr );
  if (!cxt->base) { free(cxt); return NULL; }
  CXT.n_params = 5;                 
  CXT.ALPHA = 1.;                   
  CXT.MU = 0.;                      
  CXT.SIGMA = 1.;                   
  CXT.logPDFPOLE = -UNUR_INFINITY;  
  CXT.dlogPDFPOLE = UNUR_INFINITY;  
  CXT.area = DISTR.area;            
  CXT.BD_LEFT = DISTR.BD_LEFT;      
  CXT.BD_RIGHT = DISTR.BD_RIGHT;    
  CXT.mode = DISTR.mode;            
  if (DISTR.cdf)     CXT.cdf = _unur_cdf_cxtrans;          
  if (DISTR.pdf)     CXT.pdf = _unur_pdf_cxtrans;          
  if (DISTR.logpdf)  CXT.logpdf = _unur_logpdf_cxtrans;    
  if (DISTR.dpdf)    CXT.dpdf = _unur_dpdf_cxtrans;        
  if (DISTR.dlogpdf) CXT.dlogpdf = _unur_dlogpdf_cxtrans;  
  cxt->set = distr->set;
  return cxt;
} 
const struct unur_distr *
unur_distr_cxtrans_get_distribution( const struct unur_distr *cxt )
{
  _unur_check_NULL( distr_name, cxt, NULL );
  _unur_check_distr_object( cxt, CONT, NULL );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL;
  }
  return cxt->base;
} 
int
unur_distr_cxtrans_set_alpha( struct unur_distr *cxt, double alpha )
{
  double alpha_bak;
  _unur_check_NULL( distr_name, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  CHECK_NULL( cxt->base, UNUR_ERR_NULL );
  if (alpha < 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"alpha < 0");
    return UNUR_ERR_DISTR_SET;
  }
  if (_unur_iszero(alpha) && cxt->base->data.cont.BD_LEFT < 0. ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"invalid domain");
    return UNUR_ERR_DISTR_SET;
  }
  alpha_bak = CXT.ALPHA;  
  CXT.ALPHA = alpha;
  if (_unur_distr_cxtrans_compute_domain(cxt) != UNUR_SUCCESS) {
    CXT.ALPHA = alpha_bak;
    return UNUR_ERR_DISTR_SET;
  }
  cxt->set &= ~UNUR_DISTR_SET_MODE; 
  return UNUR_SUCCESS;
} 
double
unur_distr_cxtrans_get_alpha( const struct unur_distr *cxt )
{
  _unur_check_NULL( distr_name, cxt, -UNUR_INFINITY );
  _unur_check_distr_object( cxt, CONT, -UNUR_INFINITY );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return -UNUR_INFINITY; }
  return CXT.ALPHA;
} 
int
unur_distr_cxtrans_set_rescale( struct unur_distr *cxt, double mu, double sigma )
{
  double mu_bak;
  double sigma_bak;
  _unur_check_NULL( distr_name, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  CHECK_NULL( cxt->base, UNUR_ERR_NULL );
  if (sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"sigma <= 0");
    return UNUR_ERR_DISTR_SET;
  }
  mu_bak = CXT.MU;
  CXT.MU = mu;
  sigma_bak = CXT.SIGMA;
  CXT.SIGMA = sigma;
  if (_unur_distr_cxtrans_compute_domain(cxt) != UNUR_SUCCESS) {
    CXT.MU = mu_bak;
    CXT.SIGMA = sigma_bak;
    return UNUR_ERR_DISTR_SET;
  }
  cxt->set &= ~UNUR_DISTR_SET_MODE; 
  return UNUR_SUCCESS;
} 
double
unur_distr_cxtrans_get_mu( const struct unur_distr *cxt )
{
  _unur_check_NULL( distr_name, cxt, -UNUR_INFINITY );
  _unur_check_distr_object( cxt, CONT, -UNUR_INFINITY );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return -UNUR_INFINITY; }
  return CXT.MU;
} 
double
unur_distr_cxtrans_get_sigma( const struct unur_distr *cxt )
{
  _unur_check_NULL( distr_name, cxt, -UNUR_INFINITY );
  _unur_check_distr_object( cxt, CONT, -UNUR_INFINITY );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return -UNUR_INFINITY; }
  return CXT.SIGMA;
} 
int
unur_distr_cxtrans_set_logpdfpole( struct unur_distr *cxt, double logpdfpole, double dlogpdfpole )
{
  _unur_check_NULL( distr_name, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  cxt->set |= UNUR_DISTR_SET_GENERIC; 
  CXT.logPDFPOLE = logpdfpole;
  CXT.dlogPDFPOLE = dlogpdfpole;
  return UNUR_SUCCESS;
} 
int
unur_distr_cxtrans_set_domain( struct unur_distr *cxt, double left, double right )
{
  _unur_check_NULL( NULL, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  if (_unur_isinf(CXT.ALPHA)==1) {
    if (left < CXT.MU) {
      _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left < 0");
      return UNUR_ERR_DISTR_SET;
    }
  }
  return unur_distr_cont_set_domain( cxt, left, right );
} 
int _unur_distr_cxtrans_compute_domain( struct unur_distr *cxt )
{
  double left, right;
  double left_new, right_new;
  double alpha;
  CHECK_NULL( cxt, UNUR_ERR_NULL );
  CHECK_NULL( cxt->base, UNUR_ERR_NULL );
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  left  = cxt->base->data.cont.BD_LEFT;
  right = cxt->base->data.cont.BD_RIGHT;
  alpha = CXT.ALPHA;
  if (_unur_isinf(alpha)==1) {
    left_new  = (_unur_isfinite(left)) ? exp(left) : 0.;
    right_new = exp(right);
  }
  else if (_unur_iszero(alpha)) {
    if (left < 0. ) {
      _unur_error(distr_name,UNUR_ERR_DISTR_SET,"invalid domain");
      return UNUR_ERR_DISTR_SET;
    }
    left_new  = (left<=0.) ? -UNUR_INFINITY : log(left);
    right_new = log(right);
  }
  else if (alpha > 0.) {
    left_new  = (left>=0.)  ? pow(left,alpha)  : -pow(-left,alpha);
    right_new = (right>=0.) ? pow(right,alpha) : -pow(-right,alpha);
  }
  else {
    _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  if (_unur_isnan(left_new) || _unur_isnan(right_new)) {
      _unur_error(distr_name,UNUR_ERR_DISTR_SET,"NaN in now domain boundaries");
      return UNUR_ERR_DISTR_SET;
  }
  CXT.trunc[0] = CXT.domain[0] = left_new;     
  CXT.trunc[1] = CXT.domain[1] = right_new;    
  return UNUR_SUCCESS;
} 
#define CDF(x)      ((*(cxt->base->data.cont.cdf))    ((x), cxt->base))
#define PDF(x)      ((*(cxt->base->data.cont.pdf))    ((x), cxt->base))
#define logPDF(x)   ((*(cxt->base->data.cont.logpdf)) ((x), cxt->base))
#define dPDF(x)     ((*(cxt->base->data.cont.dpdf))   ((x), cxt->base))
#define dlogPDF(x)  ((*(cxt->base->data.cont.dlogpdf))((x), cxt->base))
#define POW(x)      ((x>=0.) ? pow(x,1./alpha) : -pow(-x,1./alpha))
#define dPOW(x)     ( pow(fabs(x), 1./alpha-1.) / alpha )          
#define dlogPOW(x)  ( (1./alpha-1.)*log(fabs(x)) - log(alpha) )
#define ddPOW(x)    ( ((x>=0.)?(1.-alpha):(alpha-1.)) * (_unur_isfsame(alpha,0.5)?1.0:pow(fabs(x),1./alpha-2.)) / (alpha*alpha) )
#define rescale(x)  (CXT.SIGMA * (x) + CXT.MU)
double
_unur_cdf_cxtrans( double x, const struct unur_distr *cxt )
{
  double alpha, s, mu;
  CHECK_NULL( cxt, UNUR_INFINITY );
  CHECK_NULL( cxt->base, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.cdf, UNUR_INFINITY );
  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;
  if (_unur_isinf(alpha)==1) {
    return ((x<=0.) ? 0. : CDF(s*log(x)+mu));
  }
  if (_unur_iszero(alpha)) {
    return CDF(s*exp(x)+mu);
  }
  if (alpha > 0.) {
    return CDF(s*POW(x)+mu);
  }
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return UNUR_INFINITY;
} 
double
_unur_pdf_cxtrans( double x, const struct unur_distr *cxt )
{
  double alpha, s, mu;
  CHECK_NULL( cxt, UNUR_INFINITY );
  CHECK_NULL( cxt->base, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.pdf, UNUR_INFINITY );
  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return -UNUR_INFINITY;
    else {
      double fx = PDF(s*log(x)+mu);
      return (_unur_isfinite(fx) ? fx * s/x : _unur_pdf_at_pole(cxt));
    }
  }
  if (_unur_iszero(alpha)) {
    double ex = s * exp(x) + mu;
    if (! _unur_isfinite(ex)) {
      return 0.;
    }
    else {
      double fx = PDF(ex);
      return (_unur_isfinite(fx) ? fx * s*ex :  _unur_pdf_at_pole(cxt));
    }
  }
  if (_unur_isone(alpha)) {
    double fx = PDF(s*x+mu);
    return (_unur_isfinite(fx) ? s*fx :  _unur_pdf_at_pole(cxt));
  }
  if (alpha > 0.) {
    double phix = s * POW(x) + mu;
    if (! _unur_isfinite(phix)) {
      return 0.;
    }
    else {
      double fx = PDF(phix);
      if (_unur_isfinite(fx) && (!_unur_iszero(x) || alpha < 1.)) {
	double fcx =  fx * s * dPOW(x);
	return (_unur_isfinite(fcx) ? fcx : 0.);
      }
      else 
	return  _unur_pdf_at_pole(cxt);
    }
  }
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return UNUR_INFINITY;
} 
double
_unur_logpdf_cxtrans( double x, const struct unur_distr *cxt )
{
  double alpha, s, logs, mu;
  CHECK_NULL( cxt, UNUR_INFINITY );
  CHECK_NULL( cxt->base, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.logpdf, UNUR_INFINITY );
  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;
  logs = log(CXT.SIGMA);
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return -UNUR_INFINITY;
    else {
      double logx = log(x);
      double logfx = logPDF(s*logx+mu);
      return (_unur_isfinite(logfx) ? (logfx - logx + logs) : CXT.logPDFPOLE);
    }
  }
  if (_unur_iszero(alpha)) {
    double ex = s * exp(x) + mu;
    if (! _unur_isfinite(ex)) {
      return -UNUR_INFINITY;
    }
    else {
      double logfx = logPDF(ex);
      return (_unur_isfinite(logfx) ? (logfx + x + logs) : CXT.logPDFPOLE);
    }
  }
  if (_unur_isone(alpha)) {
    double logfx = logPDF(s*x+mu);
    return (_unur_isfinite(logfx) ? (logfx + logs) : CXT.logPDFPOLE);
  }
  if (alpha > 0.) {
    double phix = s * POW(x) + mu;
    if (! _unur_isfinite(phix)) {
      return -UNUR_INFINITY;
    }
    else {
      double logfx = logPDF(phix);
      if (_unur_isfinite(logfx) && (!_unur_iszero(x) || alpha < 1.)) {
	double logfcx =  logfx + logs + dlogPOW(x);
	return (_unur_isfinite(logfcx) ? logfcx : -UNUR_INFINITY);
      }
      else 
	return CXT.logPDFPOLE;
    }
  }
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return UNUR_INFINITY;
} 
double
_unur_dpdf_cxtrans( double x, const struct unur_distr *cxt )
{
  double alpha, s, mu;
  CHECK_NULL( cxt, UNUR_INFINITY );
  CHECK_NULL( cxt->base, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.pdf, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.dpdf, UNUR_INFINITY );
  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return 0.;
    else {
      double logx = s*log(x)+mu;
      double fx = PDF(logx);
      double dfx = dPDF(logx);
      return (_unur_isfinite(fx) ? s*(s*dfx - fx)/(x*x) : _unur_dpdf_at_pole(cxt));
    }
  }
  if (_unur_iszero(alpha)) {
    double ex = s*exp(x)+mu;
    if (! _unur_isfinite(ex)) {
      return 0.;
    }
    else {
      double fx = PDF(ex);
      double dfx = dPDF(ex);
      double dfcx = s * (dfx * s*ex*ex + fx * ex);
      if (! _unur_isfinite(fx) ) return _unur_dpdf_at_pole(cxt);
      if (! _unur_isfinite(dfcx) ) return (dfx>0 ? UNUR_INFINITY : -UNUR_INFINITY);
      return dfcx;
    }
  }
  if (_unur_isone(alpha)) {
    double fx = PDF(s*x+mu);
    double dfx = dPDF(s*x+mu);
    return (_unur_isfinite(fx) ? s*dfx : _unur_dpdf_at_pole(cxt));
  }
  if (alpha > 0.) {
    double phix = s*POW(x)+mu;
    if (! _unur_isfinite(phix)) {
      return 0.;
    }
    else {
      double fx = PDF(phix);
      double dfx = dPDF(phix);
      double dphix = dPOW(x);
      double ddphix = ddPOW(x);
      if (_unur_isfinite(fx) && (!_unur_iszero(x) || alpha <= 0.5)) {
	double dfcx = s*(dfx * s*dphix*dphix + fx * s*ddphix);
	return (_unur_isfinite(dfcx) ? dfcx : 0.);
      }
      else
	return _unur_dpdf_at_pole(cxt);
    }
  }
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return UNUR_INFINITY;
} 
double
_unur_dlogpdf_cxtrans( double x, const struct unur_distr *cxt )
{
  double alpha, s, mu;
  CHECK_NULL( cxt, UNUR_INFINITY );
  CHECK_NULL( cxt->base, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.logpdf, UNUR_INFINITY );
  CHECK_NULL( cxt->base->data.cont.dlogpdf, UNUR_INFINITY );
  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return -UNUR_INFINITY;
    else {
      double logx = s*log(x)+mu;
      double logfx = logPDF(logx);
      double dlogfx = dlogPDF(logx);
      return (_unur_isfinite(logfx) ? ((s*dlogfx-1)/x) : CXT.dlogPDFPOLE);
    }
  }
  if (_unur_iszero(alpha)) {
    double ex = s*exp(x)+mu;
    if (! _unur_isfinite(ex)) {
      return (x>1. ? -UNUR_INFINITY : UNUR_INFINITY);
    }
    else {
      double logfx = logPDF(ex);
      double dlogfx = dlogPDF(ex);
      return (_unur_isfinite(logfx) ? s*dlogfx*ex + 1 : CXT.dlogPDFPOLE);
    }
  }
  if (_unur_isone(alpha)) {
    double logfx = logPDF(x);
    return (_unur_isfinite(logfx) ? s*dlogPDF(x) : CXT.dlogPDFPOLE);
  }
  if (alpha > 0.) {
    double phix = s*POW(x)+mu;
    if (! _unur_isfinite(phix)) {
      return ((x>1. || (x>-1. && x < 0.)) ? -UNUR_INFINITY : UNUR_INFINITY);
    }
    else {
      double logfx = logPDF(phix);
      if (_unur_isfinite(logfx) && (!_unur_iszero(x) || alpha <= 1.)) {
	double dlogfcx = ((x>=0.)?1.:-1.) * (dlogPDF(phix) * s*dPOW(x) + (1./alpha-1.)/x);
	if (! _unur_isfinite(dlogfcx)) {
	  return ((x>1. || (x>-1. && x < 0.)) ? -UNUR_INFINITY : UNUR_INFINITY);
	}
	return dlogfcx;
      }
      else
	return CXT.dlogPDFPOLE;
    }
  }
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return UNUR_INFINITY;
} 
double 
_unur_pdf_at_pole( const struct unur_distr *cxt )
{
  return exp(CXT.logPDFPOLE);
} 
double
_unur_dpdf_at_pole( const struct unur_distr *cxt )
{
  double fx = _unur_pdf_at_pole(cxt);
  if (! (_unur_isfinite(CXT.logPDFPOLE) && _unur_isfinite(fx)) )
    return UNUR_INFINITY;
  else 
    return (fx * CXT.dlogPDFPOLE);
} 
#ifdef UNUR_ENABLE_LOGGING
void
_unur_distr_cxtrans_debug( const struct unur_distr *cxt, const char *genid )
{
  FILE *LOG;
  CHECK_NULL(cxt,RETURN_VOID);
  COOKIE_CHECK(cxt,CK_DISTR_CONT,RETURN_VOID);
  CHECK_NULL(cxt->base,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous univariate distribution of transformed random variable\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,cxt->name);
  fprintf(LOG,"%s:\talpha = %g\t",genid,CXT.ALPHA);
  if (_unur_isinf(CXT.ALPHA)==1)
    fprintf(LOG,"[ exponential transformation: Y = exp(Z) ]\n"); 
  else if (_unur_iszero(CXT.ALPHA))
    fprintf(LOG,"[ logarithmic transformation: Y = log(Z) ]\n"); 
  else
    fprintf(LOG,"[ power transformation: Y = Z^alpha ]\n"); 
  fprintf(LOG,"%s:\tmu = %g, sigma = %g\t[Z = (X-%g)/%g]\n",genid, CXT.MU, CXT.SIGMA, CXT.MU, CXT.SIGMA);
  fprintf(LOG,"%s:\n",genid);
  fprintf(LOG,"%s:\tvalues used at pole of underlying distribution\n",genid);
  fprintf(LOG,"%s:\t\tlogPDF  = %g\t(PDF  = %g)",genid, CXT.logPDFPOLE, _unur_pdf_at_pole(cxt));
  _unur_print_if_default(cxt,UNUR_DISTR_SET_GENERIC);
  fprintf(LOG,"\n");
  fprintf(LOG,"%s:\t\tdlogPDF = %g\t(dPDF = %g)",genid, CXT.dlogPDFPOLE, _unur_dpdf_at_pole(cxt));
  _unur_print_if_default(cxt,UNUR_DISTR_SET_GENERIC);
  fprintf(LOG,"\n");
  if (cxt->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %g\n",genid,CXT.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)",genid,CXT.BD_LEFT,CXT.BD_RIGHT);
  _unur_print_if_default(cxt,UNUR_DISTR_SET_DOMAIN);
  fprintf(LOG,"\n%s:\tarea below PDF = %g",genid,CXT.area);
  _unur_print_if_default(cxt,UNUR_DISTR_SET_PDFAREA);
  fprintf(LOG,"\n%s:\n",genid);
  fprintf(LOG,"%s: Underlying distribution:\n",genid);
  _unur_distr_cont_debug(cxt->base, genid);
} 
#endif    
