/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

int
_unur_pinv_preprocessing (struct unur_gen *gen)
{
  switch (gen->variant) {
  case PINV_VARIANT_PDF:
    if (_unur_pinv_relevant_support(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    if (_unur_pinv_approx_pdfarea(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    if (_unur_pinv_computational_domain(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    if (_unur_pinv_pdfarea(gen) != UNUR_SUCCESS) 
      return UNUR_FAILURE;
    break;
  case PINV_VARIANT_CDF:
    if (_unur_pinv_computational_domain_CDF(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
  return UNUR_SUCCESS;
} 
int
_unur_pinv_relevant_support ( struct unur_gen *gen )
{
  double fb;
  if(GEN->sleft) {
    fb = PDF(GEN->dleft);
    if (fb > 1.e-20 && fb < 1.e300) {
      GEN->bleft = GEN->dleft;        
      GEN->sleft = FALSE;
    }
  }
  if(GEN->sright) {
    fb = PDF(GEN->dright);
    if (fb > 1.e-20 && fb < 1.e300) {
      GEN->bright = GEN->dright;        
      GEN->sright = FALSE;
    }
  }
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_searchborder(gen, DISTR.center, GEN->bleft, 
					 &(GEN->dleft), &(GEN->sleft) );
    if (!_unur_isfinite(GEN->bleft)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get left boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  if(GEN->sright) {
    GEN->bright = _unur_pinv_searchborder(gen, DISTR.center, GEN->bright,
					  &(GEN->dright), &(GEN->sright) );
    if (!_unur_isfinite(GEN->bright)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get right boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_relevant_support(gen);
#endif
  return UNUR_SUCCESS;
} 
double
_unur_pinv_searchborder (struct unur_gen *gen, double x0, double bound,
			 double *dom, int *search)
{
  double x;         
  double xs, xl;    
  double fx;        
  double fs, fl;    
  double fllim;     
  double fulim;     
  fllim = PDF(x0) * PINV_PDFLLIM;
  fulim = 1.e4 * fllim;
  if (fllim <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(center) too small");
    return INFINITY;
  }
  xl = x0; 
  fl = INFINITY;
  x = _unur_arcmean(x0,bound);
  while ( (fx=PDF(x)) > fllim ) {
    if (_unur_FP_same(x,bound))
      return bound;
    xl = x; fl = fx;
    x = _unur_arcmean(x,bound);
  }
  xs = x; fs = fx;
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0");
    return INFINITY;
  }
  while (!_unur_FP_same(xs,xl)) {
    if (_unur_iszero(fs)) {
      *dom = xs;
    }
    x = xs/2. + xl/2.;
    fx = PDF(x);
    if (fx < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0");
      return INFINITY;
    }
    if (fx < fllim) {
      xs = x; fs = fx;
    }
    else {
      if (fl > fulim) {
	xl = x; fl = fx;
      }
      else {
	return x;      
      }
    }
  }
  *search = FALSE;
  return xl;
} 
int
_unur_pinv_approx_pdfarea (struct unur_gen *gen )
{
  double tol;   
  int i;        
  int res = UNUR_SUCCESS; 
  for (i=1; i<=2; i++) {
    tol = PINV_UERROR_AREA_APPROX * GEN->area;
    DISTR.center = _unur_max(DISTR.center, GEN->bleft);
    DISTR.center = _unur_min(DISTR.center, GEN->bright);
    GEN->area  = 
      _unur_lobatto_adaptive(_unur_pinv_eval_PDF, gen,
			     GEN->bleft, DISTR.center - GEN->bleft, tol, NULL)
      + _unur_lobatto_adaptive(_unur_pinv_eval_PDF, gen,
			       DISTR.center, GEN->bright - DISTR.center, tol, NULL);
    if ( !_unur_isfinite(GEN->area) || _unur_iszero(GEN->area) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot estimate area below PDF");
      res = UNUR_FAILURE;
      break;
    }
    if (GEN->area > 1.e-2) {
      break;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_pdfarea(gen,TRUE);
#endif
  return res;
} 
int
_unur_pinv_pdfarea (struct unur_gen *gen)
{
  double tol;   
  tol = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION * PINV_UTOL_CORRECTION;
  DISTR.center = _unur_max(DISTR.center, GEN->bleft);
  DISTR.center = _unur_min(DISTR.center, GEN->bright);
  GEN->aCDF = _unur_lobatto_init(_unur_pinv_eval_PDF, gen,
				 GEN->bleft, DISTR.center, GEN->bright,
				 tol, NULL, PINV_MAX_LOBATTO_IVS);
  GEN->area = _unur_lobatto_integral(GEN->aCDF);
  if ( !_unur_isfinite(GEN->area) || _unur_iszero(GEN->area) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot estimate area below PDF");
    return UNUR_FAILURE;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_pdfarea(gen,FALSE);
#endif
  return UNUR_SUCCESS;
} 
int
_unur_pinv_computational_domain (struct unur_gen *gen)
{
  double tailcut_error;    
  double range;            
  tailcut_error = GEN->u_resolution * PINV_TAILCUTOFF_FACTOR;
  tailcut_error = _unur_min( tailcut_error, PINV_TAILCUTOFF_MAX );
  tailcut_error = _unur_max( tailcut_error, 2*DBL_EPSILON );
  tailcut_error *= GEN->area * PINV_UERROR_CORRECTION;
  range = GEN->bright-GEN->bleft;
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_cut( gen, GEN->dleft, GEN->bleft, -range, tailcut_error);
    if ( !_unur_isfinite(GEN->bleft) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find left boundary for computational domain");
      return UNUR_FAILURE;
    }
  }
  if(GEN->sright) {
    GEN->bright = _unur_pinv_cut( gen, GEN->dright, GEN->bright, range, tailcut_error);
    if ( !_unur_isfinite(GEN->bright) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find right boundary for computational domain");
      return UNUR_FAILURE;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_computational_domain(gen);
#endif
  return UNUR_SUCCESS;
} 
double
_unur_pinv_cut( struct unur_gen *gen, double dom, double w, double dw, double crit )
{
  double sgn = (dw>0) ? 1. : -1.; 
  double fl,fx,fr;  
  double x = w;     
  double dx;                                           
  double xnew;      
  double df;        
  double lc;        
  double area;      
  int i,j;          
  if (_unur_iszero(fabs(dw))) return w;
  x = w;
  for (i=1; i<100; i++) {
    dx = (fabs(dw) + fabs(x-w)) * 1.e-3;
    if (x-dx < GEN->dleft)  dx = x - GEN->dleft;
    if (x+dx > GEN->dright) dx = GEN->dright - x;
    for (j=1;;j++) {
      dx = dx/2.;
      if (dx < 128.*DBL_EPSILON*fabs(dw)) {
	return x;
      }
      fx = PDF(x);
      fl = PDF(x-dx);
      fr = PDF(x+dx);
      if (! (_unur_iszero(fl) || _unur_iszero(fx) ||_unur_iszero(fr)) )
	break;
    }
    df = (fr-fl)/(2.*dx);
    if (! _unur_isfinite(df)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "numerical problems with cut-off point, PDF too steep");
      return INFINITY;
    }
    lc = fl/(fl-fx)+fr/(fr-fx) - 1;
    area = fabs(fx*fx / ((lc+1.) * df));
    if (_unur_isnan(area)) {
      _unur_warning(gen->genid,UNUR_ERR_NAN,"tail probability gives NaN --> assume 0.");
      return x;
    }
    if (sgn * df > 0.) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not monotone at boundary");
      return x;
    }
    if (fabs(area/crit-1.) < 1.e-4)
      return x;
    if (_unur_iszero(lc)) {
      xnew = x + fx/df * log(crit*fabs(df)/(fx*fx));
    }
    else {
      xnew = x + fx/(lc*df) * ( pow(crit*fabs(df)*(lc+1.)/(fx*fx),lc/(lc+1.)) - 1.);
    }
    if (! _unur_isfinite(xnew)) {
      _unur_error(gen->genid,UNUR_ERR_NAN,"numerical problems with cut-off point");
      return INFINITY;
    }
    if (sgn*dom < sgn*x) {
      return dom;
    }
    x = xnew;
  }
  return x;
} 
int
_unur_pinv_computational_domain_CDF (struct unur_gen *gen)
{
  double tailcut_error;    
  double fl, fr;
  fl = CDF(DISTR.domain[0]);
  fr = CDF(DISTR.domain[1]);
  if (_unur_FP_approx(fl,fr)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"truncated domain too narrow");
    return UNUR_FAILURE;
  }
  tailcut_error = GEN->u_resolution * PINV_TAILCUTOFF_FACTOR;
  tailcut_error = _unur_min( tailcut_error, PINV_TAILCUTOFF_MAX );
  tailcut_error = _unur_max( tailcut_error, 2*DBL_EPSILON );
  tailcut_error *= GEN->area * PINV_UERROR_CORRECTION;
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_cut_CDF( gen, GEN->dleft, DISTR.center, 0.5*tailcut_error, tailcut_error);
    if ( !_unur_isfinite(GEN->bleft) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find left boundary for computational domain");
      return UNUR_FAILURE;
    }
  }
  if(GEN->sright) {
    GEN->bright = _unur_pinv_cut_CDF( gen, GEN->dright, DISTR.center, 1.-tailcut_error, 1.-0.5*tailcut_error);
    if ( !_unur_isfinite(GEN->bright) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find right boundary for computational domain");
      return UNUR_FAILURE;
    }
  }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_computational_domain(gen);
#endif
  return UNUR_SUCCESS;
} 
double
_unur_pinv_cut_CDF( struct unur_gen *gen, double dom, double x0, double ul, double uu )
{
  double x;         
  double xs, xl;    
  double fx;        
  double fs, fl;    
  double f0, fdom;  
  double dx;        
  if (_unur_FP_same(x0,dom))
      return x0;
  if (1.-ul < 4*DBL_EPSILON) ul = 1. - 4*DBL_EPSILON;
  if (1.-uu < 2*DBL_EPSILON) ul = 1. - 2*DBL_EPSILON;
  x = x0;
  f0 = CDF(x0);
  fdom = CDF(dom);
  if (_unur_iszero(f0)) {
    for (dx=0.1; f0<ul; dx*=10.) {
      dom = x0; fdom = f0;
      x0 += dx;
      f0 = CDF(x0);
      if (!_unur_isfinite(x0))
	return INFINITY;
    }
  }
  if (_unur_isone(f0)) {
    for (dx=0.1; f0>ul; dx*=10.) {
      dom = x0; fdom = f0;
      x0 -= dx;
      f0 = CDF(x0);
      if (!_unur_isfinite(x0))
	return INFINITY;
    }
  }
  if ( (f0 < ul && fdom < ul) || (f0 > uu && fdom > uu) ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"CDF too small/large on given domain");
    return dom;
  }
  if (f0 >= ul && f0 <= uu) {
    return x0;
  }
  if ( (x0 < dom && _unur_FP_greater(f0,fdom)) ||
       (x0 > dom && _unur_FP_less(f0,fdom)) ) {
    return INFINITY;
  }
  if (x0 > dom) {
    xs = dom; fs = fdom;
    xl = x0; fl = f0; }
  else {
    xs = x0; fs = f0;
    xl = dom; fl = fdom;
  }    
  x = x0;
  while (!_unur_FP_same(xs,xl)) {
    x = _unur_arcmean(xs,xl);
    fx = CDF(x);
    if (fx >= ul && fx <= uu) {
      return x;
    }
    if (fx < ul) {
      xs = x; fs = fx;
    }
    else {
      xl = x; fl = fx;
    }
  }
  return x;
} 
double
_unur_pinv_Udiff (struct unur_gen *gen, double x, double h, double *fx)
{
  switch (gen->variant) {
  case PINV_VARIANT_PDF:
    return _unur_lobatto_eval_diff(GEN->aCDF, x, h, fx);
  case PINV_VARIANT_CDF:
    return CDF(x+h) - CDF(x);
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }
} 
