/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_par *
unur_tabl_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_tabl_par) );
  COOKIE_SET(par,CK_TABL_PAR);
  par->distr        = distr;      
  PAR->slopes        = NULL;      
  PAR->n_slopes      = 0;         
  PAR->n_stp         = 30;        
  PAR->cpoints       = NULL;      
  PAR->n_cpoints     = 0;         
  PAR->area_fract    = 0.1;       
  PAR->max_ivs       = 1000;      
  PAR->max_ratio     = 0.90;      
  PAR->guide_factor  = 1.; 
  PAR->darsfactor    = 0.99;   
  PAR->bleft     = -TABL_DEFAULT_COMPUTATION_LIMIT;
  PAR->bright    = TABL_DEFAULT_COMPUTATION_LIMIT;
  par->method   = UNUR_METH_TABL;              
  par->variant  = (TABL_VARFLAG_SPLIT_MEAN |   
		   TABL_VARIANT_IA         |   
		   TABL_VARFLAG_USEEAR     |   
		   TABL_VARFLAG_USEDARS    );  
  par->set      = 0u;                          
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = par->urng;               
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_tabl_init;
  return par;
} 
int
unur_tabl_set_variant_ia( struct unur_par *par, int use_ia )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  par->variant = (use_ia) ? (par->variant | TABL_VARIANT_IA) : (par->variant & (~TABL_VARIANT_IA));
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_useear( struct unur_par *par, int useear )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (useear)
    par->variant |= TABL_VARFLAG_USEEAR;
  else
    par->variant &= ~TABL_VARFLAG_USEEAR;
  par->set |= TABL_SET_USE_EAR;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_usedars( struct unur_par *par, int usedars )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (usedars)
    par->variant |= TABL_VARFLAG_USEDARS;
  else
    par->variant &= ~TABL_VARFLAG_USEDARS;
  par->set |= TABL_SET_USE_DARS;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_darsfactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->darsfactor = factor;
  par->set |= TABL_SET_DARS_FACTOR;
  return UNUR_SUCCESS;
} 
int 
unur_tabl_set_variant_splitmode( struct unur_par *par, unsigned splitmode )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  par->variant &= ~TABL_VARMASK_SPLIT;
  switch (splitmode) {
  case 1:
    par->variant |= TABL_VARFLAG_SPLIT_POINT;
    return UNUR_SUCCESS;
  case 2:
    par->variant |= TABL_VARFLAG_SPLIT_MEAN;
    return UNUR_SUCCESS;
  case 3:
    par->variant |= TABL_VARFLAG_SPLIT_ARC;
    return UNUR_SUCCESS;
  default:
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid variant");
    return UNUR_ERR_PAR_SET;
  }
} 
int
unur_tabl_set_max_sqhratio( struct unur_par *par, double max_ratio )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ratio = max_ratio;
  par->set |= TABL_SET_MAX_SQHRATIO;
  return UNUR_SUCCESS;
} 
double
unur_tabl_get_sqhratio( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, TABL, UNUR_INFINITY );
  return (GEN->Asqueeze / GEN->Atotal);
} 
double
unur_tabl_get_hatarea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, TABL, UNUR_INFINITY );
  return GEN->Atotal;
} 
double
unur_tabl_get_squeezearea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, TABL, UNUR_INFINITY );
  return GEN->Asqueeze;
} 
int
unur_tabl_set_max_intervals( struct unur_par *par, int max_ivs )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ivs = max_ivs;
  par->set |= TABL_SET_MAX_IVS;
  return UNUR_SUCCESS;
} 
int
unur_tabl_get_n_intervals( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, TABL, 0 );
  return GEN->n_ivs;
} 
int
unur_tabl_set_areafraction( struct unur_par *par, double fraction )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (fraction <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"area factor <= 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->area_fract = fraction;
  par->set |= TABL_SET_AREAFRACTION;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_cpoints( struct unur_par *par, int n_cpoints, const double *cpoints )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (n_cpoints <= 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (cpoints) 
    for( i=1; i<n_cpoints; i++ )
      if (cpoints[i] <= cpoints[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
  if (cpoints != NULL) {
    PAR->cpoints = cpoints;
    PAR->n_cpoints = n_cpoints;
  }
  else {
    PAR->n_stp = n_cpoints;
    par->set |= TABL_SET_N_STP;
  }
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_nstp( struct unur_par *par, int n_stp )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->n_stp = n_stp;
  par->set |= TABL_SET_N_STP;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_slopes( struct unur_par *par, const double *slopes, int n_slopes )
{
  int i;
  double lmax,rmin,rmax;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if( n_slopes <= 0 ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"number of slopes <= 0");
    return UNUR_ERR_PAR_SET;
  }
  lmax = -UNUR_INFINITY;
  for( i=0; i<n_slopes; i++ ) {
    rmin = _unur_min(slopes[2*i],slopes[2*i+1]);
    rmax = _unur_max(slopes[2*i],slopes[2*i+1]);
    if (!(lmax<=rmin || _unur_FP_same(lmax,rmin))) {
      _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"slopes (overlapping or not in ascending order)");
      return UNUR_ERR_PAR_SET;
    }
    lmax = rmax;
  }
  if (! (_unur_isfinite(slopes[0]) && _unur_isfinite(slopes[2*n_slopes-1])) ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"slopes must be bounded");
    return UNUR_ERR_PAR_SET;
  }
  PAR->slopes = slopes;
  PAR->n_slopes = n_slopes;
  par->set |= TABL_SET_SLOPES;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_guidefactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->guide_factor = factor;
  par->set |= TABL_SET_GUIDEFACTOR;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_boundary( struct unur_par *par, double left, double right )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  if (left >= right) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (left <= -UNUR_INFINITY || right >= UNUR_INFINITY) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- UNUR_INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }
  PAR->bleft = left;
  PAR->bright = right;
  par->set |= TABL_SET_BOUNDARY;
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  par->variant = (verify) ? (par->variant | TABL_VARFLAG_VERIFY) : (par->variant & (~TABL_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_tabl_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TABL, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  gen->variant = (verify) 
    ? (gen->variant | TABL_VARFLAG_VERIFY) 
    : (gen->variant & (~TABL_VARFLAG_VERIFY));
  SAMPLE = _unur_tabl_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
int
unur_tabl_set_pedantic( struct unur_par *par, int pedantic )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );
  par->variant = (pedantic) ? (par->variant | TABL_VARFLAG_PEDANTIC) : (par->variant & (~TABL_VARFLAG_PEDANTIC));
  return UNUR_SUCCESS;
} 
int 
unur_tabl_chg_truncated( struct unur_gen *gen, double left, double right )
{
  double Umin, Umax;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TABL, UNUR_ERR_GEN_INVALID );
  if (GEN->max_ivs > GEN->n_ivs) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"adaptive rejection sampling disabled for truncated distribution");
    GEN->max_ivs = GEN->n_ivs;
  }
  if (gen->variant & TABL_VARIANT_IA) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"cannot use IA for truncated distribution, switch to RH");
    gen->variant &= ~TABL_VARIANT_IA;
    SAMPLE = (gen->variant & TABL_VARFLAG_VERIFY) ? _unur_tabl_rh_sample_check : _unur_tabl_rh_sample;
  }
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain not subset of domain");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain not subset of domain");
    right = DISTR.domain[1];
  }
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  Umin = _unur_tabl_eval_cdfhat(gen,left);
  Umax = _unur_tabl_eval_cdfhat(gen,right);
  if (Umin > Umax) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
  if (_unur_FP_equal(Umin,Umax)) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return UNUR_ERR_DISTR_SET;
    }
  }
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = Umin;
  GEN->Umax = Umax;
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;
#ifdef UNUR_ENABLE_LOGGING
#endif
  return UNUR_SUCCESS;
} 
double
_unur_tabl_eval_cdfhat( struct unur_gen *gen, double x )
{
  struct unur_tabl_interval *iv;
  double Aint = 0.;
  double cdf;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_INFINITY);
  if (x <= DISTR.domain[0]) return 0.;
  if (x >= DISTR.domain[1]) return 1.;
  for (iv = GEN->iv; iv->next!=NULL; iv=iv->next) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_INFINITY); 
    if (x<iv->xmin || x<iv->xmax) break;
    Aint = iv->Acum; 
    if (iv->next == NULL) 
      return 1.;
  }
  Aint += iv->fmax * ((iv->xmin > iv->xmax) ? (x - iv->xmax) : (x - iv->xmin));
  cdf = Aint / GEN->Atotal;
  return ((cdf > 1.) ? 1. : cdf);
} 
