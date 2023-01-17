/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_par *
unur_tdr_new( const struct unur_distr* distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL; }
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of PDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_tdr_par) );
  COOKIE_SET(par,CK_TDR_PAR);
  par->distr              = distr;  
  PAR->guide_factor        = 2.;    
  PAR->c_T                 = -0.5;  
  PAR->starting_cpoints    = NULL;  
  PAR->n_starting_cpoints  = 30;    
  PAR->percentiles         = NULL;  
  PAR->n_percentiles       = 2;     
  PAR->retry_ncpoints      = 50;    
  PAR->max_ivs             = 100;   
  PAR->max_ratio           = 0.99;  
  PAR->bound_for_adding    = 0.5;   
  PAR->darsfactor          = 0.99;   
  PAR->darsrule            = 1;     
  par->method   = UNUR_METH_TDR;                 
  par->variant  = ( TDR_VARFLAG_USECENTER |      
		    TDR_VARFLAG_USEMODE   |
                    TDR_VARIANT_PS );
  par->set      = 0u;                   
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = par->urng;               
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_tdr_init;
  return par;
} 
int
unur_tdr_set_cpoints( struct unur_par *par, int n_stp, const double *stp )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }
  if (stp) 
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
  PAR->starting_cpoints = stp;
  PAR->n_starting_cpoints = n_stp;
  par->set |= TDR_SET_N_STP | ((stp) ? TDR_SET_STP : 0);
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_reinit_percentiles( struct unur_par *par, int n_percentiles, const double *percentiles )
{
  int i;
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }
  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
  if (percentiles) {
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }
  PAR->percentiles = percentiles;
  PAR->n_percentiles = n_percentiles;
  par->set |= TDR_SET_N_PERCENTILES | ((percentiles) ? TDR_SET_PERCENTILES : 0);
  return UNUR_SUCCESS;
} 
int
unur_tdr_chg_reinit_percentiles( struct unur_gen *gen, int n_percentiles, const double *percentiles )
{
  int i;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }
  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
  if (percentiles) {
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }
  GEN->n_percentiles = n_percentiles;
  GEN->percentiles = _unur_xrealloc( GEN->percentiles, n_percentiles * sizeof(double) );
  if (percentiles) {
    memcpy( GEN->percentiles, percentiles, n_percentiles * sizeof(double) );
  }
  else {
    if (n_percentiles == 2) {
      GEN->percentiles[0] = 0.25;
      GEN->percentiles[1] = 0.75;
    }
    else {
      for (i=0; i<n_percentiles; i++ )
	GEN->percentiles[i] = (i + 1.) / (n_percentiles + 1.);
    }
  }
  gen->set |= TDR_SET_N_PERCENTILES | ((percentiles) ? TDR_SET_PERCENTILES : 0);
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_reinit_ncpoints( struct unur_par *par, int ncpoints )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }
  PAR->retry_ncpoints = ncpoints;
  par->set |= TDR_SET_RETRY_NCPOINTS; 
  return UNUR_SUCCESS;
} 
int
unur_tdr_chg_reinit_ncpoints( struct unur_gen *gen, int ncpoints )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }
  GEN->retry_ncpoints = ncpoints;
  gen->set |= TDR_SET_RETRY_NCPOINTS; 
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_guidefactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->guide_factor = factor;
  par->set |= TDR_SET_GUIDEFACTOR;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_max_sqhratio( struct unur_par *par, double max_ratio )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (max_ratio < 0. || max_ratio > 1.+DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ratio = max_ratio;
  par->set |= TDR_SET_MAX_SQHRATIO;
  return UNUR_SUCCESS;
} 
double
unur_tdr_get_sqhratio( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, TDR, UNUR_INFINITY );
  return (GEN->Asqueeze / GEN->Atotal);
} 
double
unur_tdr_get_hatarea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, TDR, UNUR_INFINITY );
  return GEN->Atotal;
} 
double
unur_tdr_get_squeezearea( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, TDR, UNUR_INFINITY );
  return GEN->Asqueeze;
} 
int
unur_tdr_set_max_intervals( struct unur_par *par, int max_ivs )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ivs = max_ivs;
  par->set |= TDR_SET_MAX_IVS;
  return UNUR_SUCCESS;
} 
int
_unur_tdr_is_ARS_running( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, FALSE );
  _unur_check_gen_object( gen, TDR, FALSE );
  return (GEN->n_ivs < GEN->max_ivs) ? TRUE : FALSE;
} 
int
unur_tdr_set_usecenter( struct unur_par *par, int usecenter )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (usecenter) ? (par->variant | TDR_VARFLAG_USECENTER) : (par->variant & (~TDR_VARFLAG_USECENTER));
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_usemode( struct unur_par *par, int usemode )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (usemode) ? (par->variant | TDR_VARFLAG_USEMODE) : (par->variant & (~TDR_VARFLAG_USEMODE));
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_variant_gw( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_GW;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_variant_ps( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_PS;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_variant_ia( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_IA;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_usedars( struct unur_par *par, int usedars )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (usedars < 0 || usedars > 3) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid rule for DARS");
    return UNUR_ERR_PAR_SET;
  }
  PAR->darsrule = usedars;
  par->variant = (usedars) ? (par->variant | TDR_VARFLAG_USEDARS) : (par->variant & (~TDR_VARFLAG_USEDARS));
  par->set |= TDR_SET_USE_DARS;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_darsfactor( struct unur_par *par, double factor )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->darsfactor = factor;
  par->set |= TDR_SET_DARS_FACTOR;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_c( struct unur_par *par, double c )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  if (c > 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c > 0");
    return UNUR_ERR_PAR_SET;
  }
  if (c < -0.5) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"c < -0.5 not implemented yet");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_iszero(c) && c > -0.5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"-0.5 < c < 0 not recommended. using c = -0.5 instead.");
    c = -0.5;
  }
  PAR->c_T = c;
  par->set |= TDR_SET_C;
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (verify) ? (par->variant | TDR_VARFLAG_VERIFY) : (par->variant & (~TDR_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_tdr_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;
  gen->variant = (verify) 
    ? (gen->variant | TDR_VARFLAG_VERIFY) 
    : (gen->variant & (~TDR_VARFLAG_VERIFY));
  SAMPLE = _unur_tdr_getSAMPLE(gen);
  return UNUR_SUCCESS;
} 
int
unur_tdr_set_pedantic( struct unur_par *par, int pedantic )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDR );
  par->variant = (pedantic) ? (par->variant | TDR_VARFLAG_PEDANTIC) : (par->variant & (~TDR_VARFLAG_PEDANTIC));
  return UNUR_SUCCESS;
} 
int 
unur_tdr_chg_truncated( struct unur_gen *gen, double left, double right )
{
  double Umin, Umax;
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDR, UNUR_ERR_GEN_INVALID );
  if (GEN->max_ivs > GEN->n_ivs) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"adaptive rejection sampling disabled for truncated distribution");
    GEN->max_ivs = GEN->n_ivs;
  }
  if ((gen->variant & TDR_VARMASK_VARIANT) == TDR_VARIANT_IA) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"cannot use IA for truncated distribution, switch to PS");
    gen->variant = (gen->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_PS;
    SAMPLE = (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ps_sample_check : _unur_tdr_ps_sample;
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
  Umin = _unur_tdr_eval_cdfhat(gen,left);
  Umax = (right < DISTR.domain[1]) ? _unur_tdr_eval_cdfhat(gen,right) : 1.;
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
_unur_tdr_eval_cdfhat( struct unur_gen *gen, double x )
{
  struct unur_tdr_interval *iv;
  double Aint;
  double cdf;
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_INFINITY);
  if (x <= DISTR.domain[0]) return 0.;
  if (x >= DISTR.domain[1]) return 1.;
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    
    for (iv = GEN->iv; iv->next!=NULL; iv=iv->next) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_INFINITY); 
      if (x < iv->next->x) break;
    }
    if (iv->next == NULL)
      return 1.;
    if (x < iv->ip) {
      Aint = _unur_tdr_interval_area( gen, iv, iv->dTfx, x);
      if (!_unur_isfinite(Aint)) { 
	Aint = 0.;
      }
      cdf = (iv->prev) ? iv->prev->Acum + Aint : Aint;
    }
    else {
      Aint = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, x);
      if (!_unur_isfinite(Aint)) { 
	Aint = 0.;
      }
      cdf = iv->Acum - Aint;
      if (cdf < 0.) return 0.;
    }
    cdf /= GEN->Atotal;
    return ((cdf > 1.) ? 1. : cdf);
  case TDR_VARIANT_IA:    
  case TDR_VARIANT_PS:    
    for (iv = GEN->iv; iv->next!=NULL; iv=iv->next) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_INFINITY); 
      if (x <= iv->next->ip) break;
    }
    if (iv->next == NULL)
      return 1.;
    Aint = _unur_tdr_interval_area( gen, iv, iv->dTfx, x);
    if (!_unur_isfinite(Aint)) { 
      Aint = 0.;
    }
    cdf = ((x>iv->x) ? Aint : -Aint) + iv->Acum - iv->Ahatr;
    if (cdf < 0.) return 0.;
    cdf /= GEN->Atotal;
    return ((cdf > 1.) ? 1. : cdf);
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_INFINITY;
  }
} 
