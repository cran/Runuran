/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_par *
unur_pinv_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.pdf == NULL && DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF or CDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_pinv_par) );
  COOKIE_SET(par,CK_PINV_PAR);
  par->distr   = distr;           
  PAR->order = 5L;               
  PAR->smooth = 0L;              
  PAR->u_resolution = 1.0e-10;   
  PAR->bleft = -1.e100;          
  PAR->bright = 1.e100;          
  PAR->sleft = TRUE;             
  PAR->sright = TRUE;            
  PAR->max_ivs = PINV_DEFAULT_MAX_IVS; 
  PAR->n_extra_testpoints = 0L;  
  par->method   = UNUR_METH_PINV; 
  par->variant  = 0u;             
  if (DISTR_IN.pdf != NULL)
    par->variant |= PINV_VARIANT_PDF;  
  par->set      = 0u;                      
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_pinv_init;
  return par;
} 
int
unur_pinv_set_order( struct unur_par *par, int order)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (order<3 || order>MAX_ORDER) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"order <3 or >17");
    return UNUR_ERR_PAR_SET;
  }
  PAR->order = order;
  par->set |= PINV_SET_ORDER;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_smoothness( struct unur_par *par, int smooth)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (smooth < 0L || smooth > 2L) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothness must be 0, 1, or 2");
    return UNUR_ERR_PAR_SET;
  }
  PAR->smooth = smooth;
  par->set |= PINV_SET_SMOOTH;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_u_resolution( struct unur_par *par, double u_resolution )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (u_resolution > 1.001e-5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too large --> use 1.e-5 instead");
    u_resolution = 1.e-5;
  }
  if (u_resolution < 0.999e-15 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too small --> use 1.e-15 instead");
    u_resolution = 1.e-15;
  }
  PAR->u_resolution = u_resolution;
  par->set |= PINV_SET_U_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_extra_testpoints( struct unur_par *par, int n_points)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (n_points < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of extra test point < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->n_extra_testpoints = n_points;
  par->set |= PINV_SET_N_EXTRA_TP;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_use_upoints( struct unur_par *par, int use_upoints )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (use_upoints)
    par->variant |= PINV_VARIANT_UPOINTS;
  else
    par->variant &= ~PINV_VARIANT_UPOINTS;
  par->set |= PINV_SET_UPOINTS;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_usepdf( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (par->distr->data.cont.pdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF missing");
    return UNUR_ERR_PAR_SET;
  }
  par->variant |= PINV_VARIANT_PDF;
  par->set |= PINV_SET_VARIANT;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_usecdf( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (par->distr->data.cont.cdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF missing");
    return UNUR_ERR_PAR_SET;
  }
  par->variant &= ~PINV_VARIANT_PDF;
  par->set |= PINV_SET_VARIANT;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_boundary( struct unur_par *par, double left, double right )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (!_unur_FP_less(left,right)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (! (_unur_isfinite(left) && _unur_isfinite(right)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- UNUR_INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }
  PAR->bleft = left;
  PAR->bright = right;
  par->set |= PINV_SET_BOUNDARY;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_searchboundary( struct unur_par *par, int left, int right )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  PAR->sleft  = (left)  ? TRUE : FALSE;
  PAR->sright = (right) ? TRUE : FALSE;
  par->set |= PINV_SET_SEARCHBOUNDARY;
  return UNUR_SUCCESS;
} 
int
unur_pinv_set_max_intervals( struct unur_par *par, int max_ivs )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (max_ivs < 100 || max_ivs > 1000000) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 100 or > 1000000");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_ivs = max_ivs;
  par->set |= PINV_SET_MAX_IVS;
  return UNUR_SUCCESS;
} 
int
unur_pinv_get_n_intervals( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, PINV, 0 );
  return GEN->n_ivs;
} 
int
unur_pinv_set_keepcdf( struct unur_par *par, int keepcdf)
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );
  if (keepcdf)
    par->variant |= PINV_VARIANT_KEEPCDF;
  else
    par->variant &= ~PINV_VARIANT_KEEPCDF;
  par->set |= PINV_SET_KEEPCDF;
  return UNUR_SUCCESS;
} 
