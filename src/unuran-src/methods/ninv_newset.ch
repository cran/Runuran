/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_par *
unur_ninv_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);
  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }
  par = _unur_par_new( sizeof(struct unur_ninv_par) );
  COOKIE_SET(par,CK_NINV_PAR);
  par->distr   = distr;            
  PAR->max_iter  = 100;            
  PAR->x_resolution = 1.0e-8;      
  PAR->u_resolution = -1.;         
  PAR->s[0]      = 0.0;     
  PAR->s[1]      = 0.0;     
  PAR->table_on  = FALSE;   
  par->method   = UNUR_METH_NINV;          
  par->variant  = NINV_VARFLAG_REGULA;     
  par->set      = 0u;                      
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_ninv_init;
  return par;
} 
int
unur_ninv_set_usenewton( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   
    return UNUR_ERR_DISTR_REQUIRED;
 }
  par->variant = NINV_VARFLAG_NEWTON;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_useregula( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  par->variant = NINV_VARFLAG_REGULA;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_usebisect( struct unur_par *par )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  par->variant = NINV_VARFLAG_BISECT;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_max_iter( struct unur_par *par, int max_iter )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximal iterations");
    return UNUR_ERR_PAR_SET;
  }
  PAR->max_iter = max_iter;
  par->set |= NINV_SET_MAX_ITER;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_max_iter( struct unur_gen *gen, int max_iter )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (max_iter < 1) {
    _unur_warning(gen->genid, UNUR_ERR_PAR_SET, "maximal iterations");
    return UNUR_ERR_PAR_SET;
  }
  GEN->max_iter = max_iter;
  gen->set |= NINV_SET_MAX_ITER;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_x_resolution( struct unur_par *par, double x_resolution )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (x_resolution > 0. && x_resolution < 2.*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"x-resolution too small");
    x_resolution = 2.*DBL_EPSILON;
  }
  PAR->x_resolution = x_resolution;
  par->set |= NINV_SET_X_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_x_resolution( struct unur_gen *gen, double x_resolution )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (x_resolution > 0. && x_resolution < DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"x-resolution too small");
    x_resolution = 2.*DBL_EPSILON;
  }
  GEN->x_resolution = x_resolution;
  gen->set |= NINV_SET_X_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_u_resolution( struct unur_par *par, double u_resolution )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if (u_resolution > 0. && u_resolution < 5*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too small");
    u_resolution = 1.e-15;
  }
  PAR->u_resolution = u_resolution;
  par->set |= NINV_SET_U_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_u_resolution( struct unur_gen *gen, double u_resolution )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (u_resolution > 0. && u_resolution < 5*DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"u-resolution too small");
    u_resolution = 1.e-15;
  }
  GEN->u_resolution = u_resolution;
  gen->set |= NINV_SET_U_RESOLUTION;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_start( struct unur_par *par, double s1, double s2 )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  if ( s1 <= s2 ) {
     PAR->s[0] = s1;
     PAR->s[1] = s2;
  }
  else {
     PAR->s[0] = s2;
     PAR->s[1] = s1;
  }
  par->set |= NINV_SET_START;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_start( struct unur_gen *gen, double s1, double s2 )
{
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if ( s1 <= s2 ) {
     GEN->s[0] = s1;
     GEN->s[1] = s2;
  }
  else {
     GEN->s[0] = s2;
     GEN->s[1] = s1;
  }
  GEN->table_on = FALSE;
  _unur_ninv_compute_start(gen);
  gen->set |= NINV_SET_START;
  return UNUR_SUCCESS;
} 
int
unur_ninv_set_table( struct unur_par *par, int tbl_pnts )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );
  PAR->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  PAR->table_on = TRUE;
  return UNUR_SUCCESS;
} 
int
unur_ninv_chg_table( struct unur_gen *gen, int tbl_pnts )
{
  int result;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  GEN->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  result = _unur_ninv_create_table(gen); 
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_CHG) 
    if (result==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
#endif
  return result;
} 
int 
unur_ninv_chg_truncated( struct unur_gen *gen, double left, double right )
{
  double Umin, Umax;
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }
  Umin = (left > -UNUR_INFINITY) ? CDF(left)  : 0.;
  Umax = (right < UNUR_INFINITY) ? CDF(right) : 1.;
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
  if (gen->debug & NINV_DEBUG_CHG) 
    _unur_ninv_debug_chg_truncated( gen );
#endif
  return UNUR_SUCCESS;
} 
