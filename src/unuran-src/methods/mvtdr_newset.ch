/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_par *
unur_mvtdr_new( const struct unur_distr *distr )
{ 
  struct unur_par *par;
  _unur_check_NULL( GENTYPE,distr,NULL );
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);
  if (distr->dim < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_PROP,"dim < 2"); return NULL; }
  if ( ! ((DISTR_IN.pdf && DISTR_IN.dpdf) || (DISTR_IN.logpdf && DISTR_IN.dlogpdf)) ) { 
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"d/(log)PDF");
    return NULL;
  }
  par = _unur_par_new( sizeof(struct unur_mvtdr_par) );
  COOKIE_SET(par,CK_MVTDR_PAR);
  par->distr    = distr;      
  par->method   = UNUR_METH_MVTDR ;   
  par->variant  = 0u;                 
  par->set      = 0u;                     
  par->urng     = unur_get_default_urng(); 
  par->urng_aux = NULL;                    
  par->debug    = _unur_default_debugflag; 
  par->init = _unur_mvtdr_init;
  PAR->steps_min = 5;
  PAR->max_cones = 10000;
  PAR->bound_splitting = 1.5;
  return par;
} 
int 
unur_mvtdr_set_stepsmin( struct unur_par *par, int stepsmin )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );
  if (stepsmin < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"stepsmin < 0");
    return UNUR_ERR_PAR_SET;
  }
  PAR->steps_min = stepsmin;
  par->set |= MVTDR_SET_STEPSMIN;
  return UNUR_SUCCESS;
} 
int
unur_mvtdr_set_boundsplitting( UNUR_PAR *par, double boundsplitting )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );
  PAR->bound_splitting = boundsplitting;
  par->set |= MVTDR_SET_BOUNDSPLITTING;
  return UNUR_SUCCESS;
} 
int 
unur_mvtdr_set_maxcones( struct unur_par *par, int maxcones )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );
  PAR->max_cones = maxcones;
  par->set |= MVTDR_SET_MAXCONES;
  return UNUR_SUCCESS;
} 
int
unur_mvtdr_get_ncones( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, MVTDR, 0 );
  return GEN->n_cone;
} 
double
unur_mvtdr_get_hatvol( const struct unur_gen *gen )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  _unur_check_gen_object( gen, MVTDR, UNUR_INFINITY );
  return GEN->Htot;
} 
int
unur_mvtdr_set_verify( struct unur_par *par, int verify )
{
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MVTDR );
  par->variant = (verify) ? (par->variant | MVTDR_VARFLAG_VERIFY) : (par->variant & (~MVTDR_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
int
unur_mvtdr_chg_verify( struct unur_gen *gen, int verify )
{
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, MVTDR, UNUR_ERR_GEN_INVALID );
  if (SAMPLE == _unur_sample_cvec_error) 
    return UNUR_FAILURE;
  gen->variant = (verify) 
    ? (gen->variant | MVTDR_VARFLAG_VERIFY) 
    : (gen->variant & (~MVTDR_VARFLAG_VERIFY));
  return UNUR_SUCCESS;
} 
