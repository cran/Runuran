/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

int
_unur_mvtdr_sample_cvec( struct unur_gen *gen, double *rpoint )
{
  CONE *c;       
  double gx;     
  double U;      
  double f, h;   
  int i,j;
  double *S = GEN->S;  
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MVTDR_GEN,UNUR_ERR_COOKIE);
  while( 1 ) { 
    U = _unur_call_urng(gen->urng);      
    c = (GEN->guide)[(int) (U * GEN->guide_size)]; 
    U *= GEN->Htot;
    while (c->next!=NULL && c->Hsum < U) 
      c = c->next;
    if (GEN->has_domain)
      unur_tdr_chg_truncated(GEN_GAMMA, 0., c->beta * c->height );
    gx = unur_sample_cont(GEN_GAMMA) / (c->beta);
    _unur_mvtdr_simplex_sample(gen, S);
    for( i=0; i<GEN->dim; i++ ) rpoint[i] = GEN->center[i];
    for( j=0; j<GEN->dim; j++ ) {
      double x = gx * S[j] / c->gv[j];
      for( i=0; i<GEN->dim; i++ )
	rpoint[i] += x * (c->v[j])->coord[i];
    }
    f = PDF(rpoint);                        
    h = T_inv( c->alpha - c->beta * gx );   
    if ( (gen->variant & MVTDR_VARFLAG_VERIFY) &&
	 ((1.+UNUR_EPSILON) * h < f ) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    if( _unur_call_urng(gen->urng) * h <= f )
      return UNUR_SUCCESS;
  }
} 
int
_unur_mvtdr_simplex_sample( const struct unur_gen *gen, double *U )
{
  int dim = GEN->dim;
  if (dim == 2) {
    U[0] = _unur_call_urng(gen->urng);
    U[1] = 1. - U[0];
    return UNUR_SUCCESS;
  }
  if (dim == 3) {
    U[0] = _unur_call_urng(gen->urng);
    U[1] = _unur_call_urng(gen->urng);
    if( U[0] > U[1] ) {
      U[2] = U[0]; U[0] = U[1]; U[1] = U[2];
    }
    U[2] = 1. - U[1];
    U[1] = U[1] - U[0];
    return UNUR_SUCCESS;
  }
  if (dim >3) {
    int i,j;
    double U_aux;
    for( i=0; i<dim-1; i++ )
      U[i] = _unur_call_urng(gen->urng);
    for( i=1; i<dim-1; i++ ) {
      U_aux = U[i];
      for( j=i; j>0 && U[j-1] > U_aux; j-- )
	U[j] = U[j-1];
      U[j] = U_aux;
    }
    U[dim-1] = 1.;
    for( i=dim-1; i>0; i-- )
      U[i] -= U[i-1];
    return UNUR_SUCCESS;
  }
  _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return UNUR_FAILURE;
} 
