/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_gen *
_unur_mvtdr_init( struct unur_par *par )
{ 
  struct unur_gen *gen;
  _unur_check_NULL( GENTYPE,par,NULL );
  if ( par->method != UNUR_METH_MVTDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MVTDR_PAR,NULL);
  gen = _unur_mvtdr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_mvtdr_debug_init_start(gen);
#endif
  if (!(GEN->pdfcenter > 0.)) {
    _unur_error(gen->genid,UNUR_ERR_DISTR_DOMAIN,"center out of support of PDF");
    _unur_mvtdr_free(gen); return NULL;
  }
  if(_unur_mvtdr_create_hat(gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create hat");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_mvtdr_debug_init_finished(gen, FALSE);
#endif
    _unur_mvtdr_free(gen); return NULL;
  }
  _unur_mvtdr_max_gamma(gen);
  GEN_GAMMA = _unur_mvtdr_gammagen( gen, (double)(GEN->dim) );
  if ( GEN_GAMMA == NULL ) {
      _unur_mvtdr_free(gen); return NULL; }
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_mvtdr_debug_init_finished(gen, TRUE);
#endif
  return gen;
} 
static struct unur_gen *
_unur_mvtdr_create( struct unur_par *par )
{
  struct unur_gen *gen;
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MVTDR_PAR,NULL);
  gen = _unur_generic_create( par, sizeof(struct unur_mvtdr_gen) );
  COOKIE_SET(gen,CK_MVTDR_GEN);
  GEN->dim = gen->distr->dim; 
  gen->genid = _unur_set_genid(GENTYPE);
  SAMPLE = _unur_mvtdr_sample_cvec;
  gen->destroy = _unur_mvtdr_free;
  gen->clone = _unur_mvtdr_clone;
  GEN->n_steps = 0;   
  GEN->steps_min = _unur_max( 0, PAR->steps_min );  
  if ( (1 << (GEN->dim + GEN->steps_min)) > PAR->max_cones) {
    PAR->max_cones = 1 << (GEN->dim + GEN->steps_min);
  }
  GEN->max_gamma = UNUR_INFINITY;       
  GEN->cone = NULL;
  GEN->last_cone = NULL;
  GEN->n_cone = 0;                      
  GEN->max_cones = PAR->max_cones;      
  GEN->bound_splitting = PAR->bound_splitting;    
  GEN->vertex = NULL;
  GEN->last_vertex = NULL;
  GEN->n_vertex = 0;                    
  GEN->etable = NULL;                   
  GEN->etable_size = 0;                 
  GEN->guide = NULL;
  GEN->guide_size = 0;
  GEN->S         = malloc( GEN->dim * sizeof(double) );
  GEN->g         = malloc( GEN->dim * sizeof(double) );
  GEN->tp_coord  = malloc( GEN->dim * sizeof(double) );
  GEN->tp_mcoord = malloc( GEN->dim * sizeof(double) );
  GEN->tp_Tgrad  = malloc( GEN->dim * sizeof(double) );
  if (GEN->S==NULL || GEN->g==NULL || GEN->tp_coord==NULL || 
      GEN->tp_mcoord==NULL || GEN->tp_Tgrad==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,"");
    _unur_mvtdr_free(gen); return NULL;
  }
  GEN->center = unur_distr_cvec_get_center(gen->distr);
  GEN->pdfcenter = PDF(GEN->center);
  GEN->has_domain = (gen->distr->set & UNUR_DISTR_SET_DOMAIN) ? TRUE : FALSE;
#ifdef UNUR_ENABLE_INFO
  gen->info = _unur_mvtdr_info;
#endif
  return gen;
} 
struct unur_gen *
_unur_mvtdr_clone( const struct unur_gen *gen )
{ 
#define CLONE  ((struct unur_mvtdr_gen*)clone->datap)
  struct unur_gen *clone;
  int error = FALSE;
  size_t size;
  VERTEX *vt, **vtindex;
  CONE *c;
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MVTDR_GEN,NULL);
  clone = _unur_generic_clone( gen, GENTYPE );
  CLONE->center = unur_distr_cvec_get_center(clone->distr);
  size = GEN->dim * sizeof(double);
  CLONE->S         = malloc(size);
  CLONE->g         = malloc(size);
  CLONE->tp_coord  = malloc(size);
  CLONE->tp_mcoord = malloc(size);
  CLONE->tp_Tgrad  = malloc(size);
  vtindex = malloc(GEN->n_vertex * sizeof (VERTEX *));
  if (CLONE->S==NULL || CLONE->g==NULL || CLONE->tp_coord==NULL || 
      CLONE->tp_mcoord==NULL || CLONE->tp_Tgrad==NULL || vtindex==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,"");
    if (vtindex) free (vtindex);
    _unur_mvtdr_free(clone); return NULL;
  }
  if (GEN->S) memcpy( CLONE->S, GEN->S, size );
  if (GEN->g) memcpy( CLONE->g, GEN->g, size );
  if (GEN->tp_coord) memcpy( CLONE->tp_coord, GEN->tp_coord, size );
  if (GEN->tp_mcoord) memcpy( CLONE->tp_mcoord, GEN->tp_mcoord, size );
  if (GEN->tp_Tgrad) memcpy( CLONE->tp_Tgrad, GEN->tp_Tgrad, size );
  CLONE->vertex = NULL;  CLONE->n_vertex = 0;
  CLONE->cone = NULL;    CLONE->n_cone = 0;
  CLONE->guide = NULL;
  for (vt = GEN->vertex; vt != NULL; vt = vt->next) {
    VERTEX *vtc = _unur_mvtdr_vertex_new( clone );
    if (vtc == NULL) {
      error = TRUE; break; }
    memcpy(vtc->coord, vt->coord, size);
    vtc->index = vt->index;
    vtindex[vt->index] = vtc;
  }
  for (c = GEN->cone; c != NULL && !error; c = c->next) {
    CONE *cc, *cc_next;
    VERTEX **v;
    double *center, *gv;
    int i;
    cc = _unur_mvtdr_cone_new( clone );
    if (cc == NULL) {
      error = TRUE; break; }
    cc_next = cc->next;
    center = cc->center;
    gv = cc->gv;
    v = cc->v;
    memcpy(cc,c,sizeof(CONE));
    memcpy(center, c->center, size);
    memcpy(gv, c->gv, size);
    for (i=0; i<GEN->dim; i++)
      v[i] = vtindex[(c->v[i])->index];
    cc->next = cc_next;
    cc->center = center;
    cc->gv = gv;
    cc->v = v;
  }
  if (_unur_mvtdr_make_guide_table(clone) != UNUR_SUCCESS)
    error = TRUE;
  free (vtindex);
  if (error == TRUE) {
    _unur_mvtdr_free(clone); return NULL;
  }
  return clone;
#undef CLONE
} 
void
_unur_mvtdr_free( struct unur_gen *gen )
{ 
  VERTEX *vt, *vt_next;
  CONE *c, *c_next;
  if( !gen ) 
    return;
  if ( gen->method != UNUR_METH_MVTDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);
  SAMPLE = NULL;   
  _unur_mvtdr_etable_free(gen);
  for (vt = GEN->vertex; vt != NULL; vt = vt_next) {
    vt_next = vt->next;
    free (vt->coord);    
    free (vt);
  }
  for (c = GEN->cone; c != NULL; c = c_next) {
      c_next = c->next;
      free (c->v);        
      free (c->center);   
      free (c->gv);       
      free (c);
  }
  if (GEN->guide) free (GEN->guide);
  if (GEN->S)         free (GEN->S);
  if (GEN->g)         free (GEN->g);
  if (GEN->tp_coord)  free (GEN->tp_coord);
  if (GEN->tp_mcoord) free (GEN->tp_mcoord);
  if (GEN->tp_Tgrad)  free (GEN->tp_Tgrad);
  _unur_generic_free(gen);
} 
struct unur_gen *
_unur_mvtdr_gammagen( struct unur_gen *gen, double alpha )
{
  struct unur_distr *gammadistr;
  struct unur_par   *gammapar;
  struct unur_gen   *gammagen;
  double shape;
  shape = alpha;
  gammadistr = unur_distr_gamma(&shape,1);
  if (_unur_isfinite(GEN->max_gamma)) {
    unur_distr_cont_set_domain(gammadistr,0.,GEN->max_gamma);
  }
  gammapar = unur_tdr_new( gammadistr );
  unur_tdr_set_usedars( gammapar, TRUE );
  unur_tdr_set_max_sqhratio( gammapar, MVTDR_TDR_SQH_RATIO);
  if (! GEN->has_domain) {
    unur_tdr_set_variant_ia( gammapar );
  }
  gammagen = unur_init( gammapar );
  _unur_distr_free( gammadistr );
  if (gammagen == NULL) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,
		"Cannot create aux Gamma generator");
    return NULL;
  }
  gammagen->urng = gen->urng;
  gammagen->debug = gen->debug;
  return gammagen;
} 
int
_unur_mvtdr_create_hat( struct unur_gen *gen )
{
  int step;            
  double Hi_bound;     
  CONE *c;
  int n_splitted;
  if( _unur_mvtdr_initial_vertices(gen) != UNUR_SUCCESS ) 
    return UNUR_FAILURE;
  if( _unur_mvtdr_initial_cones(gen) != UNUR_SUCCESS ) 
    return UNUR_FAILURE;
  for( step = 1; step <= GEN->steps_min; step++ ) {
    if (_unur_mvtdr_triangulate(gen,step,TRUE) < 0)
      return UNUR_FAILURE;
  }
  for( c = GEN->cone; c != NULL; c = c->next )
    _unur_mvtdr_tp_find (gen,c);
  while( _unur_mvtdr_triangulate(gen,step,FALSE) > 0 ) {
    if (GEN->n_cone > GEN->max_cones)
      return UNUR_FAILURE;
    step++;
  }
  GEN->n_steps = step-1;
  GEN->Htot = 0.;                 
  for( c=GEN->cone; c!=NULL; c=c->next ) {
    GEN->Htot += c->Hi;           
    c->Hsum = GEN->Htot;          
  }
  while (1) {
    Hi_bound = GEN->bound_splitting * GEN->Htot / GEN->n_cone;
    GEN->Htot = 0.;
    n_splitted = 0;
    for( c=GEN->cone; c!=NULL; c=c->next ) {   
      while( Hi_bound < c->Hi && GEN->n_cone < GEN->max_cones ) { 
	if (_unur_mvtdr_cone_split(gen,c,c->level+1) != UNUR_SUCCESS)
	  return UNUR_FAILURE;
	++n_splitted;
	_unur_mvtdr_tp_find (gen,c);
	_unur_mvtdr_tp_find (gen,GEN->last_cone);
      }
      GEN->Htot += c->Hi;           
      c->Hsum = GEN->Htot;          
      if( c == GEN->last_cone ) break;
    }
    if (!n_splitted || GEN->n_cone >= GEN->max_cones) break;
  }
  if (_unur_mvtdr_make_guide_table(gen) != UNUR_SUCCESS) 
    return UNUR_FAILURE;
  if (GEN->dim > 2)
    _unur_mvtdr_etable_free(gen);
  return UNUR_SUCCESS;
} 
int
_unur_mvtdr_initial_cones( struct unur_gen *gen )
{
  int i,j,k;
  CONE *c;
  int max_c;       
  VERTEX *vt;
  VERTEX **ivtl;   
  int dim = GEN->dim;
  int error = FALSE;
  int have_negative_index = FALSE;
  int cone_out_of_domain;
  ivtl = malloc(2 * dim * sizeof(VERTEX*));
  if( ivtl==NULL ) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return UNUR_ERR_MALLOC; 
  }
  for (vt = GEN->vertex, i=0; i < 2*GEN->dim && vt!=NULL; vt = vt->next, i++) {
    if (vt->index<0) have_negative_index = TRUE;
    ivtl[i] = vt;
  }
  max_c = 1 << dim;
  for( k=0; k<max_c; k++ ) {
    if (have_negative_index) {
      cone_out_of_domain = FALSE;
      for( i=0; i < dim; i++ ) {
	if ( (!((k>>i)&1) && ivtl[i]->index < 0) ||
	     ( ((k>>i)&1) && ivtl[i+dim]->index < 0) ) {
	  cone_out_of_domain = TRUE;  break;
	}
      }
      if (cone_out_of_domain)  
	continue;
    }
    c = _unur_mvtdr_cone_new(gen);
    if (c==NULL) { error = TRUE; break; }
    c->level = 0;
    j = 0;
    for( i=0; i < dim; i++ )
      if (!((k>>i)&1))  (c->v)[j++] = ivtl[i];
    for( i=0; i < dim && j < dim; i++ )
      if ( ((k>>i)&1))  (c->v)[j++] = ivtl[i + dim];
    c->logdetf = 0.;
    c->tp = -1.;   
  }
  free (ivtl);
  if (error==TRUE) return UNUR_ERR_MALLOC;
  return UNUR_SUCCESS;
} 
CONE *
_unur_mvtdr_cone_new( struct unur_gen *gen )
{
  CONE *c; 
  c = malloc(sizeof(CONE));
  if (c==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }
  if (GEN->cone == NULL)
    GEN->last_cone = GEN->cone = c;
  else
    GEN->last_cone = GEN->last_cone->next = c;
  c->next = NULL;
  c->v      = malloc( GEN->dim * sizeof(VERTEX *));
  c->center = malloc( GEN->dim * sizeof(double));
  c->gv     = malloc( GEN->dim * sizeof(double));
  if (c->v==NULL || c->center==NULL || c->gv==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }
  c->height = UNUR_INFINITY;
  c->tp = -1.;
  c->Hi = UNUR_INFINITY;
  ++(GEN->n_cone);
  return c;
} 
int
_unur_mvtdr_cone_center( struct unur_gen *gen, CONE *c )
{
  int i,k;
  double norm;
  int dim = GEN->dim;
  norm = 0.;
  for( i=0; i<dim; i++ ) {
    c->center[i] = 0.;
    for( k=0; k<dim; k++ )
      c->center[i] += (c->v[k])->coord[i];       
    norm += c->center[i] * c->center[i];         
  }
  norm = sqrt(norm);
  for( i=0; i<dim; i++ )
    c->center[i] /= norm;
  return UNUR_SUCCESS;
} 
int
_unur_mvtdr_cone_params( struct unur_gen *gen, CONE *c )
{
  double Tf,f;                  
  double Tderf;                 
  int i;                        
  int dim = GEN->dim;           
  double *g = GEN->g;               
  double *coord = GEN->tp_coord;    
  double *mcoord = GEN->tp_mcoord;  
  double *Tgrad = GEN->tp_Tgrad;    
  double tolerance = TOLERANCE * GEN->pdfcenter / dim;
  for( i=0; i<dim; i++ ) {
    coord[i] = c->tp * c->center[i];
    mcoord[i] = coord[i] + GEN->center[i];
  }
  if( DISTR.logpdf != NULL ) {
    c->Tfp = Tf = logPDF(mcoord);
    if (! _unur_isfinite(Tf))
      return UNUR_ERR_DISTR_DOMAIN;
  }
  else {
    f = PDF(mcoord);
    if( f < tolerance )
      return UNUR_ERR_DISTR_DOMAIN;
    c->Tfp = Tf = T(f);
  }
  if( DISTR.dlogpdf != NULL ) {
    dlogPDF(Tgrad,mcoord);
  }
  else {
    dPDF(Tgrad,mcoord);
    Tderf = T_deriv(exp(Tf));
    for( i=0; i<dim; i++ )
      Tgrad[i] *= Tderf;
  }
  c->alpha = Tf - _unur_vector_scalar_product(dim,Tgrad,coord);
  c->beta = _unur_vector_norm(dim,Tgrad);
  if( c->beta < tolerance )
    return UNUR_FAILURE;
  for( i=0; i<dim; i++ )
    g[i] = - Tgrad[i] / c->beta;
  c->logai = c->logdetf;
  for( i=0; i<dim; i++ ) {
    c->gv[i] = _unur_vector_scalar_product(dim,g,(c->v[i])->coord);   
    if( c->gv[i] < tolerance )
      return UNUR_FAILURE;
    else
      c->logai -= log(c->gv[i]);
  }
  if (_unur_mvtdr_cone_height(gen,c) != UNUR_SUCCESS) 
    return UNUR_FAILURE;
  return UNUR_SUCCESS;
} 
double
_unur_mvtdr_cone_logH( struct unur_gen *gen, CONE *c )
{
  double logH;
  switch ( _unur_mvtdr_cone_params(gen,c) ) {
  case UNUR_SUCCESS:
    break;
  case UNUR_ERR_DISTR_DOMAIN:
    return -UNUR_INFINITY;
  default:
    return UNUR_INFINITY;
  }
  logH = c->alpha - GEN->dim * log(c->beta) + c->logai;
  if (_unur_isfinite(c->height)) {
    if (c->height < 1.e-50) 
      return -UNUR_INFINITY;
    else
      logH += log(_unur_SF_incomplete_gamma(c->beta*c->height,(double)GEN->dim));
  }
  return (_unur_isfinite(logH)) ? logH : UNUR_INFINITY;
} 
int
_unur_mvtdr_cone_split( struct unur_gen *gen, CONE *c, int step )
{
  CONE *newc;       
  VERTEX *newv;      
  int dim = GEN->dim;  
  int i;
  if (dim == 2)
    newv = _unur_mvtdr_vertex_on_edge(gen,c->v);
  else
    newv = _unur_mvtdr_etable_find_or_insert(gen,c->v);
  if (newv==NULL) return UNUR_FAILURE;
  newc = _unur_mvtdr_cone_new(gen);   
  if (newc==NULL) return UNUR_ERR_MALLOC;
  newc->level = step;                 
  for (i=0; i<dim-1; i++)
    newc->v[i] = c->v[i+1];           
  newc->v[dim-1] = newv;              
  newc->logdetf = c->logdetf - log(2.*newv->norm);  
  newc->tp = c->tp;                   
  c->level = step;                    
  for (i=0; i<dim-2; i++)
    c->v[i+1] = c->v[i+2];            
  (c->v)[dim-1] = newv;               
  c->logdetf = newc->logdetf;         
  GEN->n_steps = _unur_max(GEN->n_steps, step); 
  return UNUR_SUCCESS;
} 
int
_unur_mvtdr_triangulate( struct unur_gen *gen, int step, int all )
{
  int k,nc;
  CONE *c;
  int dim = GEN->dim;  
  if (dim > 2) {
    if( step % (dim-1) == 1 )
      if( _unur_mvtdr_etable_new(gen, _unur_mvtdr_number_vertices(gen, (step/(dim-1)+1)*(dim-1) ))
	  != UNUR_SUCCESS )
	return -1;
  }
  nc = GEN->n_cone;
  for( k=0, c=GEN->cone; k<nc; k++ ) {
    if( all ) {
      if (_unur_mvtdr_cone_split(gen,c,step) != UNUR_SUCCESS)
	return -1;
    }
    else if ( c->tp < 0. ) {
      if (_unur_mvtdr_cone_split(gen,c,step) != UNUR_SUCCESS)
	return -1;
      _unur_mvtdr_tp_find (gen,c);
      _unur_mvtdr_tp_find (gen,GEN->last_cone);
    }
    c = c->next;
  }
  return (GEN->n_cone - nc);
} 
int
_unur_mvtdr_cone_height( struct unur_gen *gen, CONE *c )
{
#define A(i,j)  (AA[(dim+1)*(i)+(j)])
  double *AA;
  int dim = GEN->dim;
#define ll(i)   (domain[2*(i)]   - GEN->center[(i)])
#define ur(i)   (domain[2*(i)+1] - GEN->center[(i)])
  double *domain;
  int i,j,row,ipc,ipr;
  double pc,pr,ratio;
  double sgn;
  if (! GEN->has_domain)
    return UNUR_SUCCESS;
  if (DISTR.domainrect == NULL) {
    _unur_error(gen->genid,UNUR_ERR_DISTR_DOMAIN,"no domain given");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  domain = DISTR.domainrect;
  AA = _unur_xmalloc( (dim+1)*(dim+1)*sizeof(double) );
  for( i=0, row=0; i<dim; i++ ) {
    for( j=0, sgn=0.; j<dim; j++ ) {
      if( (c->v[j])->coord[i] > 0. ) {
	sgn = +1.; break;
      }
      if( (c->v[j])->coord[i] < 0. ) {
	sgn = -1.; break;
      }
    }
    if (_unur_iszero(sgn)) continue;
    for( j=0; j<dim; j++ )
      A(row,j) = sgn * (c->v[j])->coord[i];
    A(row,dim) = (sgn > 0.) ? ur(i) : -(ll(i));
    row++;
  }
  for( j=0; j<dim; j++ )
    A(row,j) = -(c->gv[j]);
  A(row,dim) = 0.;
  while( 1 ) {
    for( j=0,pc=0.,ipc=-1; j<dim; j++) {
      if( A(row,j) < pc ) {
        ipc = j;
        pc = A(row,ipc);
      }
    }
    if( ipc == -1 ) {   
      c->height = A(row,dim);
      break;
    }
    for( i=0,pr=-1.,ipr=-1; i<row; i++ ) {
      if( A(i,ipc) <= 0. ) continue;
      ratio = A(i,dim) / A(i,ipc);
      if( pr < 0 || pr > ratio ) {
        ipr = i;
        pr = ratio;
      }
    }
    if( ipr == -1 ) {
      c->height = UNUR_INFINITY;
      break;
    }
    for( i=0; i<=row; i++ )
      if( i!= ipr )
        for( j=0; j<dim+1; j++ )
          if( j!= ipc ) 
            A(i,j) -= A(ipr,j) * A(i,ipc) / A(ipr,ipc);
    for( i=0; i<=row; i++ )
      if( i != ipr )
        A(i,ipc) /= -A(ipr,ipc);
    for( j=0; i<dim; i++ )
      if( j != ipc )
        A(ipr,j) /= A(ipr,ipc);
    A(ipr,ipc) = 1./A(ipr,ipc);
  }
  free (AA);
  if (_unur_isnan(c->height))  c->height = UNUR_INFINITY;
  return UNUR_SUCCESS;
#undef A
#undef ll
#undef ur
} 
int
_unur_mvtdr_max_gamma( struct unur_gen *gen )
{
  double max, tmp;
  CONE *c;
  if (!GEN->has_domain) {
    GEN->max_gamma = UNUR_INFINITY;
  }
  else {
    max = 0.;
    for (c = GEN->cone; c != NULL; c = c->next) {
      tmp = c->height * c->beta;
      max = _unur_max(max, tmp);
    }
    GEN->max_gamma = (max > 0.) ? max : UNUR_INFINITY;
  }
  return UNUR_SUCCESS;
} 
double
_unur_mvtdr_tp_min (double t, void *p )
{
  TP_ARG *a = p;
  (a->c)->tp = a->t = t;
  a->logH = _unur_mvtdr_cone_logH (a->gen, a->c);
  switch (_unur_isinf(a->logH)) {
  case -1:
    a->logH = UNUR_INFINITY;
    a->status = MVTDR_CONE_DOMAIN;
    break;
  case 1:
    a->status = MVTDR_CONE_INVALID;
    break;
  case 0:
  default:
    a->status = MVTDR_CONE_OK;
  }
  if( a->status != MVTDR_CONE_OK )
    (a->c)->tp = -1.;
  return a->logH;
} 
double
_unur_mvtdr_tp_min_aux(double t, void *p)
{
  return (- _unur_mvtdr_tp_min(t, p) );
}
int
_unur_mvtdr_tp_find( struct unur_gen *gen, CONE *c )
{
  struct unur_funct_generic tpaux;
  TP_ARG a[3];    
  int i;
  _unur_mvtdr_cone_center(gen,c);
  for (i=0; i<3; i++) { a[i].c = c; a[i].gen = gen; }
  switch (_unur_mvtdr_tp_search(gen,a)) {
  case UNUR_SUCCESS:
    break;
  case UNUR_ERR_DISTR_DOMAIN:
    c->tp = 0.;
    c->Hi = 0.;
    c->height = 0.;
    return UNUR_ERR_DISTR_DOMAIN;
  case UNUR_FAILURE:
  default:
    return UNUR_FAILURE;
  }
  switch( _unur_mvtdr_tp_bracket(gen,a) ) {      
  case TP_BRACKET:                 
    tpaux.f = _unur_mvtdr_tp_min_aux;
    tpaux.params = a+1;
    c->tp = _unur_util_brent( tpaux, a[0].t, a[2].t, a[1].t, FIND_TP_TOL);
    c->Hi = exp(a[1].logH);
    break;                         
  case TP_LEFT:                    
    c->tp = a[0].t;
    c->Hi = exp(a[0].logH);
    break;
  case TP_MIDDLE:                  
    _unur_mvtdr_tp_min(a[1].t, a+1);  
    c->Hi = exp(a[1].logH);
    break;
  case TP_RIGHT:                   
    c->tp = a[2].t;
    c->Hi = exp(a[2].logH);
    break;
  default:                         
    c->tp = -1.;
    return UNUR_FAILURE;
  }
  return UNUR_SUCCESS;
} 
int
_unur_mvtdr_tp_search( struct unur_gen *gen ATTRIBUTE__UNUSED, TP_ARG *a )
{
#define N_STEPS  (10)  
  int i;                           
  int is_unbounded_cone = FALSE;  
  double start = 1.;               
  a[0].t = 0.;     
  a[1].t = start;  
  a[2].t = -1.;    
  for( i=1; i <= N_STEPS; i++ ) {
    _unur_mvtdr_tp_min(a[1].t, a+1);  
    if (a[1].status == MVTDR_CONE_OK) {
      return UNUR_SUCCESS;
    }
    else if (a[1].status == MVTDR_CONE_DOMAIN) {
      break;
    }
    else {
      is_unbounded_cone = TRUE;
      a[0].t = a[1].t;
      a[1].t *= 2.;
    }
  }
  a[0].t = 0.;        
  a[1].t = start/2.;  
  a[2].t = 1.;        
  for( i=0;; i++ ) {
    _unur_mvtdr_tp_min(a[1].t, a+1);  
    if (a[1].status == MVTDR_CONE_OK) {
      return UNUR_SUCCESS;
    }
    else if (a[1].status == MVTDR_CONE_DOMAIN) {
      if (a[1].t < 1.e-20 ) {
	return (is_unbounded_cone) ? UNUR_FAILURE : UNUR_ERR_DISTR_DOMAIN;
      }
      a[2].t = a[1].t;
      a[1].t /= 10.;
    }
    else {
      if (i > N_STEPS) {
	return UNUR_FAILURE;
      }
      is_unbounded_cone = TRUE;
      a[2].t = a[1].t;
      a[1].t /= 2.;
    }
  }
#undef N_STEPS
} 
int 
_unur_mvtdr_tp_bracket( struct unur_gen *gen ATTRIBUTE__UNUSED, TP_ARG *a )
{
#define N_STEPS  (10)  
  int i;                 
  double tleft, tright;  
  tleft = a[0].t;
  a[0].t = a[1].t / 2.;
  for( i=1; i <= _unur_max(1,N_STEPS); i++ ) {
    _unur_mvtdr_tp_min(a[0].t, a);  
    if( a[0].status != MVTDR_CONE_OK ) {
      tleft = a[0].t;                     
      a[0].t += (a[1].t - a[0].t) / 2.;   
    }
    else if( a[0].logH <= a[1].logH ) {
      a[2].t = a[1].t; a[2].logH = a[1].logH; a[2].status = MVTDR_CONE_OK;
      a[1].t = a[0].t; a[1].logH = a[0].logH; a[1].status = MVTDR_CONE_OK;
      a[0].t = tleft + (a[0].t - tleft)*0.5;
    }
    else  
      break;
  }
  if( a[0].status != MVTDR_CONE_OK )
    return TP_MIDDLE;
  if( a[0].logH <= a[1].logH )
    return TP_LEFT;
  if( a[2].t < 0. )
    a[2].t = 1.1 * a[1].t;
  tright = -1.;    
  tleft = a[1].t;
  for( i=1; i <= _unur_max(1,N_STEPS); i++ ) {
    _unur_mvtdr_tp_min(a[2].t, a+2);  
    if( a[2].status != MVTDR_CONE_OK ) {
      tright = a[2].t;
      a[2].t = (tleft + a[2].t) * 0.5;   
    }
    else if( a[2].logH <= a[1].logH ) {
      tleft = a[2].t;
      a[2].t = (tright < 0.) ? a[2].t * 2. : (tright + a[2].t) * 0.5;
    }
    else    
      break;
  }
  if( a[2].status != MVTDR_CONE_OK )
    return TP_MIDDLE; 
  if( a[2].logH <= a[1].logH )
    return TP_RIGHT;
  return TP_BRACKET;
#undef N_STEPS
} 
int
_unur_mvtdr_initial_vertices( struct unur_gen *gen )
{
  VERTEX *vt;
  int i,k;
  double d;
  double *domain;
  domain = ( (GEN->has_domain && DISTR.domainrect != NULL) 
	     ? DISTR.domainrect : NULL );
  for ( d=1.; d > -2.; d -= 2.) {
    for( k=0; k<GEN->dim; k++ ) {
      vt = _unur_mvtdr_vertex_new(gen);
      if (vt==NULL) return UNUR_FAILURE;
      for( i=0; i<GEN->dim; i++ ) {
	(vt->coord)[i] = (i==k) ? d : 0.;
      }
      vt->norm = 1.;
      if ( domain ) {
	if ( (d<0 && _unur_FP_equal(GEN->center[k],domain[2*k]))   ||
	     (d>0 && _unur_FP_equal(GEN->center[k],domain[2*k+1])) )
	  vt->index = -vt->index -1;
      }
    }
  }
  return UNUR_SUCCESS;
} 
VERTEX *
_unur_mvtdr_vertex_new( struct unur_gen *gen )
{
  VERTEX *v;
  v = malloc(sizeof(VERTEX));
  if (v==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }
  if (GEN->vertex == NULL) {
    GEN->last_vertex = GEN->vertex = v;
  }
  else {
    GEN->last_vertex = GEN->last_vertex->next = v;
  }
  v->next = NULL;
  v->coord = malloc(GEN->dim * sizeof(double));
  if (v->coord==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }
  v->index = GEN->n_vertex;
  ++(GEN->n_vertex);
  return GEN->last_vertex;
} 
VERTEX *
_unur_mvtdr_vertex_on_edge( struct unur_gen *gen, VERTEX **vl )
{
  int i;
  VERTEX *newv;          
  newv = _unur_mvtdr_vertex_new(gen);
  if (newv==NULL) return NULL;
  for( i=0; i<GEN->dim; i++ )
    newv->coord[i] =
      0.5 * ( ((vl[0])->coord)[i] + ((vl[1])->coord)[i] );
  newv->norm = _unur_vector_norm(GEN->dim, newv->coord);
  for( i=0; i<GEN->dim; i++ )
    newv->coord[i] /= newv->norm;
  return newv;
} 
int
_unur_mvtdr_number_vertices( struct unur_gen *gen, int level )
{
  if (level < 0 || GEN->dim < 2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return -1;  
  }
  switch (GEN->dim) {
  case 2: {
    return (1 << (level+2));
  }
  case 3: {
    static int nv[]={ 6, 13, 18, 40, 66,142,258,538,1026,2098,4098,8290,16386,32962,65538,131458,262146};
    return nv[_unur_min(level,16)];
  }
  case 4: {
    static int nv[]={ 8, 19, 25, 32, 80,128,192,456, 824,1408,3120,5968,11008,23264,45600, 87552};
    return nv[_unur_min(level,15)];
  }
  case 5: {
    static int nv[]={10, 26, 33, 41, 50,140,220,321, 450,1186,2158,3636, 5890,13970,27130};
    return nv[_unur_min(level,14)];
  }
  case 6: {
    static int nv[]={12, 34, 42, 51, 61, 72,224,348, 501, 681, 912,2660, 4896, 8254};
    return nv[_unur_min(level,13)];
  }
  case 7: {
    static int nv[]={14, 43, 52, 62, 73, 85, 98,336, 518, 743, 985,1289, 1666};
    return nv[_unur_min(level,12)];
  }
  case 8: {
    static int nv[]={16, 53, 63, 74, 86, 99,113,128, 480, 736,1059};
    return nv[_unur_min(level,10)];
  }
  case 9: {
    static int nv[]={18, 64, 75, 87,100,114,129,145, 162, 660};
    return nv[_unur_min(level,9)];
  }
  case 10: {
    static int nv[]={20, 76, 88,101,115,130,146,163, 181, 200};
    return nv[_unur_min(level,9)];
  }
  case 11: {
    static int nv[]={22, 89,102,116,131,147,164,182, 201, 221, 242};
    return nv[_unur_min(level,10)];
  }
  default: { 
    static int nv[]={24,103,117,132,148,165,183,202, 222, 243, 265, 288};
    return nv[_unur_min(level,11)];
  }
  }
} 
#define _unur_mvtdr_etable_hash(x,y)  ( (3*((x)+(y))/2) % GEN->etable_size )
int
_unur_mvtdr_etable_new( struct unur_gen *gen, int size )
{
  int n;
  _unur_mvtdr_etable_free(gen);
  GEN->etable_size = size;
  GEN->etable = _unur_xmalloc( size * sizeof(E_TABLE*) );
  if (GEN->etable==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return UNUR_ERR_MALLOC; }
  for (n = 0; n< size; n++) 
    GEN->etable[n] = NULL;
  return UNUR_SUCCESS;
} 
void
_unur_mvtdr_etable_free( struct unur_gen *gen )
{
  int i;
  E_TABLE *et, *et_next;
  if( GEN->etable == NULL )
    return;
  for( i=0; i<GEN->etable_size; i++ ) {
    for (et = GEN->etable[i]; et != NULL; et = et_next) {
      et_next = et->next;
      free (et);
    }
  }
  free( GEN->etable );
  GEN->etable = NULL;
  GEN->etable_size = 0;
} 
VERTEX *
_unur_mvtdr_etable_find_or_insert( struct unur_gen *gen, VERTEX **vidx )
{
  E_TABLE *pet, *pet_last;  
  int idx[2];               
  int hidx;                 
  CHECK_NULL(GEN->etable,NULL);
  idx[0] = vidx[0]->index;
  idx[1] = vidx[1]->index;
  hidx = _unur_mvtdr_etable_hash(idx[0],idx[1]);
  pet = pet_last = *(GEN->etable + hidx);
  while( pet != NULL ) {
    if( pet->index[0] == idx[0] && pet->index[1] == idx[1] )
      break;   
    pet_last = pet;
    pet =  pet->next;
  }
  if( pet == NULL ) {
    pet = malloc( sizeof(E_TABLE) );
    if (pet==NULL) {
      _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return NULL; }
    pet->next = NULL;
    if (pet_last == NULL)
      *(GEN->etable + hidx) = pet;
    else
      pet_last->next = pet;
    pet->index[0] = idx[0];
    pet->index[1] = idx[1];
    pet->vertex = _unur_mvtdr_vertex_on_edge(gen,vidx);
  }
  return pet->vertex;
} 
int
_unur_mvtdr_make_guide_table( struct unur_gen *gen )
{
  int j;
  CONE *c;
  GEN->guide_size = GEN->n_cone * GUIDE_TABLE_SIZE;
  GEN->guide = malloc (GEN->guide_size * sizeof(CONE*));
  if (GEN->guide==NULL) {
    _unur_error(gen->genid,UNUR_ERR_MALLOC,""); return UNUR_ERR_MALLOC; }
  for( j = 0; j < GEN->guide_size ; j++ )
    GEN->guide[j] = NULL;
  for( c=GEN->cone, j=0; c!=NULL && j<GEN->guide_size; j++ ) {
    while( c->Hsum / GEN->Htot < (double) j / GEN->guide_size )
      c=c->next;
    (GEN->guide)[j] = c;
    if( c == GEN->last_cone ) break;
  }
  if( j<GEN->guide_size )
    for( ; j<GEN->guide_size; j++ )
      (GEN->guide)[j] = GEN->last_cone;
  return UNUR_SUCCESS;
} 
