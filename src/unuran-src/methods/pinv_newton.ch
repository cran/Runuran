/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

int
_unur_pinv_create_table( struct unur_gen *gen )
{
  double utol;               
  double maxerror;            
  double h;                  
  int i;                     
  int iter;                  
  int cont;                  
  int use_linear;            
  int right_bd;              
  double chebyshev[MAX_ORDER+1]; 
  double xval[MAX_ORDER+1];  
  int n_decr_h = 0;          
  int n_incr_h = 0;          
  int n_use_linear = 0;      
  int k;                     
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);
  _unur_pinv_chebyshev_points(GEN->order,chebyshev);
  utol = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION;
  h = (GEN->bright-GEN->bleft)/128.;
  if (_unur_pinv_interval( gen, 0, GEN->bleft, 0.) != UNUR_SUCCESS) 
    return UNUR_ERR_GEN_CONDITION;
  i = 0;                
  cont = TRUE;          
  use_linear = FALSE;   
  for (iter=0; cont ; iter++) {
    if (iter >= PINV_MAX_ITER_IVS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "maximum number of iterations exceeded");
      return UNUR_ERR_GEN_CONDITION;
    }
    if(GEN->iv[i].xi+h > GEN->bright) {
      h = GEN->bright - GEN->iv[i].xi;
      cont = FALSE;    
      right_bd = TRUE; 
    }
    else {
      right_bd = FALSE;
    }
    for(k=0; k<=GEN->order; k++)
      xval[k] = GEN->iv[i].xi + h * chebyshev[k];
    if (!use_linear) {
      if (_unur_pinv_newton_create(gen,&(GEN->iv[i]),xval) != UNUR_SUCCESS)
	use_linear = TRUE;
    }
    if (use_linear) {
      ++n_use_linear;
      if (_unur_pinv_linear_create(gen,&(GEN->iv[i]),xval) != UNUR_SUCCESS) {
	if (i==0) { 
	  GEN->bleft = GEN->iv[i].xi + h;
	  GEN->iv[i].xi = GEN->bleft;
	  continue;  
	}
	else if (right_bd) { 
	  GEN->bright = GEN->iv[i].xi;
	  cont = FALSE;
	  break;  
	}
	else {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PDF too close to 0 on relevant part of domain --> abort");
	  return UNUR_ERR_GEN_CONDITION;
	}
      }
    }
    maxerror = _unur_pinv_newton_maxerror(gen,&(GEN->iv[i]),xval);
    if (!(maxerror <= utol)) {
      h *= (maxerror > 4.*utol) ? 0.81 : 0.9;
      cont = TRUE;  
      ++n_decr_h;
    }
    else {
      if ( _unur_pinv_interval( gen, i+1, GEN->iv[i].xi+h, 
				GEN->iv[i].cdfi +(GEN->iv)[i].ui[GEN->order-1])
      	   != UNUR_SUCCESS )
	return UNUR_ERR_GEN_CONDITION;
      if (maxerror < 0.3*utol) {
	h *= (maxerror < 0.1*utol) ? 2. : 1.2;
	++n_incr_h;
      }
      i++;
      use_linear = FALSE;
    }
  }
  _unur_pinv_lastinterval(gen);
  GEN->Umax = GEN->iv[GEN->n_ivs].cdfi;
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_create_table(gen,iter,n_incr_h,n_decr_h,n_use_linear);
#endif
  return UNUR_SUCCESS;
}  
int 
_unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx )
{
  struct unur_pinv_interval *iv;
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  if (i >= GEN->max_ivs) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"maximum number of intervals exceeded");
    return UNUR_ERR_GEN_CONDITION;
  }
  iv = GEN->iv+i;     
  iv->xi = x;         
  iv->cdfi = cdfx;    
  COOKIE_SET(iv,CK_PINV_IV);
  iv->ui = _unur_xmalloc( GEN->order * sizeof(double) );
  iv->zi = _unur_xmalloc( GEN->order * sizeof(double) );
  GEN->n_ivs = i;
  _unur_lobatto_find_linear(GEN->aCDF,x);
  return UNUR_SUCCESS;
} 
int 
_unur_pinv_lastinterval( struct unur_gen *gen )
{
  double *ui, *zi;
  int i;
  struct unur_pinv_interval *last_iv = GEN->iv + GEN->n_ivs;
  ui = last_iv->ui;
  zi = last_iv->zi;
  GEN->iv = _unur_xrealloc( GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  for (i=0; i<GEN->order; i++) {
    ui[i] = 0.;
    zi[i] = 0.;
  }
  return UNUR_SUCCESS;
} 
int
_unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
			  double *xval)
{
  double fx = -1.;       
  double *ui = iv->ui;   
  double *zi = iv->zi;   
  double xi, dxi;        
  double area;           
  int i,k;               
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);
  for(i=0; i<GEN->order; i++) {
    xi = xval[i];
    dxi = xval[i+1]-xval[i];
    area = _unur_pinv_Udiff(gen, xi, dxi, &fx);
    if (_unur_iszero(area)) return UNUR_ERR_SILENT;
    ui[i] = (i>0) ? (ui[i-1]+area) : area;
    zi[i] = dxi/area;
  }
  for(k=1; k<GEN->order; k++) {
    for(i=GEN->order-1; i>k; i--) {
      zi[i] = (zi[i]-zi[i-1]) / (ui[i]-ui[i-(k+1)]);
    }
    zi[k] = (zi[k]-zi[k-1]) / ui[k];
  }
  for (i=0; i<GEN->order; i++) {
    if (!_unur_isfinite(zi[i])) 
      return UNUR_ERR_SILENT;
  }
  return UNUR_SUCCESS;
} 
int
_unur_pinv_linear_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
			  double *xval)
{
  double *ui = iv->ui;   
  double *zi = iv->zi;   
  double x0, x1;         
  double area;           
  int i;                 
  x0 = xval[0];
  x1 = xval[GEN->order];
  for (i=1; i<GEN->order; i++) {
    ui[i] = zi[i] = 0.;
  }
  area = _unur_pinv_Udiff(gen, x0, x1-x0, NULL);
  ui[0] = area;
  zi[0] = (x1 - x0) / area;
  ui[GEN->order-1] = area;
  if (!_unur_isfinite(zi[0]))
    return UNUR_ERR_GEN_CONDITION;
  return UNUR_SUCCESS;
} 
int
_unur_pinv_chebyshev_points (int order, double *pt)
{
  int i;
  double phi = M_PI*0.5/(order+1); 
  pt[0] = 0.;
  for(i=1; i<order; i++)
    pt[i] = sin(i*phi) * sin((i+1)*phi)/cos(phi);
  pt[order] = 1.;
  return UNUR_SUCCESS;
} 
double
_unur_pinv_newton_eval ( double q, double ui[], double zi[], int order )
{
  int k;
  double chi;
  chi = zi[order-1];
  for (k=order-2; k>=0; k--)
    chi = chi*(q-ui[k])+zi[k];
  return (chi*q);
} 
double
_unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv,
			    double xval[])
{
  double x0 = iv->xi;    
  double *ui = iv->ui;   
  double *zi = iv->zi;   
  double maxerror = 0.;  
  double uerror;         
  double x;              
  double u;              
  int is_linear;         
  double testu[MAX_ORDER];  
  int i;                 
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);
  is_linear = _unur_iszero(zi[1]) ? TRUE : FALSE;
  if (is_linear && zi[0] < 0.)
    return 1000.;
  if (is_linear) 
    _unur_pinv_linear_testpoints(GEN->order,ui,testu);
  else
    _unur_pinv_newton_testpoints(GEN->order,ui,testu);
  for(i=0; i<GEN->order; i++){
    x = _unur_pinv_newton_eval(testu[i], ui, zi, GEN->order);
    if (!is_linear)
      if (! (xval[i] <= x0+x && x0+x <=xval[i+1]) )
	return 1000.;
    if (i==0 || xval==NULL)
      u = _unur_pinv_Udiff(gen, x0, x, NULL);
    else
      u = ui[i-1] + _unur_pinv_Udiff(gen, xval[i], x+x0-xval[i], NULL);
    if (!_unur_isfinite(u))
      return INFINITY;
    uerror = fabs(u - testu[i]);
    if (uerror>maxerror) maxerror = uerror;
  }
  return maxerror;
} 
int
_unur_pinv_newton_testpoints (int order, double ui[], double utest[])
{
  int k,j,i;
  double sum, qsum,x;
  for(k=0; k<order; k++) {
    x = (k>0) ? 0.5*(ui[k-1]+ui[k]) : 0.5*ui[k];
    for(j=1; j<=2; j++) {
      sum = 1./x;
      qsum = sum*sum;
      for(i=0; i<order; i++){
	sum += 1./(x-ui[i]);
	qsum += 1./((x-ui[i])*(x-ui[i]));
      }
      x += sum/qsum;
    }
    utest[k] = x;
  }
  return UNUR_SUCCESS;
} 
int
_unur_pinv_linear_testpoints (int order, double ui[], double utest[])
{
  int k;
  double dx = ui[order-1] / order;  
  for(k=0; k<order; k++)
    utest[k] = (k+0.5) * dx;
  return UNUR_SUCCESS;
} 
