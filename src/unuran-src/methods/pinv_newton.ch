/* Copyright (c) 2000-2012 Wolfgang Hoermann and Josef Leydold */
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
  int smooth;                
  int use_linear;            
  int right_bd;              
  int use_upoints;           
  double chebyshev[3][MAX_ORDER+1]; 
  double xval[MAX_ORDER+1];  
  int n_decr_h = 0;          
  int n_incr_h = 0;          
  int n_use_linear = 0;      
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);
  utol = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION;
  h = (GEN->bright-GEN->bleft)/128.;
  if (_unur_pinv_interval( gen, 0, GEN->bleft, 0.) != UNUR_SUCCESS) 
    return UNUR_ERR_GEN_CONDITION;
  for (smooth=0; smooth<=GEN->smooth; ++smooth)
    _unur_pinv_chebyshev_points(chebyshev[smooth],GEN->order,smooth);
  i = 0;                
  cont = TRUE;          
  use_linear = FALSE;   
  use_upoints = FALSE;  
  for (iter=0; cont ; iter++) {
    if (iter >= PINV_MAX_ITER_IVS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "maximum number of iterations exceeded");
      return UNUR_ERR_GEN_CONDITION;
    }
    if(!_unur_FP_less(GEN->iv[i].xi+h,GEN->bright)) {
      h = GEN->bright - GEN->iv[i].xi;
      cont = FALSE;    
      right_bd = TRUE; 
    }
    else {
      right_bd = FALSE;
    }
    smooth = GEN->smooth;
    use_linear = FALSE;
    switch (smooth) {
    case 2:
      _unur_pinv_newton_cpoints(xval, GEN->order, GEN->iv+i, h, chebyshev[smooth], smooth, use_upoints);
      if (_unur_pinv_newton_create(gen,GEN->iv+i,xval) == UNUR_SUCCESS)
	break;
      smooth = 1;
    case 1:
      if (GEN->order % 2 == 1) {
	_unur_pinv_newton_cpoints(xval, GEN->order, GEN->iv+i, h, chebyshev[smooth], smooth, use_upoints);
	if (_unur_pinv_newton_create(gen,GEN->iv+i,xval) == UNUR_SUCCESS)
	  break;
      }
      smooth = 0;
    case 0:
    default:
      _unur_pinv_newton_cpoints(xval, GEN->order, GEN->iv+i, h, chebyshev[smooth], smooth, use_upoints);
      if (_unur_pinv_newton_create(gen,GEN->iv+i,xval) == UNUR_SUCCESS)
	break;
    case -1:
      use_linear = TRUE;
    }
    if (use_linear) {
      ++n_use_linear;
      if (_unur_pinv_linear_create(gen,GEN->iv+i,xval) != UNUR_SUCCESS) {
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
    maxerror = _unur_pinv_newton_maxerror(gen,&(GEN->iv[i]),xval,use_linear);
    if (!(maxerror <= utol)) {
      h *= (maxerror > 4.*utol) ? 0.81 : 0.9;
      cont = TRUE;  
      ++n_decr_h;
      use_upoints = FALSE;
      continue;
    }
    if (gen->variant & PINV_VARIANT_UPOINTS) {
      if (!use_upoints && !use_linear ) {
	use_upoints = TRUE;
	cont = TRUE;
	continue;
      }
      use_upoints = FALSE;
    }
    if ( _unur_pinv_interval( gen, i+1, GEN->iv[i].xi+h, 
			      GEN->iv[i].cdfi +(GEN->iv)[i].ui[GEN->order-1])
	 != UNUR_SUCCESS )
      return UNUR_ERR_GEN_CONDITION;
    if (maxerror < 0.3*utol) {
      h *= (maxerror < 0.1*utol) ? 2. : 1.2;
      ++n_incr_h;
    }
    i++;
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
_unur_pinv_chebyshev_points (double *pt, int order, int smooth)
{
  int i,k;
  int n_pt = (order+1)/(smooth+1); 
  double phi = M_PI*0.5/n_pt;     
  pt[0] = 0.;
  for (i=1; i<=order; i++) {
    k = i / (smooth+1);
    if (k<n_pt-1)
      pt[i] = sin(k*phi) * sin((k+1)*phi)/cos(phi);
    else
      pt[i] = 1.;
  }
  return UNUR_SUCCESS;
} 
int
_unur_pinv_newton_cpoints (double *xval, int order, struct unur_pinv_interval *iv,
			   double h, double *chebyshev, int smooth, int use_upoints)
{
  int k;
  if (use_upoints) {
    double hu = iv->ui[order-1];
    for(k=0; k<=order; k++)
      xval[k] = (k % (smooth+1))
       	? xval[k-1]
 	: iv->xi + _unur_pinv_newton_eval(hu * chebyshev[k], iv->ui, iv->zi, order);
  }
  else {
    for(k=0; k<=order; k++) {
      xval[k] = (k % (smooth+1))
      	? xval[k-1]
      	: iv->xi + h * chebyshev[k];
    }
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
    if (!_unur_FP_same(xval[i],xval[i+1])) {
      dxi = xval[i+1]-xval[i];
      area = _unur_pinv_Udiff(gen, xi, dxi, &fx);
      if (_unur_iszero(area)) return UNUR_ERR_SILENT;
      ui[i] = (i>0) ? (ui[i-1]+area) : area;
      zi[i] = dxi/area;
    }
    else {
      ui[i] = (i>0) ? ui[i-1] : 0.;
      zi[i] = 1./PDF(xi);
    }
  }
  for(i=GEN->order-1; i>=1; i--) {
    if (!_unur_FP_same(zi[i],zi[i-1]))
      zi[i] = (i>1) 
	? (zi[i]-zi[i-1]) / (ui[i]-ui[i-2])
	: (zi[1]-zi[0]) / ui[1];
    else
      zi[i] = (DISTR.dpdf != NULL) ? (-0.5 * dPDF(xval[i]) * pow(zi[i],3)) : INFINITY;
  }
  for(k=2; k<GEN->order; k++) {
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
double
_unur_pinv_newton_eval ( double q, double *ui, double *zi, int order )
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
			    double *xval, int is_linear)
{
  double x0 = iv->xi;    
  double *ui = iv->ui;   
  double *zi = iv->zi;   
  double maxerror = 0.;  
  double uerror;         
  double x;              
  double u;              
  double testu[MAX_ORDER];  
  int i;                 
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);
  if (is_linear && zi[0] < 0.)
    return 1001.;
  if (is_linear) 
    _unur_pinv_linear_testpoints(testu,ui,GEN->order);
  else
    _unur_pinv_newton_testpoints(testu,ui,GEN->order);
  for(i=0; i<GEN->order; i++) {
    x = _unur_pinv_newton_eval(testu[i], ui, zi, GEN->order);
    if (!is_linear) {
      if (! (xval[i] <= x0+x && x0+x <= xval[i+1]) )
	if (! _unur_FP_same(xval[i], xval[i+1]))
	  return 1002.;
    }
    if (i==0 || xval==NULL)
      u = _unur_pinv_Udiff(gen, x0, x, NULL);
    else
      u = ui[i-1] + _unur_pinv_Udiff(gen, xval[i], x+x0-xval[i], NULL);
    if (!_unur_isfinite(u))
      return INFINITY;
    uerror = fabs(u - testu[i]);
    if (uerror>maxerror) maxerror = uerror;
  }
  if (GEN->order == 3 && GEN->smooth==1 && xval!=NULL &&
      ! _unur_pinv_cubic_hermite_is_monotone(gen,ui,zi,xval))
    return 1003.;
  return maxerror;
} 
int
_unur_pinv_newton_testpoints (double *utest, double *ui, int order)
{
  int k,j,i;
  double sum, qsum,x;
  for(k=0; k<order; k++) {
    if ( (k==0 && _unur_iszero(ui[0])) ||
    	 (k>0  && _unur_FP_same(ui[k-1],ui[k])) ) {
      utest[k] = ui[k]; 
      continue;
    }
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
_unur_pinv_linear_testpoints (double *utest, double *ui, int order)
{
  int k;
  double dx = ui[order-1] / order;  
  for(k=0; k<order; k++)
    utest[k] = (k+0.5) * dx;
  return UNUR_SUCCESS;
} 
int
_unur_pinv_cubic_hermite_is_monotone(struct unur_gen *gen, double *ui, double *zi, double *xval)
{
  double f0,f1,dq;
  if (_unur_iszero(ui[2])) return TRUE;
  dq = (xval[2] - xval[0]) / ui[2];
  f0 = zi[0];
  f1 = 1./PDF(xval[2]);
  return (f0 <= 3. * dq && f1 <= 3. * dq) ? TRUE : FALSE;
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
