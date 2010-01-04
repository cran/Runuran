/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#define MAX_FLAT_COUNT  (40)
double
_unur_ninv_newton( const struct unur_gen *gen, double U )
{ 
  double x;           
  double fx;          
  double dfx;         
  double fxabs;       
  double xtmp, fxtmp; 
  double xold;        
  double fxtmpabs;    
  double damp;        
  double step;        
  int i;              
  int flat_count;     
  double rel_u_resolution; 
  int x_goal, u_goal; 
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);
  rel_u_resolution = ( (GEN->u_resolution > 0.) ? 
                       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
                       INFINITY );
  if (GEN->table_on) {
    if ( _unur_FP_same(GEN->CDFmin,GEN->CDFmax) ) {
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (U - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }
    if (_unur_FP_is_infinity(GEN->table[i+1])) {
      x  = GEN->table[i];
      fx = GEN->f_table[i];
    }
    else {
      x  = GEN->table[i+1];
      fx = GEN->f_table[i+1];
    }
  }
  else { 
    x  = GEN->s[0];
    fx = GEN->CDFs[0];
  }
  if ( x < DISTR.trunc[0] ){
    x  = DISTR.trunc[0];
    fx = GEN->Umin;    
  }
  else if ( x > DISTR.trunc[1] ){
    x  = DISTR.trunc[1];
    fx = GEN->Umax;    
  }
  fx   -= U;
  dfx   = PDF(x);
  fxabs = fabs(fx);
  xold  = x;    
  damp = 2.;          
  step = 1.;
  for (i=0; i < GEN->max_iter; i++) {
    flat_count = 0;
    while (_unur_iszero(dfx)) {   
      if (_unur_iszero(fx))  
	break; 
      if (fx > 0.) {         
        xtmp = x - step; 
	xtmp = _unur_max( xtmp, DISTR.domain[0] );
      }
      else {
        xtmp  = x + step;
	xtmp = _unur_min( xtmp, DISTR.domain[1] );
      }
      fxtmp    = CDF(xtmp) - U;
      fxtmpabs = fabs(fxtmp);
      if ( fxtmpabs < fxabs ) {       
        step = 1.;     
        x     = xtmp;
        fx    = fxtmp;
      }
      else if ( fxtmp*fx < 0. ) {     
        step /= 2.;                      
      } 
      else {                          
        step *= 2.;    
        x     = xtmp;
        fx    = fxtmp;
      }  
      dfx   = PDF(x);
      fxabs = fabs(fx);
      if (flat_count < MAX_FLAT_COUNT)
	flat_count++;
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
		    "Newton's method cannot leave flat region");
	x = _unur_max( x, DISTR.trunc[0]);
	x = _unur_min( x, DISTR.trunc[1]);
	return x;
      }
    }   
    step = 1.;   
    if (_unur_iszero(fx))  
      break;
    if (_unur_isfinite(dfx)) {
      do {    
	damp /= 2.;
	xtmp = x - damp * fx/dfx;
	xtmp = _unur_min( xtmp, DISTR.trunc[1] );
	xtmp = _unur_max( xtmp, DISTR.trunc[0] );
	fxtmp = CDF(xtmp) - U;
      } while (fabs(fxtmp) > fxabs * (1.+UNUR_SQRT_DBL_EPSILON));   
    }
    else {
      xtmp = 0.5*(x + xold);
      fxtmp = CDF(xtmp) - U;
    }
    damp  = 2.;       
    xold  = x;        
    x     = xtmp;     
    fx    = fxtmp;    
    dfx   = PDF(x);   
    fxabs = fabs(fx); 
    if ( GEN->x_resolution > 0. ) {
      if ( _unur_iszero(fx) ||                             
           fabs(x-xold) < GEN->x_resolution * (fabs(x) + GEN->x_resolution) ) {
 	x_goal = TRUE;
      }
      else
        x_goal = FALSE;
    }
    else {
      x_goal = TRUE;
    }
    if ( GEN->u_resolution > 0. ) {
      if ( fabs(fx) < 0.9 * rel_u_resolution ) {    
      	u_goal = TRUE;
      }
      else if ( _unur_FP_same(xold, x) ) {
        _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
                      "sharp peak or pole: accuracy goal in u cannot be reached");
        u_goal = TRUE;
      }
      else
        u_goal = FALSE;
    }
    else {
      u_goal = TRUE;
    }
    if (x_goal && u_goal)
      break;
  }  
  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");
  x = _unur_max( x, DISTR.trunc[0]);
  x = _unur_min( x, DISTR.trunc[1]);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample(gen, U, x, fx, i);
#endif
  return x;
} 
#undef MAX_FLAT_COUNT
