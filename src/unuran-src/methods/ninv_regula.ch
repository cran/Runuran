/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#define MAX_STEPS (100)
#define STEPFAC  (0.4)
#define I_CHANGE_TO_BISEC (50)
double 
_unur_ninv_regula( const struct unur_gen *gen, double u )
{ 
  double x1, x2, a, xtmp;
  double f1, f2,fa, ftmp;
  double length;         
  double lengthabs;      
  double lengthsgn;      
  double dx;             
  int count_nosc = 0;    
  int i;                 
  double min_step_size;  
  double rel_u_resolution; 
  CHECK_NULL(gen, INFINITY);  COOKIE_CHECK(gen, CK_NINV_GEN, INFINITY);
  rel_u_resolution = ( (GEN->u_resolution > 0.) ? 
		       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
		       INFINITY );
  if ( _unur_ninv_bracket( gen, u, &x1, &f1, &x2, &f2 ) 
       != UNUR_SUCCESS )
    return x2;
  a = x1; fa = f1; 
  for (i=0; TRUE; i++) {
    if ( f1*f2 < 0.) { 
      count_nosc = 0;   
      if ( fabs(f1) < fabs(f2) ) {
	xtmp = x1; ftmp = f1;
	x1 = x2;   f1 = f2;
	x2 = xtmp; f2 = ftmp;
      }
      a = x1; fa = f1;
    }
    else {
      count_nosc++;
    }
    length = x2 - a;                       
    lengthabs = fabs(length);              
    lengthsgn = (length < 0.) ? -1. : 1.;
    if (_unur_ninv_accuracy( gen, GEN->x_resolution, rel_u_resolution,
    			     x2, f2, a, fa ))
      break;
    if (i >= GEN->max_iter)
      break;
    dx = (_unur_FP_same(f1,f2)) ? length/2. : f2*(x2-x1)/(f2-f1) ;  
    if (GEN->u_resolution < 0.) 
      min_step_size = fabs(x2) * GEN->x_resolution;
    else
      min_step_size = lengthabs * DBL_EPSILON; 
    if ( fabs(dx) < min_step_size ) {
      dx = lengthsgn * 0.99 * min_step_size;
      while ( x2 == x2 - dx ){ 
        if ( dx != 2.*dx)    
          dx = 2.*dx;
        else
          dx = length/2.;    
      }
    }
    if ( count_nosc > 1 || i > I_CHANGE_TO_BISEC ||
	 (lengthabs-GEN->x_resolution*fabs(x2))/(dx*lengthsgn) <= 1. )
      dx = length/2.; 
    x1 = x2;       f1 = f2;
    x2 = x2-dx;    f2 = CDF(x2) - u; 
  }  
  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");
  x2 = _unur_max( x2, DISTR.trunc[0]);
  x2 = _unur_min( x2, DISTR.trunc[1]);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample( gen,u,x2,f2,i );
#endif
  return x2;
} 
double 
_unur_ninv_bisect( const struct unur_gen *gen, double u )
{ 
  double x1, x2, mid=0.;   
  double f1, f2, fmid;     
  int i;                   
  double rel_u_resolution; 
  CHECK_NULL(gen, INFINITY);  COOKIE_CHECK(gen, CK_NINV_GEN, INFINITY);
  rel_u_resolution = ( (GEN->u_resolution > 0.) ?
		       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
		       INFINITY );
  if ( _unur_ninv_bracket( gen, u, &x1, &f1, &x2, &f2 )
       != UNUR_SUCCESS )
    return x2;
  for (i=0; i<GEN->max_iter; i++) {
    mid = x1 + (x2-x1)/2.;
    fmid = CDF(mid) - u;
    if (f1*fmid <= 0) {
      x2 = mid; f2 = fmid;
      if (_unur_ninv_accuracy( gen, GEN->x_resolution, rel_u_resolution,
			       mid, fmid, x1, f1 ))
	break;
    }
    else {
      x1 = mid; f1 = fmid;
      if (_unur_ninv_accuracy( gen, GEN->x_resolution, rel_u_resolution,
			       mid, fmid, x2, f2 ))
	break;
    }
  }
  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");
  mid = _unur_max( mid, DISTR.trunc[0]);
  mid = _unur_min( mid, DISTR.trunc[1]);
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample( gen,u,mid,CDF(mid)-u,i );
#endif
  return mid;
} 
int
_unur_ninv_bracket( const struct unur_gen *gen, double u, 
		    double *xl, double *fl, double *xu, double *fu )
{ 
  int i;                 
  double x1, x2, xtmp;   
  double f1, f2, ftmp;   
  double step;           
  int step_count;        
  if (GEN->table_on) {
    if ( _unur_FP_same(GEN->CDFmin, GEN->CDFmax) ) {
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (u - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }
    if ( ! _unur_FP_is_minus_infinity(GEN->table[i]) ){
      x1 = GEN->table[i];
      f1 = GEN->f_table[i]; 
    }
    else{
      x1 = GEN->table[i+1] + (GEN->table[i+1] - GEN->table[i+2]);
      f1 = CDF(x1);
    }
    if( ! _unur_FP_is_infinity(GEN->table[i+1]) ){
      x2 = GEN->table[i+1];
      f2 = GEN->f_table[i+1];
    }
    else{
      x2 = GEN->table[i] + (GEN->table[i] - GEN->table[i-1]);
      f2 = CDF(x2);
    }
  }
  else { 
    x1 =  GEN->s[0];      
    f1 =  GEN->CDFs[0];
    x2 =  GEN->s[1];         
    f2 =  GEN->CDFs[1];
  }
  if ( x1 >= x2 ) { 
    xtmp = x1; ftmp = f1;
    x1   = x2; f1   = f2;
    x2 = xtmp + fabs(xtmp)*DBL_EPSILON;
    f2 = CDF(x2); 
  }
  if ( x1 < DISTR.trunc[0] || x1 >= DISTR.trunc[1] ){
    x1 = DISTR.trunc[0];
    f1 = GEN->Umin;    
  }
  if ( x2 > DISTR.trunc[1] || x2 <= DISTR.trunc[0] ){
    x2 = DISTR.trunc[1];
    f2 = GEN->Umax;    
  }
  f1 -= u;  f2 -= u;
  step = (GEN->s[1]-GEN->s[0]) * STEPFAC;
  step_count = 0;
  while ( f1*f2 > 0. ) {
    if ( f1 > 0. ) {     
      x2  = x1;  
      f2  = f1;
      x1 -= step;   
      f1  = CDF(x1) - u;
    }
    else {         
      x1  = x2;
      f1  = f2;
      x2 += step;
      f2  = CDF(x2) - u;
    }
    if (step_count < MAX_STEPS) {
      ++step_count;
      step *= 2.;
      if( step_count > 20 && step < 1.) step = 1.; 
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
                  "Regula Falsi cannot find interval with sign change");
      *xu = (f1>0.) ? DISTR.trunc[0] : DISTR.trunc[1];
      return UNUR_ERR_GEN_SAMPLING;
    }
  }
  *xl = x1; *xu = x2;
  *fl = f1; *fu = f2;
  return UNUR_SUCCESS;
} 
int
_unur_ninv_accuracy( const struct unur_gen *gen,
		     double x_resol, double u_resol,
		     double x0, double f0, double x1, double f1 )
{ 
  int x_goal, u_goal;     
  if ( x_resol > 0. ) {
    if ( _unur_iszero(f0) ||          
  	 fabs(x1-x0) < x_resol * (fabs(x0) + x_resol) ) {
      x_goal = TRUE;
    }
    else if ( _unur_FP_same(f0,f1) ) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
  		    "flat region: accuracy goal in x cannot be reached");
      x_goal = TRUE;
    }
    else
      x_goal = FALSE;
  }
  else {
    x_goal = TRUE;
  }
  if ( GEN->u_resolution > 0. ) {
    if (fabs(f0) < 0.9 * u_resol) {
      u_goal = TRUE;
    }
    else if ( _unur_FP_same(x0,x1) ) {
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
  return (x_goal && u_goal);
} 
#undef MAX_STEPS
#undef STEPFAC
#undef I_CHANGE_TO_BISEC
