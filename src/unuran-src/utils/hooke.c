/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "hooke_source.h"
#define HOOKE_SUBITERMAX 10
static double best_nearby(struct unur_funct_vgeneric faux,
                          double *delta, double *point, 
                          double prevbest, int dim);
static double best_nearby(struct unur_funct_vgeneric faux, 
                          double *delta, double *point, 
                          double prevbest, int dim)
{
	   double	   *z;
	   double	   minf, ftmp;
	   int		   i;
	   z=(double *) malloc( dim*sizeof(double));
	   minf = prevbest;
	   for (i = 0; i < dim; i++)
		   z[i] = point[i];
	   for (i = 0; i < dim; i++) {
		   z[i] = point[i] + delta[i];
		   ftmp = faux.f(z, faux.params);
		   if (ftmp < minf)
			   minf = ftmp;
		   else {
			   delta[i] = 0.0 - delta[i];
			   z[i] = point[i] + delta[i];
			   ftmp = faux.f(z, faux.params);
			   if (ftmp < minf)
				   minf = ftmp;
			   else
				   z[i] = point[i];
		   }
	   }
	   for (i = 0; i < dim; i++)
		   point[i] = z[i];
           free(z);		   
	   return (minf);
} 
int _unur_hooke(struct unur_funct_vgeneric faux, 
           int dim, double *startpt, double *endpt, 
           double rho, double epsilon, long itermax)
{
           double  *delta, *xbefore, *newx;
	   double  newf, fbefore, steplength, tmp;
	   int	   i, keep;
	   int	   iters, isubiters;
	   delta   = (double *) malloc( dim*sizeof(double));
	   xbefore = (double *) malloc( dim*sizeof(double));
	   newx    = (double *) malloc( dim*sizeof(double));
	   for (i = 0; i < dim; i++) {
		   newx[i] = xbefore[i] = startpt[i];
		   delta[i] = fabs(startpt[i] * rho);
		   if (_unur_iszero(delta[i])) delta[i] = rho;
	   }
	   steplength = rho;
	   iters = 0;
	   fbefore = faux.f(newx, faux.params);
	   newf = fbefore;
	   while ((iters < itermax) && (steplength > epsilon)) {
		   iters++;
		   for (i = 0; i < dim; i++) {
			   newx[i] = xbefore[i];
		   }
		   newf = best_nearby(faux, delta, newx, fbefore, dim);
		   keep = 1;
		   isubiters=0;
		   while ((newf < fbefore) && (keep == 1) ) {
			   for (i = 0; i < dim; i++) {
				   if (newx[i] <= xbefore[i])
					   delta[i] = 0.0 - fabs(delta[i]);
				   else
					   delta[i] = fabs(delta[i]);
				   tmp = xbefore[i];
				   xbefore[i] = newx[i];
				   newx[i] = newx[i] + newx[i] - tmp;
			   }
			   fbefore = newf;
			   newf = best_nearby(faux, delta, newx, fbefore, dim);
			   if (newf >= fbefore)
				   break;
			   keep = 0;
			   for (i = 0; i < dim; i++) {
			     if (fabs(newx[i] - xbefore[i]) >  (0.5 * fabs(delta[i]))) {
			       keep = 1; 
			       break;
			     }
			   }
			   if ( isubiters++ >= HOOKE_SUBITERMAX) break;  
		   }
		   if ((steplength >= epsilon) && (newf >= fbefore)) {
			   steplength = steplength * rho;
			   for (i = 0; i < dim; i++) {
				   delta[i] *= rho;
			   }
		   }
	   }
	   for (i = 0; i < dim; i++)
		   endpt[i] = xbefore[i];
 	   free(delta); free(xbefore); free(newx);		   
	   return (iters);
} 
