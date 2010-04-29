/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distributions.c                                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2009 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "Runuran.h"

/*---------------------------------------------------------------------------*/

#define DISTRIBUTION(name)  do { \
    if ( strcmp( distribution, #name)==0 ) \
      return unur_distr_##name (params,n_params); \
    } while(0)


/*****************************************************************************/
/* Continuous distribution                                                   */

UNUR_DISTR *
_Runuran_get_std_cont( const char *distribution, const double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Create UNU.RAN object for special continuous distribution.           */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   distribution ... name of special distribution                      */
     /*   params       ... vector of parameter values                        */
     /*   n_params     ... number of parameters                              */
     /*----------------------------------------------------------------------*/
{
  switch (distribution[0]) {

  case 'b':
    DISTRIBUTION (beta);
    break;
    
  case 'c':
    DISTRIBUTION (cauchy);
    DISTRIBUTION (chi);
    DISTRIBUTION (chisquare);
    break;

  case 'e':
    DISTRIBUTION (exponential);
    DISTRIBUTION (extremeI);
    DISTRIBUTION (extremeII);
    break;

  case 'F':
    DISTRIBUTION (F);
    break;
      
  case 'g':
    DISTRIBUTION (gamma);
    DISTRIBUTION (ghyp);
    DISTRIBUTION (gig);
    DISTRIBUTION (gig2);
    break;

  case 'h':
    DISTRIBUTION (hyperbolic);
    break;
      
  case 'i':
    DISTRIBUTION (ig);
    break;

  case 'l':
    DISTRIBUTION (laplace);
    DISTRIBUTION (logistic);
    DISTRIBUTION (lomax);
    DISTRIBUTION (lognormal);
    break;

  case 'n':
    DISTRIBUTION (normal);
    break;

  case 'p':
    DISTRIBUTION (pareto);
    DISTRIBUTION (powerexponential);
    break;

  case 'r':
    DISTRIBUTION (rayleigh);
    break;

  case 's':
    DISTRIBUTION (slash);
    DISTRIBUTION (student);
    break;

  case 't':
    DISTRIBUTION (triangular);
    break;

  case 'w':
    DISTRIBUTION (weibull);
    break;
  }

  /* unknown distribution */
  return NULL;
} /* _Runuran_get_std_cont() */
 

/*****************************************************************************/
/* Discrete distribution                                                     */

UNUR_DISTR *
_Runuran_get_std_discr( const char *distribution, const double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* Create UNU.RAN object for special discrete distribution.             */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   distribution ... name of special distribution                      */
     /*   params   ... vector of parameter values                            */
     /*   n_params ... number of parameters                                  */
     /*----------------------------------------------------------------------*/
{
  switch (distribution[0]) {
  case 'b':
    DISTRIBUTION (binomial);
    break;

  case 'g':
    DISTRIBUTION (geometric);
    break;

  case 'h':
    DISTRIBUTION (hypergeometric);
    break;

  case 'l':
    DISTRIBUTION (logarithmic);
    break;

  case 'n':
    DISTRIBUTION (negativebinomial);
    break;

  case 'p':
    DISTRIBUTION (poisson);
    break;
  }

  /* unknown distribution */
  return NULL;
} /* _Runuran_get_std_discr() */
 
/*---------------------------------------------------------------------------*/
