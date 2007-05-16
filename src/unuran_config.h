/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      unuran_config.h                                              *
 *                                                                           *
 *   compiler switches, compile time options and default values              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
#ifndef UNURAN_CONFIG_H_SEEN
#define UNURAN_CONFIG_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *  Logging and debugging.                                                   *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Default name of log file.                                                 */
/*                                                                           */

#define UNUR_LOG_FILE "unuran.log"

/*---------------------------------------------------------------------------*/
/* Set default flag for debugging of generators:                             */
/*                                                                           */
/*   UNUR_DEBUG_OFF    ... switch off debugging information                  */
/*   UNUR_DEBUG_INIT   ... pameters and structure of generator only          */
/*   UNUR_DEBUG_SETUP  ... information for setup step                        */
/*   UNUR_DEBUG_ADAPT  ... trace adaptive steps                              */ 
/*   UNUR_DEBUG_SAMPLE ... trace sampling                                    */
/*   UNUR_DEBUG_ALL    ... write all available debugging information         */
/*                                                                           */
/* Detailed discription of possible flags in file `./utils/debug.h'          */
/*                                                                           */
/* Debugging information is written into the log file.                       */
/* It only works if the library is built with the ----enable-logging         */
/* configure flag.                                                           */

#define UNUR_DEBUGFLAG_DEFAULT   UNUR_DEBUG_INIT


/*****************************************************************************
 *  Compile time parameters for generators.                                  *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Distribution objects.                                                     */

/* Maximal number of parameters for the PDF of a distribution.               */
/* (It must be at least 5!)                                                  */
#define UNUR_DISTR_MAXPARAMS  5

/* Maximal size of automatically created probability vectors.                */
#define UNUR_MAX_AUTO_PV    100000


/*****************************************************************************
 *  Interface for uniform random number generators.                          *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Default generators:                                                       */
/* We set the default generators in R_init_Runuran()                         */
#define UNUR_URNG_DEFAULT      (NULL)
#define UNUR_URNG_AUX_DEFAULT  (NULL)


/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/
