#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.cmv                                                     ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.cmv", 
         ## add slots for continuous multivariate distributions
         representation = representation(
                 ndim = "integer",    # dimensions of distribution
                 pdf  = "function"    # PDF of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 ndim = as.integer(1),
                 pdf  = NULL
                 ),
         ## superclass
         contains = "unuran.distr",
         ## seal this class
         sealed = TRUE )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.cmv",
          function(.Object, dim=1, pdf=NULL, mode=NULL) {
                  ## dim .... dimension of distribution
                  ## pdf .... probability density function (PDF)

                  ## Check entries
                  ndim <- as.integer(dim)
                  if (ndim < 1 || ndim > 100000)
                          stop("invalid argument 'dim'", call.=FALSE)
                  if(! (is.function(pdf) || is.null(pdf)) )
                          stop("invalid argument 'pdf'", call.=FALSE)
                  if(! (is.numeric(mode) || is.null(mode)) )
                          stop("invalid argument 'mode'", call.=FALSE)
                  if( (! is.null(mode)) && length(mode)!=ndim ) 
                          stop("argument 'mode' must have length 'dim'", call.=FALSE)
                  
                  ## Store informations (if provided)
                  .Object@ndim <- ndim
                  if (is.function(pdf))  .Object@pdf <- pdf

                  ## We need an evironment for evaluating R expressions
                  .Object@env <- new.env()
                  
                  ## Create UNUR_DISTR object
                  .Object@distr <-.Call("Runuran_cmv_init",
                                        .Object, .Object@env,
                                        .Object@ndim, .Object@pdf, mode,
                                        PACKAGE="Runuran")

                  ## Check UNU.RAN object
                  if (is.null(.Object@distr)) {
                          stop("Cannot create UNU.RAN distribution object", call.=FALSE)
                  }

                  ## return new UNU.RAN object
                  .Object
          } )

## End ----------------------------------------------------------------------
