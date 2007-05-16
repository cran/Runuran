#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.cont                                                    ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.cont", 
         ## add slots for continuous univariate distributions
         representation = representation(
                 cdf  = "function",    # CDF of distribution
                 pdf  = "function",    # PDF of distribution
                 dpdf = "function"     # derivative of PDF of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 cdf  = NULL,
                 pdf  = NULL,
                 dpdf = NULL
                 ),
         ## superclass
         contains = "unuran.distr",
         ## seal this class
         sealed = TRUE )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.cont",
          function(.Object, cdf=NULL, pdf=NULL, dpdf=NULL, islog=TRUE, lb=-Inf, ub=Inf) {
                  ## cdf .... cumulative distribution function (CDF)
                  ## pdf .... probability density function (PDF)
                  ## dpdf ... derivative of PDF
                  ## islog .. whether CDF and PDF are given as logarithms
                  ##          (the dpdf is then the derative of log(pdf)!)
                  ## lb ..... lower bound of domain
                  ## ub ..... upper bound of domain

                  ## Check entries
                  if(! (is.numeric(lb) && is.numeric(ub) && lb < ub) )
                          stop("invalid domain ('lb','ub')", call.=FALSE)

                  if(! (is.function(cdf) || is.null(cdf)) )
                          stop("invalid argument 'cdf'", call.=FALSE)
                  if(! (is.function(pdf) || is.null(pdf)) )
                          stop("invalid argument 'pdf'", call.=FALSE)
                  if(! (is.function(dpdf) || is.null(dpdf)) )
                          stop("invalid argument 'dpdf'", call.=FALSE)

                  if(!is.logical(islog))
                          stop("argument 'islog' must be boolean", call.=FALSE)

                  ## Store informations (if provided)
                  if (is.function(cdf))  .Object@cdf  <- cdf
                  if (is.function(pdf))  .Object@pdf  <- pdf
                  if (is.function(dpdf)) .Object@dpdf <- dpdf

                  ## We need an evironment for evaluating R expressions
                  .Object@env <- new.env()
                  
                  ## Create UNUR_DISTR object
                  .Object@distr <-.Call("Runuran_cont_init",
                                        .Object, .Object@env,
                                        .Object@cdf, .Object@pdf, .Object@dpdf, islog,
                                        c(lb,ub),
                                        PACKAGE="Runuran")

                  ## Check UNU.RAN object
                  if (is.null(.Object@distr)) {
                          stop("Cannot create UNU.RAN distribution object", call.=FALSE)
                  }

                  ## return new UNU.RAN object
                  .Object
          } )

## End ----------------------------------------------------------------------
