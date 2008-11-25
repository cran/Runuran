#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.discr                                                   ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.discr", 
         ## add slots for discrete distributions
         representation = representation(
                 pmf  = "function"     # PMF of distribution
                 ),
         ## defaults for slots
         prototype = list(
                 pmf  = NULL
                 ),
         ## superclass
         contains = "unuran.distr",
         ## seal this class
         sealed = TRUE )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.discr",
          function(.Object, pv=NULL, pmf=NULL, mode=NA, lb=0, ub=Inf, sum=NA, name=NA) {
                  ## pv ..... probability vector (PV)
                  ## pmf .... probability mass function (PMF)
                  ## mode ... mode of distribution
                  ## lb ..... lower bound of domain
                  ## ub ..... upper bound of domain
                  ## sum .... sum over PV / PMF
                  ## name ... name of distribution
                  
                  ## Check entries
                  if(! (is.numeric(lb) && is.numeric(ub) && lb < ub) )
                          stop("invalid domain ('lb','ub')", call.=FALSE)
                  if (!(is.numeric(pv) || is.null(pv) ))
                          stop("invalid argument 'pv'", call.=FALSE)
                  if(! (is.function(pmf) || is.null(pmf)) )
                          stop("invalid argument 'pmf'", call.=FALSE)

                  ## Store informations (if provided)
                  if (!is.na(name))     .Object@name <- name
                  if (is.function(pmf)) .Object@pmf  <- pmf
                  ## (There is no need to store the PV)

                  ## We need an evironment for evaluating R expressions
                  .Object@env <- new.env()
                  
                  ## Create UNUR_DISTR object
                  .Object@distr <-.Call("Runuran_discr_init",
                                        .Object, .Object@env,
                                        pv, .Object@pmf,
                                        mode, c(lb,ub), sum, name,
                                        PACKAGE="Runuran")

                  ## Check UNU.RAN object
                  if (is.null(.Object@distr)) {
                          stop("Cannot create UNU.RAN distribution object", call.=FALSE)
                  }

                  ## return new UNU.RAN object
                  .Object
          } )

## Shortcut
unuran.discr.new <- function(pv=NULL, pmf=NULL, mode=NA, lb=0, ub=Inf, sum=NA, name=NA) {
        new("unuran.discr", pv=pv, pmf=pmf, mode=mode, lb=lb, ub=ub, sum=sum,name=name)
}

## End ----------------------------------------------------------------------
