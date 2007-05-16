#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran                                                         ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran", 
         ## slots:
         representation( 
                        unur       = "externalptr",   # pointer to UNU.RAN object
                        distr      = "unuran.distr",   # pointer to S4 distribution object
                        distr.str  = "character",     # distribution
                        method.str = "character"      # generation method
                        ),
         ## defaults for slots
         prototype = list(
                 unur       = NULL,
                 distr      = NULL,
                 distr.str  = character(),
                 method.str = "auto"
                 ),
         ## misc
         sealed = TRUE
         )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran",  
          function(.Object, distr=NULL, method="auto") {

                  ## Check entries
                  if (is.null(distr)) {
                          stop("no distribution given", call.=FALSE) }
                  if (!is.character(method)) {
                          stop("'method' must be a character string", call.=FALSE) }

                  ## Store informations 
                  .Object@distr.str <- ifelse(is.character(distr), distr, "[S4 class]")
                  .Object@method.str <- method

                  ## Create UNU.RAN object
                  if (is.character(distr)) {
                          .Object@unur <-.Call("Runuran_init", .Object, distr, method, PACKAGE="Runuran")
                  } else { if (class(distr)=="unuran.discr" ||
                               class(distr)=="unuran.cont"  ||
                               class(distr)=="unuran.cmv") {
                          .Object@unur <-.Call("Runuran_init", .Object, distr@distr, method, PACKAGE="Runuran")
                          .Object@distr <- distr
                  } else {
                          stop("'distr' must be a character string or a Runuran distribution object", call.=FALSE)
                  }}

                  ## Check UNU.RAN object
                  if (is.null(.Object@unur)) {
                          stop("Cannot create UNU.RAN object", call.=FALSE)
                  }
                  
                  ## return new UNU.RAN object
                  .Object
          } )

## Shortcut
unuran.new <- function(distr,method="auto") {
        new("unuran",distr,method)
}

## Validity -----------------------------------------------------------------

## Sampling -----------------------------------------------------------------

## unuran.sample
## ( We avoid using a method as this has an expensive overhead. )
unuran.sample <- function(unr,n=1) { 
        .Call("Runuran_sample", unr, n, PACKAGE="Runuran")
}

## r
##    method alias for unuran.sample  (slow!!)
if(!isGeneric("r"))
        setGeneric("r", function(unur,...) standardGeneric("r"))

setMethod("r", "unuran",
          function(unur,n=1) {
                  .Call("Runuran_sample", unur, n, PACKAGE="Runuran")
          } )

## Printing -----------------------------------------------------------------

## print strings of UNU.RAN object
setMethod( "print", "unuran",
          function(x, ...) {
                  cat("\nObject is UNU.RAN object:\n")
                  cat("\tmethod: ",x@method.str,"\n")
                  cat("\tdistr:  ",x@distr.str,"\n\n")
} )

setMethod( "show", "unuran",
          function(object) { print(object) } )


## End ----------------------------------------------------------------------
