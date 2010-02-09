#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007-2010, Josef Leydold and Wolfgang Hoermann                    ##
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
                        data       = "list",          # list for packed data
                        dom        = "numeric",       # domain of distribution
                        distr      = "unuran.distr",  # pointer to S4 distribution object
                        distr.str  = "character",     # distribution
                        method.str = "character"      # generation method
                        ),
         ## defaults for slots
         prototype = list(
           unur       = NULL,
           data       = NULL,
           dom        = NULL,
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

## ur
## ( We avoid using a method as this has an expensive overhead. )
ur <- function(unr,n=1) { 
        .Call("Runuran_sample", unr, n, PACKAGE="Runuran")
}

## unuran.sample: deprecated name or ur()
unuran.sample <- function(unr,n=1) { 
        .Call("Runuran_sample", unr, n, PACKAGE="Runuran")
}

## Quantile -----------------------------------------------------------------

## uq
uq <- function(unr,U) { 
        .Call("Runuran_quantile", unr, U, PACKAGE="Runuran")
}

## PDF & PMF -- [EXPERIMENTAL] ----------------------------------------------

## ud
ud <- function(distr,x) {
  if ( ! (is(distr,"unuran.cont") || is(distr,"unuran.discr")) )
    stop("argument 'distr' must be UNU.RAN distribution object")
  .Call("Runuran_PDF", distr@distr, x, PACKAGE="Runuran")
}

## Packing ------------------------------------------------------------------

## We have the following two situtations for
## slots 'unur', 'data', and 'dom':
##
## 1. Runuran object is NOT PACKED:
##    'unur' ... contains pointer to UNU.RAN object
##    'data' ... is set to NULL
##    'dom'  ... is ignored
##
## 2. Runuran object is PACKED:
##    'unur' ... contains NULL pointer '(nil)'
##    'data' ... contains list of coefficients
##    'dom'  ... contains domain of distribution
##
## Otherwise, if 'unur' points to '(nil)' and 'data' is still set to NULL
## then the Runuran object is broken.
## (This happens in particular when an unpacked Runuran object is saved in
## a workspace and restored in a new R session.)
##
## It should not happen that both 'unur' and 'data' contain non-NULL objects.
##
## ..........................................................................

## unuran.packed
##    Some Runuran object can be packed such that all data are stored as 
##    R object (and thus can be copied and saved within R)
if(!isGeneric("unuran.packed"))
  setGeneric("unuran.packed", function(unr) standardGeneric("unuran.packed"))

setMethod("unuran.packed", "unuran", 
          function(unr) {
            if (is.null(unr@data)) return(FALSE)
            return(TRUE)
          } )

if(!isGeneric("unuran.packed<-"))
  setGeneric("unuran.packed<-", function(unr, value) standardGeneric("unuran.packed<-"))

setReplaceMethod("unuran.packed", "unuran", 
                 function(unr, value) {
                   value <- as.logical(value)
                   is.packed <- !is.null(unr@data)

                   if (value && is.packed) {
                     warning("[UNU.RAN - warning] object already PACKED", call.=FALSE)
                     
                   }
                   if (!value && is.packed) {
                     ## we cannot unpack object data
                     stop("[UNU.RAN - error] Cannot unpack 'unuran' object", call.=FALSE)
                   }
                   if (value && !is.packed) {
                     ## pack data
                     .Call("Runuran_pack", unr, PACKAGE="Runuran")
                   }
                   ## otherwise: nothing to do
                   
                   return (unr)
                 } )


## Printing -----------------------------------------------------------------

## print strings of UNU.RAN object
setMethod( "print", "unuran",
          function(x, ...) {
                  cat("\nObject is UNU.RAN object:\n")
                  cat("\tmethod: ",x@method.str,"\n")
                  cat("\tdistr:  ",x@distr.str,"\n\n")
                  cat(.Call("Runuran_print", x, FALSE, PACKAGE="Runuran"))
                  cat("")
} )

setMethod( "show", "unuran",
          function(object) { print(object) } )

## unuran.details
## (print for information and hints)
unuran.details <- function(unr, show=TRUE, return.list=FALSE) {
  if (isTRUE(show)) {
    cat("\nObject is UNU.RAN object:\n")
    cat("\tmethod: ",unr@method.str,"\n")
    cat("\tdistr:  ",unr@distr.str,"\n\n")
    info <- .Call("Runuran_print", unr, TRUE, PACKAGE="Runuran")
    cat(info)
  }
  if (isTRUE(return.list)) {
    data <- .Call("Runuran_performance", unr, PACKAGE="Runuran")
    invisible(data)
  }
}


## Mixture ------------------------------------------------------------------

## UNU.RAN meta method for sampling from a mixture of distributions

mixt.new <- function (prob, comp, inversion=FALSE) {

  ## Check arguments
  if (length(prob) != length(comp))
    stop ("'prob' and 'comp' must have same length")
  if (! is.numeric(prob))
    stop ("invalid argment 'prob'")
  if (! is.list(comp))
    stop ("invalid argment 'comp'")

  if (! is(comp[[1]],"unuran")) {
    if (is(comp[[1]],"unuran.distr")) {
      stop ("invalid argment: 'comp' must be list of generators not of distributions")
    } else {
      stop ("invalid argment: 'comp' must be list of generators")
    }
  }
  
  ## Create empty "unuran" object
  ## Unfortunately we cannot use new(), since then method 'initialize'
  ## is called which does not work at all in for this task.
  ## Thus the following is copied from function new().
  ClassDef <- getClass("unuran", where = topenv(parent.frame()))
  .Object <- .Call("R_do_new_object", ClassDef, PACKAGE = "base")

  ## Store informations 
  .Object@distr.str <- "mixture of distributions"
  .Object@method.str <- "mixt"

  ## Create UNU.RAN object
  .Object@unur <- .Call("Runuran_mixt", .Object, prob, comp, inversion, PACKAGE="Runuran")
  if (is.null(.Object@unur)) {
    stop("Cannot create UNU.RAN object", call.=FALSE)
  }

  ## Return new UNU.RAN object
  .Object
}

## End ----------------------------------------------------------------------
