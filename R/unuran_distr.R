#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Class: unuran.distr                                                   ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

## Initialize global variables ----------------------------------------------

## Class --------------------------------------------------------------------

setClass( "unuran.discr", 
         ## slots:
         representation( 
                        distr   = "externalptr"    # pointer to UNU.RAN distribution object
                        ),
         ## defaults for slots
         prototype = list(
                 unur   = NULL
                 ),
         ## misc
         sealed = TRUE
         )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "unuran.discr",
          function(.Object, pv=NULL) {
                  ## pv ... probability vector

                  ## Check entries
                  if (is.null(pv)) {
                          stop("no probability vector given", call.=FALSE) }
                  if (!is.double(pv)) {
                          stop("'pv' must be a double vector", call.=FALSE) }

                  ## Store informations 
                  ## (currently nothing to do)

                  ## Create UNUR_DISTR object
                  .Object@distr <-.Call("Runuran_discr_init", pv, PACKAGE="Runuran")

                  ## Check UNU.RAN object
                  if (is.null(.Object@distr)) {
                          stop("Cannot create UNU.RAN distribution object", call.=FALSE)
                  }

                  ## return new UNU.RAN object
                  .Object
          } )

## Validity -----------------------------------------------------------------

## Sampling -----------------------------------------------------------------

## Printing -----------------------------------------------------------------

## print strings of UNU.RAN object
setMethod( "print", "unuran.discr",
          function(x, ...) {
                  cat("\nObject is UNU.RAN distribution object\n\n")
} )

setMethod( "show", "unuran.discr",
          function(object) { print(object) } )

## End ----------------------------------------------------------------------
