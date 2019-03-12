## --------------------------------------------------------------------------
##
## miscellaneous auxiliary functions
##
## --------------------------------------------------------------------------

## --- error and warning ----------------------------------------------------
## Remark: same as .rvgt.stop() from package 'rvgt'

.Runuran.stop <- function(..., calledby=NULL) {
    if (is.null(calledby)) {
        calledby <- 
            match.call(sys.function(sys.parent(1L)),
                       sys.call(sys.parent(1L)),
                       TRUE, parent.frame(2L))
    }

    switch(class(calledby),
           "call" = 
               stop(deparse(calledby),":\n   ",..., call.=FALSE),
           
           "character" = 
               stop(calledby,":\n   ",..., call.=FALSE),
           
           stop("Internal error"))
}

##.Runuran.warning <- function(..., calledby=NULL) {
##    if (is.null(calledby)) {
##        calledby <- 
##            match.call(sys.function(sys.parent(1L)),
##                       sys.call(sys.parent(1L)),
##                       TRUE, parent.frame(2L))
##    }
##
##    switch(class(calledby),
##           "call" = 
##               warning(deparse(calledby),":\n   ",..., call.=FALSE),
##           
##           "character" = 
##               warning(calledby,":\n   ",..., call.=FALSE),
##           
##           stop("Internal error"))
##}

##.Runuran.message <- function(..., calledby=NULL) {
##    if (is.null(calledby)) {
##        calledby <- 
##            match.call(sys.function(sys.parent(1L)),
##                       sys.call(sys.parent(1L)),
##                       TRUE, parent.frame(2L))
##    }
##
##    switch(class(calledby),
##           "call" = 
##               message(deparse(calledby),":\n   ",...),
##           
##           "character" = 
##               message(calledby,":\n   ",...),
##           
##           stop("Internal error"))
##}

## --- End ------------------------------------------------------------------
