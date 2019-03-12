## --------------------------------------------------------------------------
##
## Handle options for package "Runuran"
##
## --------------------------------------------------------------------------
##
## This is based on 'igraph.options' in package 'igraph'
## by Gabor Csardi <csardi.gabor@gmail.com>
##
## --------------------------------------------------------------------------

## --- Defaults for options -------------------------------------------------

unuran.error.level.default = "warning"

## --- Current list of options ----------------------------------------------

.Runuran.Options <- list(
    error.level = unuran.error.level.default
)

## --- Callback: display unuran errors --------------------------------------

## default: "warning"
.unuran.error.level <- unuran.error.level.default

.Runuran.options.set.error.level <- function(level, calledby) {

    env <- asNamespace("Runuran")

    if (! is.na(pmatch(level, "default"))) {
        level <- unuran.error.level.default
    }

    if (! is.na(pmatch(level, "all"))) {
        level <- "all"
        .Call(C_Runuran_set_error_level, 3L)
    } else if (! is.na(pmatch(level, "warning"))) {
        level <- "warning"
        .Call(C_Runuran_set_error_level, 2L)
    } else if (! is.na(pmatch(level, "error"))) {
        level <- "error"
        .Call(C_Runuran_set_error_level, 1L)
    } else if (! is.na(pmatch(level, "none"))) {
        level <- "none"
        .Call(C_Runuran_set_error_level, 0L)
    } else {
        .Runuran.stop("Invalid value for option 'error.level'. ",
                      "Possible values: \"default\", \"all\", \"warning\", \"error\", \"none\".",
                      calledby=calledby)
    }

    ## store current level 
    assign(".unuran.error.level", level, env)

    ## return level 
    level
}

## --- List of callback functions for setting option values -----------------

.Runuran.options.callbacks <- list(
    ## whether UNU.RAn warnings and errors should be displayed
    error.level = .Runuran.options.set.error.level
)

## ==========================================================================
##'
##' Set or return options for Runuran library
##' 
## --------------------------------------------------------------------------
##'
##' @description
##'
##' Library \pkg{Runuran} has some parameters which (usually) affect the
##' behavior of its functions. These can be set for the whole session
##' via \code{Runuran.options}.
##' 
## --------------------------------------------------------------------------
##'
##' @details
##'
##' The function provides a tool to control the behavior of library
##' \pkg{Runuran}. A list may be given as the only argument, or any
##' number of arguments may be in the \code{name=value} form.
##' If no arguments are specified then the function returns the
##' current settings of all parameters. 
##' If a single option name is given as character string, then its
##' value is returned (or \code{NULL} if it does not exist).
##' Option \emph{values} may be abbreviated.
##'
##' Currently used parameters in alphabetical order:
##' \describe{
##'   \item{error.level}{
##'     verbosity level of error messages and warnings from the
##'     underlying UNU.RAN library. It has no effect on messages
##'     from the routines in this package.
##'     Warnings are useful for analysing possible problems with the
##'     selected combinations of distribution and method.
##'     However, they can produce quite a lot of output if the
##'     conditions of the method is appropriate for a distribution or
##'     the distribution has properties like very heavy tails or
##'     very high peaks.
##'
##'     The following levels can be set:
##'     \describe{
##'       \item{\code{"default"}:}{
##'         same as \code{"warning"}.
##'       }
##'       \item{\code{"none"}:}{
##'         all error messages and warnings are suppressed.
##'       }
##'       \item{\code{"error"}:}{
##'         only show error messages.
##'       }
##'       \item{\code{"warning"}:}{
##'         show error messages and some of the warnings.
##'       }
##'       \item{\code{"all"}:}{
##'         show all error messages and warnings.
##'       }
##'     }
##'   }
##' }
## 
## --------------------------------------------------------------------------
##'
##' @author Josef Leydold \email{josef.leydold@@wu.ac.at}
##'
## --------------------------------------------------------------------------
##' 
##' @examples
##'
##' ## save current options
##' oldval <- Runuran.options()
##'
##' ## show current options
##' Runuran.options("error.level")
##'
##' ## suppress all UNU.RAN error messages and warnings
##' Runuran.options(error.level="none")
##'
##' ## restore Runuran options
##' Runuran.options(oldval)
##'
## --------------------------------------------------------------------------
##
##  Arguments:
##
##' @param \dots
##'        A list may be given as the only argument, or any number of
##'        arguments may be in the \code{name=value} form, a character
##'        string for the name of a parameter, or no argument at all
##'        may be given.
##'
## --------------------------------------------------------------------------
##'
##' @return
##'
##' \code{Runuran.options} returns a list with the updated values of the
##' parameters. If the argument list is not empty, the returned list
##' is invisible. If a single character string is given, then the
##' value of the corresponding parameter is returned (or \code{NULL}
##' if the parameter is not used).
##'
### --------------------------------------------------------------------------
### @export
### --------------------------------------------------------------------------

Runuran.options <- function(...) {
    ## ----------------------------------------------------------------------
    ## 
    ## ----------------------------------------------------------------------
    ## ... : list of options
    ## ----------------------------------------------------------------------

    ## current list of options
    current <- .Runuran.Options
    
    ## no arguments --> return list of option values
    if (nargs() == 0)
        return(current)

    ## read arguments
    temp <- list(...)

    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list =
                   ## case: options given as list
                   temp <- arg,

               character =
                   ## case: ask for value of option 'arg'
                   if( isTRUE(arg %in% names(.Runuran.Options)))
                       return(.Runuran.Options[arg])
                   else
                       return(NULL),
               .Runuran.stop("Invalid argument ", sQuote(arg),"."))
    }

    ## no options given
    if (length(temp) == 0)
        return(current)

    ## no 'key=value' pair given
    if (is.null(names(temp)))
        .Runuran.stop("Options must be given by name.")

    ## check for invalid option names
    for (cn in names(temp)) {
        if(is.na(match(cn, names(.Runuran.Options))))
            .Runuran.stop("Invalid option ", sQuote(cn), ".")
    }

    ## execute callback function on all options (if available)
    cb <- intersect(names(.Runuran.options.callbacks), names(temp))
    for (cn in cb) {
        temp[[cn]] <- .Runuran.options.callbacks[[cn]](temp[[cn]], match.call())
    }

    ## update option values
    current[names(temp)] <- temp
    assign(".Runuran.Options", current, envir = asNamespace("Runuran"))

    ## return updated list
    invisible(current)
}

## --- End ------------------------------------------------------------------
