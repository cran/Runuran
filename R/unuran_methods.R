#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Wrapper for special UNU.RAN sampling methods                          ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

#############################################################################
## Sampling methods for discrete univariate Distributions                   #
#############################################################################

## -- DGT: Guide Table Method -----------------------------------------------
##
## Type: Inversion
##
## Generate discrete random variates from a given probability vector
## useing the Guide-Table Method for discrete inversion
##

urdgt <- function (n, probvector, from = 0, by = 1) {
        distr <- new("unuran.discr",pv=probvector)
        unr <- new("unuran", distr, "DGT")
        if (from==0 && by==1)
                unuran.sample(unr,n)
        else
                from + by * unuran.sample(unr,n)
}

## -- DAU: Alias-Urn Method ------------------------------------------------
##
## Type: Patchwork
##
## Generate discrete random variates from a given probability vector
## useing the Alias-Urn Method
##

urdau <- function (n, probvector, from = 0, by = 1) {
        distr <- new("unuran.discr",pv=probvector)
        unr <- new("unuran", distr, "DAU")
        if (from==0 && by==1)
                unuran.sample(unr,n)
        else
                from + by * unuran.sample(unr,n)
}

## -- End -------------------------------------------------------------------


