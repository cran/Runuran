#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Wrapper for UNU.RAN sampling methods                                  ##
##                                                                         ##
##   Interface to the UNU.RAN library for                                  ##
##   Universal Non-Uniform RANdom variate generators                       ##
##                                                                         ##
#############################################################################

#############################################################################
##                                                                          #
## Auxiliary functions                                                      #
##                                                                          #
#############################################################################

## -- Numerical derivative --------------------------------------------------

numerical.derivative <- function (x, func, lb=-Inf, ub=Inf, xmin=1, delta=1.e-7) {
        ## x      ... argument
        ## func   ... function for which derivative has to bee computed
        ## lb, ub ... domain of 'func' (not implemented yet)
        ## xmin   ... (not implemented yet)
        ## delta  ... delta for computing difference
        h <- pmax(x*delta,delta)
        df <- (func(x+h)-func(x-h))/(2*h)
}


#############################################################################
##                                                                          #
## Sampling methods for continuous univariate Distributions                 #
##                                                                          #
#############################################################################

## -- ARS: Adaptive Rejection Sampling (TDR with T=log) ---------------------
##
## Type: Rejection
##
## Generate continuous random variates from a given PDF
##

ars.new <- function (logpdf, dlogpdf=NULL, lb=-Inf, ub=Inf, ...) {

        ## check arguments
        if (missing(logpdf) || !is.function(logpdf))
                stop ("argument 'logpdf' missing or invalid")

        ## internal version of logPDF
        f <- function(x) logpdf(x, ...) 
        
        ## derivative of the PDF
        if (is.null(dlogpdf)) {
                ## use numerical derivative
                df <- function(x) {
                        numerical.derivative(x,f)
                }
        }
        else {
                if (! is.function(dlogpdf) )
                        stop ("argument 'dlogpdf' invalid")
                else df <- function(x) dlogpdf(x,...)
	}	

        ## S4 class for continuous distribution
        dist <- new("unuran.cont", pdf=f, dpdf=df, lb=lb, ub=ub, islog=TRUE)

        ## create and return UNU.RAN object
        unuran.new(dist, "ars")
} 


## -- ITDR: Inverse Transformed Density Rejection ---------------------------
##
## Type: Rejection
##
## Generate continuous random variates from a given PDF
##

itdr.new <- function (pdf, dpdf, pole, lb=0, ub=Inf, islog=FALSE, ...) {

        ## check arguments 
        if (missing(pdf) || !is.function(pdf))
                stop ("argument 'pdf' missing or invalid")
        if (missing(dpdf) || !is.function(dpdf))
                stop ("argument 'dpdf' missing or invalid")
        if (missing(pole) || !is.numeric(pole)) 
                stop ("argument 'pole' missing or invalid")

        ## internal versions of PDF and its derivative
        f <- function(x) pdf(x, ...) 
        df <- function(x) dpdf(x, ...) 
        
        ## S4 class for continuous distribution
        dist <- new("unuran.cont", pdf=f, dpdf=df, lb=lb, ub=ub, islog=islog, mode=pole)

        ## create and return UNU.RAN object
        unuran.new(dist, "itdr")
}


## -- PINV:  Polynomial interpolation based INVersion -----------------------
##
## Type: Inversion
##
## Generate continuous random variates from a given PDF or CDF
##

pinv.new <- function (pdf, cdf, lb=-Inf, ub=Inf, islog=FALSE, center=0, uresolution=1.e-10, ...) {

        ## check arguments
        if (missing(pdf) && missing(cdf))
                stop ("argument 'pdf' or 'cdf' required")
        if (!is.numeric(center))
                stop ("argument 'center' invalid")
        if (!is.numeric(uresolution))
                stop ("argument 'uresolution' invalid")

        ## use PDF or CDF ?
        usefunc <- if (missing(pdf)) "usecdf;" else "usepdf;"
                
        ## create internal version of PDF and CDF
        if (!missing(pdf)) {
                ## PDF given
                if (!is.function(pdf))
                        stop ("argument 'pdf' must be of class 'function'")
                PDF <- function(x) pdf(x, ...)
        }
        else {
                PDF <- NULL;
        }

        if (!missing(cdf)) {
                ## CDF given
                if (!is.function(cdf))
                        stop ("argument 'cdf' must be of class 'function'")
                CDF <- function(x) cdf(x, ...)
        }
        else {
                CDF <- NULL
        }

        ## S4 class for continuous distribution
        dist <- new("unuran.cont", pdf=PDF, cdf=CDF, lb=lb, ub=ub, center=center, islog=islog)

        ## create and return UNU.RAN object
        method <- paste("pinv;",usefunc,"u_resolution=",uresolution, sep="")
        return (unuran.new(dist,method))
}


## -- SROU: Simple Ratio-Of-Uniforms Method ---------------------------------
##
## Type: Rejection
##
## Generate continuous random variates from a given PDF
##

srou.new <- function (pdf, mode, area, lb=-Inf, ub=Inf, islog=FALSE, r=1, ...) {

        ## check arguments
        if (missing(pdf) || !is.function(pdf))
                stop ("argument 'pdf' missing or invalid")
        if (missing(mode) || !is.numeric(mode))
                stop ("argument 'mode' missing or invalid")
        if (missing(area) || !is.numeric(area))
                stop ("argument 'area' missing or invalid")

        ## internal version of PDF
        f <- function(x) pdf(x, ...) 

        ## S4 class for continuous distribution
        dist <- new("unuran.cont", pdf=f, lb=lb, ub=ub, islog=islog, mode=mode, area=area)

        ## create and return UNU.RAN object
        method <- paste("srou; r=",r, sep="")
        unuran.new(dist, method)
}


## -- TDR: Transformed Density Rejection ------------------------------------
##
## Type: Rejection
##
## Generate continuous random variates from a given PDF
##

tdr.new <- function (pdf, dpdf=NULL, lb=-Inf, ub=Inf, islog=FALSE, ...) {

        ## check arguments
        if (missing(pdf) || !is.function(pdf))
                stop ("argument 'pdf' missing or invalid")

        ## internal version of PDF
        f <- function(x) pdf(x, ...) 

        ## we also need the derivative of the PDF
        if (is.null(dpdf)) {
                ## use numerical derivative
                df <- function(x) {
                        numerical.derivative(x,f)
                }
        }
        else {
                if (! is.function(dpdf) )
                        stop ("argument 'dpdf' invalid")
                else df <- function(x) dpdf(x,...)
	}	

        ## S4 class for continuous distribution
        dist <- new("unuran.cont", pdf=f, dpdf=df, lb=lb, ub=ub, islog=islog)

        ## create and return UNU.RAN object
        unuran.new(dist, "tdr")
}


#############################################################################
##                                                                          #
## Sampling methods for discrete univariate Distributions                   #
##                                                                          #
#############################################################################

## -- DARI: Discrete Automatic Rejection Inversion --------------------------
##
## Type: Rejection
##
## Generate discrete random variates from a given probability vector
## using Discrete Automatic Rejection Inversion.
##

dari.new <- function (pmf, lb=0, ub=Inf, mode=NA, sum=1, ...) {

        ## check arguemngts
        if (missing(pmf) || !is.function(pmf))
                stop ("argument 'pmf' missing or invalid")

        ## internal version of PMF
        f <- function(x) pmf(x, ...)
        
        ## S4 class for discrete distribution
        distr <- new("unuran.discr",pmf=f,lb=lb,ub=ub,mode=mode,sum=sum)

        ## create and return UNU.RAN object
        new("unuran", distr, "DARI")
}


## -- DAU: Alias-Urn Method ------------------------------------------------
##
## Type: Patchwork
##
## Generate discrete random variates from a given probability vector
## using the Alias-Urn Method
##
## Remark: we do not pass the domain to UNU.RAN
##

dau.new <- function (pv, from=1) {

        ## check arguments
        if (missing(pv) || !is.numeric(pv))
                stop ("argument 'pv' missing or invalid")

        ## S4 class for discrete distribution
        distr <- new("unuran.discr",pv=pv,lb=from)

        ## create and return UNU.RAN object
        new("unuran", distr, "DAU")
}

## -- DGT: Guide Table Method -----------------------------------------------
##
## Type: Inversion
##
## Generate discrete random variates from a given probability vector
## using the Guide-Table Method for discrete inversion
##
## Remark: we do not pass the domain to UNU.RAN
##

dgt.new <- function (pv, from=1) {

        ## check arguments
        if (missing(pv) || !is.numeric(pv))
                stop ("argument 'pv' missing or invalid")

        ## S4 class for discrete distribution
        distr <- new("unuran.discr",pv=pv,lb=from)

        ## create and return UNU.RAN object
        new("unuran", distr, "DGT")
}


#############################################################################
##                                                                          #
## Sampling methods for continuous multivariate Distributions               #
##                                                                          #
#############################################################################


## -- HITRO: Hit-and-Run sampler with Ratio-of-Uniforms ---------------------
##
## Type: MCMC
##
## Generate continuous random variates from a given PDF
##

hitro.new <- function (dim=1, pdf, ll=NULL, ur=NULL, mode=NULL, center=NULL, thinning=1, burnin=0, ...) {

        ## check arguments
        if (missing(pdf) || !is.function(pdf))
                stop ("argument 'pdf' missing or invalid")

        ## internal version of PDF
        f <- function(x) pdf(x, ...) 

        ## S4 class for continuous multivariate distribution
        dist <- new("unuran.cmv", dim=dim, pdf=f, mode=mode, center=center, ll=ll, ur=ur)

        ## create and return UNU.RAN object
        method <- paste("hitro;thinning=",thinning,";burnin=",burnin, sep="")
        unuran.new(dist, method)
}

## -- VNROU: Multivariate Naive Ratio-Of-Uniforms method --------------------
##
## Type: Rejection / Ratio-of-Uniforms
##
## Generate continuous random variates from a given PDF
##

vnrou.new <- function (dim=1, pdf, ll=NULL, ur=NULL, mode=NULL, center=NULL, ...) {

        ## check arguments
        if (missing(pdf) || !is.function(pdf))
                stop ("argument 'pdf' missing or invalid")

        ## internal version of PDF
        f <- function(x) pdf(x, ...) 

        ## S4 class for continuous multivariate distribution
        dist <- new("unuran.cmv", dim=dim, pdf=f, mode=mode, center, ll=ll, ur=ur)

        ## create and return UNU.RAN object
        unuran.new(dist, "VNROU")
}

## -- End -------------------------------------------------------------------
