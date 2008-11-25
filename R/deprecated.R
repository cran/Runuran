#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007, Josef Leydold and Wolfgang Hoermann                         ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   DEPRECATED functions!                                                 ##
##                                                                         ##
#############################################################################

#############################################################################
##                                                                          #
## Special sampling methods                                                 #
##                                                                          #
#############################################################################

## -- HINV: Hermite interpolation for approximate INVersion -----------------
## DEPRECATED!
## use function 'uq' instead.

uqhinv <- function (unr, U) { uq(unr,U) }

## -- DAU: Alias-Urn Method ------------------------------------------------
## DEPRECATED!
## use function 'ur(dau.new(...),n)' instead

urdau <- function (n, probvector, from = 0, by = 1) {
        ## create UNU.RAN object
        unr <- dau.new(pv=probvector, from=0)
        ## draw sample
        if (from==0 && by==1)
                unuran.sample(unr,n)
        else
                from + by * unuran.sample(unr,n)
}

## -- DGT: Guide Table Method -----------------------------------------------
## DEPRECATED!
## use function 'ur(dgt.new(...),n)' instead

urdgt <- function (n, probvector, from = 0, by = 1) {
        ## create UNU.RAN object
        unr <- dgt.new(pv=probvector, from=0)
        ## draw sample
        if (from==0 && by==1)
                unuran.sample(unr,n)
        else
                from + by * unuran.sample(unr,n)
}

## -- HITRO: Hit-and-Run sampler with Ratio-of-Uniforms ---------------------
## DEPRECATED!
## use function 'ur(hitro.new(...),n)' instead

urhitro <- function (n, dim=1, pdf, mode=NULL, center=NULL, ll=NULL, ur=NULL, thinning=1, burnin=0, ...) {
        ## create UNU.RAN object
        unr <- hitro.new(dim=dim, pdf=pdf, mode=mode, center=center, ll=ll, ur=ur,
                         thinning=thinning, burnin=burnin, ...)
        ## draw sample
        unuran.sample(unr,n)
}

## -- TDR: Transformed Density Rejection ------------------------------------
## DEPRECATED!
## use function 'ur(tdr.new(...),n)' instead.

urtdr <- function (n, pdf, dpdf=NULL, lb=-Inf, ub=Inf, islog=TRUE, ...) {
        ## create UNU.RAN object
        unr <- tdr.new(pdf=pdf,dpdf=dpdf,lb=lb,ub=ub,islog=islog,...)
        ## draw sample
        unuran.sample(unr,n)
}


#############################################################################
##                                                                          #
## Create distribution objects                                              #
##                                                                          #
#############################################################################

## -- Contiuous Multivariate Distributions ----------------------------------
## DEPRECATED!
## use function 'unuran.cmv.new(...)' instead.

unuran.cmv <- function(dim=1, pdf=NULL, mode=NULL, center=NULL, ll=NULL, ur=NULL) {
        new("unuran.cmv", dim=dim, pdf=pdf, mode=mode, center=center, ll=ll, ur=ur)
}

## -- Contiuous Distributions -----------------------------------------------
## DEPRECATED!
## use function 'unuran.cont.new(...)' instead.

unuran.cont <- function(cdf=NULL, pdf=NULL, dpdf=NULL, islog=TRUE,
                        mode=NA, center=NA, lb=-Inf, ub=Inf, area=NA) {
        new("unuran.cont", cdf=cdf, pdf=pdf, dpdf=dpdf, islog=islog,
            mode=mode, center=center, lb=lb, ub=ub, area=area)
}

## -- Discrete Distributions ------------------------------------------------
## DEPRECATED!
## use function 'unuran.discr.new(...)' instead.

unuran.discr <- function(pv=NULL, pmf=NULL, mode=NA, lb=0, ub=Inf, sum=NA) {
        new("unuran.discr", pv=pv, pmf=pmf, mode=mode, lb=lb, ub=ub, sum=sum)
}

## -- End -------------------------------------------------------------------
