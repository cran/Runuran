#############################################################################
##                                                                         ## 
## Test Runuran distribution classes                                       ## 
##                                                                         ## 
#############################################################################
##                                                                         ##
## Interface for classes                                                   ##
##                                                                         ##
##  - "unuran.cont"         ... univariate continuous distributions        ##
##  - "unuran.discr"        ... univariate discrete distributions          ##
##  - "unuran.cmv"          ... multivariate continuous distributions      ##
##                                                                         ##
## These are extensions of the virtual class "unuran.distr".               ##
##                                                                         ##
## Functions / Methods:                                                    ##
##                                                                         ##
##  - new()                 ... create new instance of class               ##
##  - unuran.cont.new()     ...   shortcut for "unuran.cont"               ##
##  - unuran.discr.new()    ...   shortcut for "unuran.discr"              ##
##  - unuran.cmv.new()      ...   shortcut for "unuran.cmv"                ##
##                                                                         ##
##  - ud()  ... density     ["unuran.cont" and "unuran.discr" only]        ##
##  - up()  ... CDF         ["unuran.cont" and "unuran.discr" only]        ##
##  - uq()  ... quantile    [not implemented]                              ##
##  - ur()  ... rng         [not implemented]                              ##
##                                                                         ##
#############################################################################

## --------------------------------------------------------------------------
##
## class "unuran.cont"
##
## --------------------------------------------------------------------------

context("[distr-cont] - unuran.cont.new")

## -- ud() and up(): normal -------------------------------------------------

test_that("[distr-cont-ud] unuran.cont.new() + ud()", {
    ## test: ud() and underlying PDF must give same results.
    distr <- unuran.cont.new(pdf=dnorm, lb=-Inf, ub=Inf)

    ## special values
    x <- c(NA, NaN, -Inf, Inf)
    expect_identical(ud(distr,x), dnorm(x))

    ## generic case
    x <- runif(1000, -5, 5)
    expect_equal(ud(distr,x), dnorm(x))
})

test_that("[distr-cont-up] unuran.cont.new() + up()", {
    ## test: up() and underlying CDF must give same results.
    distr <- unuran.cont.new(cdf=pnorm, lb=-Inf, ub=Inf)

    ## special values
    x <- c(NA, NaN, -Inf, Inf)
    expect_identical(up(distr,x), pnorm(x))

    ## generic case
    x <- runif(1000, -5, 5)
    expect_equal(up(distr,x), pnorm(x))
})

## -- INVALID: unuran.cont.new ----------------------------------------------

test_that("[distr-cont-i01] unuran.cont.new() with invalid arguments", {

    ## test: both 'lb' and 'ub' must be given, numeric, and 'lb'<'ub'
    msg <- "domain \\('lb','ub'\\) missing or invalid"
    expect_error(unuran.cont.new(), msg)
    expect_error(unuran.cont.new(lb=0), msg)
    expect_error(unuran.cont.new(ub=1), msg)
    expect_error(unuran.cont.new(lb=1, ub=0), msg)
    expect_error(unuran.cont.new(lb="a", ub="b"), msg)
    expect_error(unuran.cont.new(lb=0, ub="a"), msg)

    ## test: 'cdf' must be an R function
    msg <- "invalid argument 'cdf'"
    expect_error(unuran.cont.new(cdf=1, lb=0, ub=1), msg)
    
    ## test: 'pdf' must be an R function
    msg <- "invalid argument 'pdf'"
    expect_error(unuran.cont.new(pdf=1, lb=0, ub=1), msg)

    ## test: 'dpdf' must be an R function
    msg <- "invalid argument 'dpdf'"
    expect_error(unuran.cont.new(dpdf=1, lb=0, ub=1), msg)

    ## test: 'islog' must be a boolean
    msg <- "argument 'islog' must be boolean"
    expect_error(unuran.cont.new(islog=1, lb=0, ub=1), msg)

    ## test: 'mode' must be numeric
    msg <- "invalid argument 'mode'"
    expect_error(unuran.cont.new(mode="invalid", lb=0, ub=1), msg)

    ## test: 'center' must be numeric
    msg <- "invalid argument 'center'"
    expect_error(unuran.cont.new(center="invalid", lb=0, ub=1), msg)

    ## test: 'area' must be numeric and strictly positive
    msg <- "invalid argument 'area'"
    expect_error(unuran.cont.new(area="invalid", lb=0, ub=1), msg)
    expect_error(unuran.cont.new(area=0, lb=0, ub=1), msg)
    expect_error(unuran.cont.new(area=-1, lb=0, ub=1), msg)

    ## test: 'name' must be a character string
    msg <- "invalid argument 'name'"
    expect_error(unuran.cont.new(name=1, lb=0, ub=1), msg)
})

## -- INVALID: ud, up, uq, ur -----------------------------------------------

test_that("[distr-cont-i02] invalid call to ud", {
    ## test: ud() requires PDF
    msg <- "\\[UNU\\.RAN - error\\] UNU\\.RAN object does not contain \\(log\\)PDF"
    distr <- unuran.cont.new(cdf=pnorm, dpdf=function(x){-x^2}, lb=-Inf, ub=Inf)
    expect_warning(x <- ud(distr,0), msg)
    expect_equal(x, NA_real_)

    ## invalid x
    distr <- unuran.cont.new(pdf=dnorm, lb=-Inf, ub=Inf)
    msg <- "\\[UNU\\.RAN - error\\] argument invalid: 'x' must be numeric"
    expect_error(ud(distr,"invalid"), msg)
})

test_that("[distr-cont-i03] invalid call to up", {
    ## test: up() requires CDF
    msg <- "\\[UNU\\.RAN - error\\] UNU\\.RAN object does not contain CDF"
    distr <- unuran.cont.new(pdf=dnorm, dpdf=function(x){-x^2}, lb=-Inf, ub=Inf)
    expect_error(up(distr,0), msg)

    ## invalid x
    distr <- unuran.cont.new(cdf=pnorm, lb=-Inf, ub=Inf)
    msg <- "\\[UNU\\.RAN - error\\] argument invalid: 'x' must be numeric"
    expect_error(ud(distr,"invalid"), msg)
})

test_that("[distr-cont-i04] invalid call to uq", {
    ## test: uq() not implemented
    msg <- "argument 'unr' must be UNU.RAN object; method not implemented for distribution objects"
    distr <- unuran.cont.new(pdf=dnorm, dpdf=function(x){-x^2}, lb=-Inf, ub=Inf)
    expect_error(uq(distr,0.5), msg)
})

test_that("[distr-cont-i05] invalid call to ur", {
    ## test: ur() not implemented
    ## (however there is a weird error message in order to save overhead)
    msg <- "no slot of name \"unur\" for this object of class \"unuran.cont\""
    distr <- unuran.cont.new(pdf=dnorm, dpdf=function(x){-x^2}, lb=-Inf, ub=Inf)
    expect_error(ur(distr,1), msg)
})

## --------------------------------------------------------------------------
##
## class "unuran.discr"
##
## --------------------------------------------------------------------------

context("[distr-discr] - unuran.discr.new")

## -- ud() and up(): binomial -----------------------------------------------

test_that("[distr-discr-ud-binom] unuran.discr.new() + ud()", {
    ## test: ud() and underlying PMF must give same results.

    size <- 100; prob <- 0.4 
    distr <- unuran.discr.new(pmf=function(x){dbinom(x,size,prob)}, lb=0, ub=size)

    ## special values
    x <- c(NA, NaN, -Inf, Inf, 0, 1e300, -1e300)
    expect_identical(ud(distr,x), dbinom(x,size,prob))

    ## generic case
    x <- 0:size
    expect_equal(ud(distr,x), dbinom(x,size,prob))
})

test_that("[distr-discr-up-binom] unuran.discr.new() + up()", {
    ## test: up() and underlying CDF must give same results.

    size <- 100; prob <- 0.4 
    distr <- unuran.discr.new(cdf=function(x){pbinom(x,size,prob)}, lb=0, ub=size)

    ## special values
    x <- c(NA, NaN, -Inf, Inf, 0, 1e300, -1e300)
    expect_identical(up(distr,x), pbinom(x,size,prob))

    ## generic case
    x <- 0:size
    expect_equal(up(distr,x), pbinom(x,size,prob))
})

## -- ud() and up(): geometric ----------------------------------------------

test_that("[distr-discr-ud-geom] unuran.discr.new() + ud()", {
    ## test: ud() and underlying PMF must give same results.

    prob <- 0.4
    distr <- unuran.discr.new(pmf=function(x){dgeom(x,prob)}, lb=0, ub=Inf)

    ## special values
    x <- c(NA, NaN, -Inf, Inf, 0, 1e300, -1e300)
    expect_identical(ud(distr,x), dgeom(x,prob))

    ## generic case
    x <- 0:1000
    expect_equal(ud(distr,x), dgeom(x,prob))
})

test_that("[distr-discr-up-geom] unuran.discr.new() + up()", {
    ## test: up() and underlying CDF must give same results.

    prob <- 0.4
    distr <- unuran.discr.new(cdf=function(x){pgeom(x,prob)}, lb=0, ub=Inf)

    ## special values
    x <- c(NA, NaN, -Inf, Inf, 0, 1e300, -1e300)
    expect_identical(up(distr,x), pgeom(x,prob))

    ## generic case
    x <- 0:1000
    expect_equal(up(distr,x), pgeom(x,prob))
})

## -- INVALID: unuran.discr.new ---------------------------------------------

test_that("[distr-discr-i01] unuran.discr.new() with invalid arguments", {

    ## test: both 'lb' and 'ub' must be given, numeric, and 'lb'<'ub'
    msg <- "domain \\('lb','ub'\\) missing or invalid"
    expect_error(unuran.discr.new(), msg)
    expect_error(unuran.discr.new(lb=0), msg)
    expect_error(unuran.discr.new(ub=1), msg)
    expect_error(unuran.discr.new(lb=1, ub=0), msg)
    expect_error(unuran.discr.new(lb="a", ub="b"), msg)
    expect_error(unuran.discr.new(lb=0, ub="a"), msg)

    ## test: if 'pv' is given then, 'lb' must be given and numeric
    expect_error(unuran.discr.new(pv=1:3,ub=0), msg)
    
    ## test: 'cdf' must be an R function
    msg <- "invalid argument 'cdf'"
    expect_error(unuran.discr.new(cdf=1, lb=0, ub=1), msg)
    
    ## test: 'pmf' must be an R function
    msg <- "invalid argument 'pmf'"
    expect_error(unuran.discr.new(pmf=1:10, lb=0, ub=10), msg)

    ## test: 'pv' must be a numeric array
    msg <- "invalid argument 'pv'"
    expect_error(unuran.discr.new(pv=dbinom, lb=0, ub=10), msg)
    expect_error(unuran.discr.new(pv=c(NA,1:3), lb=0, ub=10), msg)

    ## test: 'mode' must be numeric
    msg <- "invalid argument 'mode'"
    expect_error(unuran.discr.new(mode="invalid", lb=0, ub=1), msg)

    ## test: 'sum' must be numeric and strictly positive
    msg <- "invalid argument 'sum'"
    expect_error(unuran.discr.new(sum="invalid", lb=0, ub=1), msg)
    expect_error(unuran.discr.new(sum=0, lb=0, ub=1), msg)
    expect_error(unuran.discr.new(sum=-1, lb=0, ub=1), msg)

    ## test: 'name' must be a character string
    msg <- "invalid argument 'name'"
    expect_error(unuran.discr.new(name=1, lb=0, ub=1), msg)
})

## -- INVALID: ud, up, uq, ur -----------------------------------------------

test_that("[distr-discr-i02] invalid call to ud", {
    ## test: ud() requires PMF
    msg <- "\\[UNU\\.RAN - error\\] UNU\\.RAN object does not contain \\(log\\)PMF"
    distr <- unuran.discr.new(cdf=function(x){pbinom(x,10,0.5)}, lb=0, ub=10)
    expect_warning(x <- ud(distr,1), msg)
    expect_equal(x, NA_real_)

    ## invalid x
    distr <- unuran.discr.new(pmf=function(x){dbinom(x,10,0.5)}, lb=0, ub=10)
    msg <- "\\[UNU\\.RAN - error\\] argument invalid: 'x' must be numeric"
    expect_error(ud(distr,"invalid"), msg)
})

test_that("[distr-discr-i03] invalid call to up", {
    ## test: up() requires CDF
    msg <- "\\[UNU\\.RAN - error\\] UNU\\.RAN object does not contain CDF"
    distr <- unuran.discr.new(pmf=function(x){dbinom(x,10,0.5)}, lb=0, ub=10)
    expect_error(up(distr,1), msg)

    ## invalid x
    distr <- unuran.discr.new(cdf=function(x){pbinom(x,10,0.5)}, lb=0, ub=10)
    msg <- "\\[UNU\\.RAN - error\\] argument invalid: 'x' must be numeric"
    expect_error(ud(distr,"invalid"), msg)
})

test_that("[distr-discr-i04] invalid call to uq", {
    ## test: uq() not implemented
    msg <- "argument 'unr' must be UNU.RAN object; method not implemented for distribution objects"
    distr <- unuran.discr.new(pmf=function(x){dbinom(x,10,0.5)}, lb=0, ub=10)
    expect_error(uq(distr,0.5), msg)
})

test_that("[distr-discr-i05] invalid call to ur", {
    ## test: ur() not implemented
    ## (however there is a weird error message in order to save overhead)
    msg <- "no slot of name \"unur\" for this object of class \"unuran.discr\""
    distr <- unuran.discr.new(pmf=function(x){dbinom(x,10,0.5)}, lb=0, ub=10)
    expect_error(ur(distr,1), msg)
})


## --------------------------------------------------------------------------
##
## class "unuran.cmv"
##
## --------------------------------------------------------------------------

context("[distr-cmv] - unuran.cmv.new")

## -- ud() ------------------------------------------------------------------

## -- INVALID: unuran.cmv.new ---------------------------------------------..

test_that("[distr-cmv-i01] unuran.cmv.new() with invalid arguments", {

    ## test; 'dim' must be an integer, 1 <= 'dim' <= 100000
    msg <- "invalid argument 'dim'"
    expect_error(unuran.cmv.new(dim=0), msg)
    expect_error(unuran.cmv.new(dim=100001), msg)
    expect_error(unuran.cmv.new(dim=2.5), msg)
    expect_error(unuran.cmv.new(dim="invalid"), msg)

    ## test: 'll' must be a numeric array of size 'dim'
    msg <- "invalid argument 'll'"
    expect_error(unuran.cmv.new(dim=2, ll="a"), msg)
    msg <- "argument 'll' must have length 'dim'"
    expect_error(unuran.cmv.new(dim=2, ll=1),   msg)
    expect_error(unuran.cmv.new(dim=2, ll=1:3), msg)

    ## test: 'ur' must be a numeric array of size 'dim'
    msg <- "invalid argument 'ur'"
    expect_error(unuran.cmv.new(dim=2, ur="a"), msg)
    msg <- "argument 'ur' must have length 'dim'"
    expect_error(unuran.cmv.new(dim=2, ur=1),   msg)
    expect_error(unuran.cmv.new(dim=2, ur=1:3), msg)

    ## test: 'll' < 'ur'
    msg <- "arguments 'll' and 'ur' invalid: condition 'll' < 'ur' violated"
    expect_error(unuran.cmv.new(dim=2, ll=c(1,1), ur=c(1,2)), msg)

    ## test: 'pdf' must be an R function
    msg <- "invalid argument 'pdf'"
    expect_error(unuran.cmv.new(pdf=1), msg)

    ## test: 'mode' must be a numeric array of size 'dim'
    msg <- "invalid argument 'mode'"
    expect_error(unuran.cmv.new(dim=2, mode="invalid"), msg)
    msg <- "argument 'mode' must have length 'dim'"
    expect_error(unuran.cmv.new(dim=2, mode=1), msg)
    expect_error(unuran.cmv.new(dim=2, mode=1:3), msg)

    ## test: 'center' must be a numeric array of size 'dim'
    msg <- "invalid argument 'center'"
    expect_error(unuran.cmv.new(dim=2, center="invalid"), msg)
    msg <- "argument 'center' must have length 'dim'"
    expect_error(unuran.cmv.new(dim=2, center=1), msg)
    expect_error(unuran.cmv.new(dim=2, center=1:3), msg)

    ## test: 'name' must be a character string
    msg <- "invalid argument 'name'"
    expect_error(unuran.cmv.new(name=1), msg)
})

## -- INVALID: ud, up, uq, ur -----------------------------------------------

test_that("[distr-cmv-i02] invalid call to ud", {
    ## test: not implemented
    msg <- "method not implemented for objects of class 'unuran.cmv'"
    distr <- unuran.cmv.new(dim=2, pdf=function(x){exp(-sum(x^2))})
    expect_error(ud(distr,0), msg)
})

test_that("[distr-cmv-i03] invalid call to up", {
    ## test: not implemented
    msg <- "method not implemented for objects of class 'unuran.cmv'"
    distr <- unuran.cmv.new(dim=2, pdf=function(x){exp(-sum(x^2))})
    expect_error(up(distr,0), msg)
})

test_that("[distr-cmv-i04] invalid call to uq", {
    ## test: not implemented
    msg <- "argument 'unr' must be UNU.RAN object; method not implemented for distribution objects"
    distr <- unuran.cmv.new(dim=2, pdf=function(x){exp(-sum(x^2))})
    expect_error(uq(distr,0), msg)
})

test_that("[distr-cmv-i05] invalid call to ur", {
    ## test: ur() not implemented
    ## (however there is a weird error message in order to save overhead)
    msg <- "no slot of name \"unur\" for this object of class \"unuran.cmv\""
    distr <- unuran.cmv.new(dim=2, pdf=function(x){exp(-sum(x^2))})
    expect_error(ur(distr,1), msg)
})

## -- End -------------------------------------------------------------------

