
#############################################################################
##                                                                         ##
##   Tests for class 'unuran'                                              ##
##                                                                         ##
#############################################################################

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e6

## break points for chi^2 GoF test
nbins <- as.integer(sqrt(samplesize))
breaks <- (0:nbins)/nbins

## level of significance
alpha <- 1.e-3

## --- Load library ---------------------------------------------------------

library(Runuran)


#############################################################################
##                                                                          #
##  Auxiliary routines                                                      #
##                                                                          #
#############################################################################

## -- Auxiliary routines ------------------------------------------------------

## Test whether there is an error -------------------------------------------
is.error <- function (expr) { is(try(expr), "try-error") }


## --- Continuous distributions ---------------------------------------------

## Create an object
unr <- new("unuran", "normal()")

## Print object
unr
print(unr)
unuran.details(unr)

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)
ur(unr)
ur(unr,10)

## Run a chi-square GoF test
chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )

## Create an object
unr <- unuran.new("normal()")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)

## Run a chi-square GoF test
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## another example (for testing print)
unr <- new("unuran", "normal()", "arou")
unr
print(unr)
unuran.details(unr)


## --- Continuous distributions - S4 distribution object --------------------

## use PDF
gausspdf <- function (x) { exp(-0.5*x^2) }
gaussdpdf <- function (x) { -x*exp(-0.5*x^2) }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, lb=-Inf, ub=Inf, center=0.1)
unr <- 0; unr <- unuran.new(gauss, "tdr")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## use PDF
gausspdf <- function (x) { exp(-0.5*x^2) }
gaussdpdf <- function (x) { -x*exp(-0.5*x^2) }
gauss <- unuran.cont.new(pdf=gausspdf, dpdf=gaussdpdf, lb=-Inf, ub=Inf, center=0.1)
unr <- 0; unr <- unuran.new(gauss, "tdr")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## use logPDF
gausspdf <- function (x) { -0.5*x^2 }
gaussdpdf <- function (x) { -x }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, islog=TRUE, lb=-Inf, ub=Inf, mode=0)
unr <- 0; unr <- unuran.new(gauss, "tdr")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## use logPDF (use ARS to test 'print' function)
gausspdf <- function (x) { -0.5*x^2 }
gaussdpdf <- function (x) { -x }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, islog=TRUE, lb=-Inf, ub=Inf)
unr <- 0; unr <- unuran.new(gauss, "ars")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))


## --- Discrete distributions -----------------------------------------------

## Create an object
unr <- 0;
unr <- new("unuran", "binomial(20,0.5)", "dgt")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)


## --- Discrete distributions - S4 distribution object ----------------------

## use PV
pv <- dbinom(0:100,100,0.3)
binom <- new("unuran.discr",pv=pv,lb=0)
unr <- 0; unr <- unuran.new(binom, "dgt")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pbinom(x,100,0.3),plot=FALSE)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## use PMF
pmf <- function(x) dbinom(x,100,0.3)
binom <- new("unuran.discr",pmf=pmf,lb=0,ub=100)
unr <- 0; unr <- unuran.new(binom, "dgt")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pbinom(x,100,0.3),plot=FALSE)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## use PMF
pmf <- function(x) dbinom(x,100,0.3)
binom <- new("unuran.discr",pmf=pmf,lb=0,ub=100)
unr <- 0; unr <- unuran.new(binom, "dari")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pbinom(x,100,0.3),plot=FALSE)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))


## --- Continuous Multivariate distributions --------------------------------

mvpdf <- function (x) { exp(-sum(x^2)) }
mvd <- new("unuran.cmv", dim=2, pdf=mvpdf, mode=c(0,0))
unr <- 0; unr <- unuran.new(mvd, "hitro")
x <- unuran.sample(unr, 10)
x

unr <- 0; unr <- unuran.new(mvd, "vnrou")
x <- unuran.sample(unr, 10)
x


## --- quantile function ----------------------------------------------------

## test U-error
unr <- unuran.new("normal()","hinv; u_resolution=1.e-13")

## single double as argument
Tmax <- 0
for (U in (0:20)/20) {
        T <- pnorm( uq(unr,U) ) - U
        print (T)
        Tmax <- max(abs(T),Tmax)
}
cat("Max. error =",Tmax,"\n")
if (Tmax > 1.e-13) stop ("Max. error exceeds limit")

## vector argument
U <- (0:20)/20
T <- pnorm( uq(unr,U) ) - U
print (T)
Tmax <- max(abs(T))
cat("Max. error =",Tmax,"\n")
if (Tmax > 1.e-13) stop ("Max. error exceeds limit")

## special arguments
for (U in c(-1,-0.001,0,0.5,1.,1.001,NA,NaN) ) {
        cat ("U =",U,"\tX =",uq(unr,U),"\n")
}

U <- c(-1,-0.001,0,0.5,1.,1.001,NA,NaN)
T <- uq(unr,U)
rbind(U,T)

uq(unr,numeric())

## test whether 'uq' throws an error when UNU.RAN object does not implement
## an inversion method
unr <- unuran.new("normal()","tdr")
if( ! is.error( uq(unr,0.5) ) )
   stop("'uq' does not detect UNU.RAN object with non-inversion method")

## test whether 'uq' detects invalid arguments
if( ! is.error( uq(1,0.5) ) )
   stop("'uq' does not detect invalid argument 'unr'")
        
if( ! is.error( uq(unr,"a") ) )
   stop("'uq' does not detect invalid argument 'U'")
        

## --- End ------------------------------------------------------------------
