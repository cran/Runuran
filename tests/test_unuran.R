
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

## --- Continuous distributions ---------------------------------------------

## Create an object
unr <- new("unuran", "normal()")

## Print object
unr
print(unr)

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)

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

## --- Continuous distributions - S4 distribution object --------------------

## use PDF
gausspdf <- function (x) { exp(-0.5*x^2) }
gaussdpdf <- function (x) { -x*exp(-0.5*x^2) }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, islog=FALSE)
unr <- 0; unr <- unuran.new(gauss, "tdr")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## use logPDF
gausspdf <- function (x) { -0.5*x^2 }
gaussdpdf <- function (x) { -x }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, islog=TRUE)
unr <- 0; unr <- unuran.new(gauss, "tdr")
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
binom <- new("unuran.discr",pv=pv)
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

## --- End ------------------------------------------------------------------
