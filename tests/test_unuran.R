
#############################################################################
##                                                                         ##
##   Tests for class 'unuran'                                              ##
##                                                                         ##
#############################################################################

## --- Test Parameters ------------------------------------------------------

samplesize <- 1e4

## --- Load library ---------------------------------------------------------

library(Runuran)

## --- Continuous distributions ---------------------------------------------

## Create an object
unr <- new("unuran", "normal()")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)

## Run a chi-square GoF test
chisq.test( hist(pnorm(x),plot=FALSE)$density )


## Create an object
unr <- unuran.new("normal()")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)

## Run a chi-square GoF test
chisq.test( hist(pnorm(x),plot=FALSE)$density )

## --- Discrete distributions -----------------------------------------------

## Create an object
unr <- new("unuran", "binomial(20,0.5)", "dgt")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)

## --- End ------------------------------------------------------------------
