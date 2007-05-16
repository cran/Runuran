#############################################################################
##                                                                         ##
##   Tests for methods                                                     ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- CONT: Chi^2 goodness-of-fit test -------------------------------------

## TDR (Transformed Density Rejection)
tdr.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*x^2) }
        dpdf <- function (x) { -x*exp(-0.5*x^2) }
        dist <- new("unuran.cont", pdf=pdf, dpdf=dpdf, islog=FALSE)
        gen <- unuran.new(dist, "tdr")
        unuran.sample(gen,n)
}
unur.test.cont("tdr.norm", rfunc=tdr.norm, pfunc=pnorm)

urtdr.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*x^2) }
        dpdf <- function (x) { -x*exp(-0.5*x^2) }
        urtdr(n, pdf=pdf, dpdf=dpdf, islog=FALSE)
}
unur.test.cont("urtdr.norm", rfunc=urtdr.norm, pfunc=pnorm)

urtdr.norm.wod <- function (n) {
        pdf <- function (x) { exp(-0.5*x^2) }
        urtdr(n, pdf=pdf, islog=FALSE)
}
unur.test.cont("urtdr.norm.wod", rfunc=urtdr.norm.wod, pfunc=pnorm)

tdr.norm.wl <- function (n) {
        logpdf <- function (x) { -0.5*x^2 }
        dlogpdf <- function (x) { -x }
        dist <- new("unuran.cont", pdf=logpdf, dpdf=dlogpdf, islog=TRUE)
        gen <- unuran.new(dist, "tdr")
        unuran.sample(gen,n)
}
unur.test.cont("tdr.norm.wl", rfunc=tdr.norm.wl, pfunc=pnorm)

urtdr.norm.wl <- function (n) {
        logpdf <- function (x) { -0.5*x^2 }
        dlogpdf <- function (x) { -x }
        urtdr(n, pdf=logpdf, dpdf=dlogpdf)
}
unur.test.cont("urtdr.norm.wl", rfunc=urtdr.norm.wl, pfunc=pnorm)

urtdr.norm.wlwod <- function (n) {
        logpdf <- function (x) { -0.5*x^2 }
        dlogpdf <- function (x) { -x }
        urtdr(n, pdf=logpdf)
}
unur.test.cont("urtdr.norm.wlwod", rfunc=urtdr.norm.wlwod, pfunc=pnorm)


## --- DISCR: Chi^2 goodness-of-fit test ------------------------------------

## DGT (Discrete Guide Table method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
dgt.binom <- function (n,lb=0,ub=size) {
        dist <- new("unuran.discr", pv=binom.probs)
        gen <- unuran.new(dist, "dgt")
        unuran.sample(gen,n)
}
unur.test.discr("dgt.binom", rfunc=dgt.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("dgt.binom", rfunc=dgt.binom, pv=binom.probs, domain=c(0,size))

urdgt.binom <- function (n,lb=0,ub=size) {
        urdgt(n, probvector=binom.probs)
}
unur.test.discr("urdgt.binom", rfunc=urdgt.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("urdgt.binom", rfunc=urdgt.binom, pv=binom.probs, domain=c(0,size))

## DAU (Discrete Alias-Urn method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
dau.binom <- function (n,lb=0,ub=size) {
        dist <- new("unuran.discr", pv=binom.probs)
        gen <- unuran.new(dist, "dau")
        unuran.sample(gen,n)
}
unur.test.discr("dau.binom", rfunc=dau.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("dau.binom", rfunc=dau.binom, pv=binom.probs, domain=c(0,size))

urdau.binom <- function (n,lb=0,ub=size) {
        urdau(n, probvector=binom.probs)
}
unur.test.discr("urdau.binom", rfunc=urdau.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("urdau.binom", rfunc=urdau.binom, pv=binom.probs, domain=c(0,size))

## --- CMV: Chi^2 goodness-of-fit test --------------------------------------

samplesize <- 1.e4

## HITRO (Hit-and-Run + Ratio-of-Uniforms)
hitro.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*sum(x^2)) }
        dist <- new("unuran.cmv", dim=2, pdf=pdf)
        gen <- unuran.new(dist, "hitro; thinning=10")
        unuran.sample(gen,n)
}
unur.test.cmv("hitro.norm", rfunc=hitro.norm, pfunc=pnorm)

urhitro.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*sum(x^2)) }
        urhitro(n, dim=2, pdf=pdf, mode=c(0,0), thinning=10)
}
unur.test.cmv("urhitro.norm", rfunc=urhitro.norm, pfunc=pnorm)

## VNROU (Naive  multivariate Ratio-of-Uniforms method)
vnrou.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*sum(x^2)) }
        dist <- new("unuran.cmv", dim=2, pdf=pdf)
        gen <- unuran.new(dist, "vnrou")
        unuran.sample(gen,n)
}
unur.test.cmv("vnrou.norm", rfunc=vnrou.norm, pfunc=pnorm)

## -- Print statistics ------------------------------------------------------

unur.test.statistic()

## -- End -------------------------------------------------------------------
