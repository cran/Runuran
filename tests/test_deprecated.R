#############################################################################
##                                                                         ##
##   Tests for deprecated function                                         ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- CONT: Call function --------------------------------------------------

## test deprecated function 'uqhinv'
unr <- unuran.new("normal()","hinv")
uqhinv(unr,(0:20)/20)


## --- CONT: Chi^2 goodness-of-fit test -------------------------------------

## TDR (Transformed Density Rejection)
urtdr.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*x^2) }
        dpdf <- function (x) { -x*exp(-0.5*x^2) }
        urtdr(n, pdf=pdf, dpdf=dpdf, islog=FALSE)
}
unur.test.cont("urtdr.norm", rfunc=urtdr.norm, pfunc=pnorm)
rm(urtdr.norm)

urtdr.norm.wod <- function (n) {
        pdf <- function (x) { exp(-0.5*x^2) }
        urtdr(n, pdf=pdf, islog=FALSE)
}
unur.test.cont("urtdr.norm.wod", rfunc=urtdr.norm.wod, pfunc=pnorm)
rm(urtdr.norm.wod)

urtdr.norm.wl <- function (n) {
        logpdf <- function (x) { -0.5*x^2 }
        dlogpdf <- function (x) { -x }
        urtdr(n, pdf=logpdf, dpdf=dlogpdf)
}
unur.test.cont("urtdr.norm.wl", rfunc=urtdr.norm.wl, pfunc=pnorm)
rm(urtdr.norm.wl)

urtdr.norm.wlwod <- function (n) {
        logpdf <- function (x) { -0.5*x^2 }
        dlogpdf <- function (x) { -x }
        urtdr(n, pdf=logpdf)
}
unur.test.cont("urtdr.norm.wlwod", rfunc=urtdr.norm.wlwod, pfunc=pnorm)
rm(urtdr.norm.wlwod)

## test arguments passed to PDF
urtdr.norm.param <- function (n) {
        pdf <- function (x,a) { exp(a*x^2) }
        urtdr(n, pdf=pdf, islog=FALSE, a=-1/2)
}
unur.test.cont("urtdr.norm.param", rfunc=urtdr.norm.param, pfunc=pnorm)
rm(urtdr.norm.param)

urtdr.norm.R <- function (n) {
        urtdr(n, pdf=dnorm, islog=FALSE)
}
unur.test.cont("urtdr.norm.R", rfunc=urtdr.norm.R, pfunc=pnorm)
rm(urtdr.norm.R)

urtdr.t.R <- function (n,df) {
        urtdr(n, pdf=dt, islog=FALSE,df=df)
}
unur.test.cont("urtdr.t.R", rfunc=urtdr.t.R, pfunc=pt, df=8)
rm(urtdr.t.R)

## --- DISCR: Chi^2 goodness-of-fit test ------------------------------------

## DGT (Discrete Guide Table method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
urdgt.binom <- function (n,lb=0,ub=size) {
        urdgt(n, probvector=binom.probs)
}
unur.test.discr("urdgt.binom", rfunc=urdgt.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("urdgt.binom", rfunc=urdgt.binom, pv=binom.probs, domain=c(0,size))
rm(urdgt.binom)
rm(size,prob,binom.pmf,binom.probs)

## DAU (Discrete Alias-Urn method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
urdau.binom <- function (n,lb=0,ub=size) {
        urdau(n, probvector=binom.probs)
}
unur.test.discr("urdau.binom", rfunc=urdau.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("urdau.binom", rfunc=urdau.binom, pv=binom.probs, domain=c(0,size))
rm(urdau.binom)
rm(size,prob,binom.pmf,binom.probs)


## --- CMV: Chi^2 goodness-of-fit test --------------------------------------

samplesize <- 1.e4

## HITRO (Hit-and-Run + Ratio-of-Uniforms)
urhitro.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*sum(x^2)) }
        urhitro(n, dim=2, pdf=pdf, mode=c(0,0), thinning=10)
}
unur.test.cmv("urhitro.norm", rfunc=urhitro.norm, pfunc=pnorm)
rm(urhitro.norm)

## -- Print statistics ------------------------------------------------------

unur.test.statistic()

## -- End -------------------------------------------------------------------
