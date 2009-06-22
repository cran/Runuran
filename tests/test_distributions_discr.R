#############################################################################
##                                                                         ##
##   Tests for wrapper functions for special distributions                 ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Discrete distributions                                                ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Remark: You must use named arguments when calling the test routines!  ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- Chi^2 goodness-of-fit test -------------------------------------------

## Binomial distribution - (replacement for rbinom) -------------------------
for (i in 1:n.rep.domains) {
        d <- sort(rbinom(2,size=200,prob=0.29))
        d[2] <- d[2]+1   ## protect agains domains of length 1
        unur.test.discr("binom", size=200, prob=0.3, domain=d)
}
for (i in 1:n.rep.params) {
        s <- as.integer(runif(1,min=10,max=1000))
        p <- runif(1,min=0.01,max=0.99)
        unur.test.discr("binom", size=s, prob=p, domain=c(0,s))
}

## Geometric distribution - (replacement for rgeom) -------------------------
for (i in 1:n.rep.domains) {
        d <- sort(rbinom(2,size=200,prob=0.29))
        d[2] <- d[2]+1   ## protect agains domains of length 1
        unur.test.discr("geom", prob=0.1, domain=d)
}
for (i in 1:n.rep.params) {
        p <- runif(1,min=0.03,max=0.99)
        unur.test.discr("geom", prob=p, domain=c(0,1000))
}
for (i in 1:n.rep.params) {
        p <- runif(1,min=0.001,max=0.02)
        unur.test.discr("geom", prob=p, domain=c(0,1000))
}


## -- Print statistics ------------------------------------------------------

unur.test.statistic()

## -- End -------------------------------------------------------------------

detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
