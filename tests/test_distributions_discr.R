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

## --- Run tests? -----------------------------------------------------------

unur.run.tests <- FALSE   ## <-- change to FALSE if tests should not be performed

if (!unur.run.tests) {
        cat("\nRunuran tests not performed!\n\n")
        quit(save="no",status=0,runLast=FALSE)
}

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e5

## level of significance
alpha <- 1.e-3

## number of repetitions
n.rep.domains <- 2
n.rep.params  <- 5

## --- Global counters ------------------------------------------------------

pvals <- numeric(0)
n.warns <- 0

## --- Load library ---------------------------------------------------------

library(Runuran)

## --- DISCR: Function for running chi^2 goodness-of-fit test ---------------

unur.test.discr <- function (distr, domain, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  Repeat test once if it fails
        ##  (we do not check the validity of the algorithm
        ##   but the correct implementation.)
        ##
        ##  distr  ... name of distribution (as used for [rpqd]... calls
        ##  domain ... domain of distribution
        ##  '...'  ... list of parameters of distribution
        ##

        ## -- domain 
        lb <- domain[1]
        ub <- domain[2]
        
        ## -- text for calling function
        cat("ur",distr,"(",paste("*",...,sep=","),
            ") domain=(",signif(lb),",",signif(ub),"): ",sep="")

        ## -- function for generating a random sample
        rfuncname <- paste("ur",distr,sep="")
        if (!exists(rfuncname))
                stop("undefined function '",rfuncname,"'")
        rfunc <- match.fun(rfuncname, descend=FALSE)
        ## -- function for computing probability vector
        dfuncname <- paste("d",distr,sep="")
        if (!exists(dfuncname))
                stop("undefined function '",dfuncname,"'")
        dfunc <- match.fun(dfuncname, descend=FALSE)

        ## -- run test and repeat once when failed the first time
        for (i in 1:2) {
                
                ## -- create probability vector
                probs <- dfunc(lb:ub,...)

                ## -- random sample
                x <- rfunc(samplesize,lb=lb,ub=ub,...)

                ## -- test domain
                too.small <- length(x[x<lb])
                too.large <- length(x[x>ub])
                if (too.small > 1 || too.small > 1) {
                        too.small <- 100*too.small/samplesize
                        too.large <- 100*too.large/samplesize
                        stop("X out of domain (",too.small,"%|",too.large,"%)", call.=FALSE)
                }

                ## -- random distribution
                hits <- hist(x,breaks=(lb:(ub+1)-0.5),plot=FALSE)$counts

                ## -- ignore some of the classes
                probs <- probs[hits>5]
                hits <- hits[hits>5]

                ## -- run unur.chiq.test --
                pval <- chisq.test(hits,p=probs,rescale.p=TRUE)$p.value

                ## -- check p-value
                if (pval > alpha) { # test passed
                        message("chisq test PASSed with p-value=",signif(pval))
                        break
                }

                ## -- test failed
                if (i>1) { # second try failed
                        stop("chisq test FAILED!  p-value=",signif(pval), call.=FALSE)
                }
                else {
                        warning("first chisq test FAILed with p-value=",signif(pval),
                                call.=FALSE,immediate.=TRUE)
                        assign("n.warns", n.warns+1, env = .GlobalEnv)

                }
        }

        ## -- update list of p-values
        assign("pvals", append(pvals, pval), env = .GlobalEnv)
}

#############################################################################
## Discrete univariate Distributions                                        #
#############################################################################

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

cat("\nTests for discrete distributions\n\n",
    "\tnumber of tests = ",length(pvals),"  (number of warnings = ",n.warns,")\n\n")
summary(pvals)


## -- End -------------------------------------------------------------------
