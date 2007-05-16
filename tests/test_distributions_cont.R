#############################################################################
##                                                                         ##
##   Tests for wrapper functions for special distributions                 ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Continuous univariate distributions                                   ##
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

## --- Function for running chi^2 goodness-of-fit test ----------------------

unur.test.cont <- function (distr, domain=NULL, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  Repeat test once if it fails
        ##  (we do not check the validity of the algorithm
        ##   but the correct implementation.)
        ##
        ##  distr  ... name of distribution (as used for [rpqd]... calls
        ##  domain ... domain of distribution
        ##  '...'  ... list of parameters of distribution
        ##

        ## -- domain ?
        have.domain = ifelse( is.null(domain), FALSE, TRUE )
        lb <- ifelse( have.domain, domain[1], -Inf)
        ub <- ifelse( have.domain, domain[2], Inf)
        
        ## -- text for calling function
        cat("ur",distr,"(",paste("*",...,sep=","),")",
                      ifelse( have.domain,paste(" domain=(",signif(lb),",",signif(ub),"): ",sep=""),": "),
                      sep="")

        ## -- function for generating a random sample
        rfuncname <- paste("ur",distr,sep="")
        if (!exists(rfuncname))
                stop("undefined function '",rfuncname,"'")
        rfunc <- match.fun(rfuncname, descend=FALSE)
        ## -- function for computing CDF
        pfuncname <- paste("p",distr,sep="")
        if (!exists(pfuncname))
                stop("undefined function '",pfuncname,"'")
        pfunc <- match.fun(pfuncname, descend=FALSE)

        ## -- run test and repeat once when failed the first time
        for (i in 1:2) {

                ## -- random sample
                if (have.domain)
                        x <- rfunc(samplesize,lb=lb,ub=ub,...)
                else
                        x <- rfunc(samplesize,...)
                ## -- test domain (if given)
                if (have.domain) {
                        too.small <- length(x[x<lb])
                        too.large <- length(x[x>ub])
                        if (too.small > 1 || too.small > 1) {
                                too.small <- 100*too.small/samplesize
                                too.large <- 100*too.large/samplesize
                                stop("X out of domain (",too.small,"%|",too.large,"%)", call.=FALSE)
                        }
                }
                
                ## -- transform into uniform distribution
                u <- pfunc(x,...)
                ## -- make histogram of with equalsized bins (classified data)
                nbins <- as.integer(sqrt(samplesize))
                breaks <- (0:nbins)/nbins
                if (have.domain) {
                        u.lb = pfunc(lb,...)
                        u.ub = pfunc(ub,...)
                        breaks <- u.lb + breaks*(u.ub-u.lb)
                        breaks[length(breaks)] <- u.ub
                }
                h <- hist(u,plot=F,breaks=breaks)$count
                ## -- run unur.chiq.test --
                pval <- chisq.test(h)$p.value
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
## Continuous univariate Distributions                                      #
#############################################################################

## Beta distribution - (replacement for rbeta) ------------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("beta", shape1=runif(1,0.1,10), shape2=runif(1,0.1,10), 
                       domain=sort(runif(2)))
for (i in 1:n.rep.params)
        unur.test.cont("beta", shape1=runif(1,0.1,10), shape2=runif(1,0.1,10))

## Burr family of distributions ---------------------------------------------

## Cauchy distribution - (replacement for rcauchy) --------------------------
unur.test.cont("cauchy")
for (i in 1:n.rep.domains)
        unur.test.cont("cauchy", domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("cauchy", location=rcauchy(1), scale=rgamma(1,shape=2))

## Chi distribution ---------------------------------------------------------

## Chi^2 distribution - (replacement for rchisq) ----------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("chisq", df=3, domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("chisq", df=rgamma(1,shape=2))

## Exponential distribution - (replacement for rexp) ------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("exp", domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("exp", rate=runif(1,0.1,10))

## Extreme value type I (Gumbel-type) distribution --------------------------

## Extreme value type II (Frechet-type) distribution ------------------------

## F distribution  - (replacement for rf) -----------------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("f", df1=runif(1,0.1,10), df2=runif(1,0.1,10), 
                       domain=sort(runif(2)))
for (i in 1:n.rep.params)
        unur.test.cont("f", df1=runif(1,0.1,10), df2=runif(1,0.1,10))

## Gamma distribution  - (replacement for rgamma) ---------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("gamma", shape=runif(1,0.1,10), domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("gamma", shape=runif(1,0.1,10), scale=runif(1,0.1,10))

## Generalized inverse Gaussian ---------------------------------------------

## Hyperbolic distribution --------------------------------------------------

## Laplace (double exponential) distribution --------------------------------

## Lognormal distribution  - (replacement for rlnorm) -----------------------
unur.test.cont("lnorm")
for (i in 1:n.rep.domains)
        unur.test.cont("lnorm", domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("lnorm", meanlog=rnorm(1), sdlog=rgamma(1,shape=2))

## Logistic distribution - (replacement for rlogistic) ----------------------
unur.test.cont("logis")
for (i in 1:n.rep.domains)
        unur.test.cont("logis", domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("logis", location=rcauchy(1), scale=rgamma(1,shape=2))

## Lomax distribution (Pareto distribution of second kind) ------------------

## Normal (Gaussian) distribution - (replacement for rnorm) -----------------
unur.test.cont("norm")
for (i in 1:n.rep.domains)
        unur.test.cont("norm", domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("norm", mean=rcauchy(1), sd=rgamma(1,shape=2))

## Pareto distribution ------------------------------------------------------

## Planck distribution ------------------------------------------------------

## Powerexponential (Subbotin) distribution ---------------------------------

## Rayleigh distribution ----------------------------------------------------

## Student's t distribution - (replacement for rt) --------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("t", df=2, domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("t", df=rgamma(1,shape=2))

## Weibull distribution - (replacement for rweibull) ------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("weibull", shape=runif(1,0.1,10), domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("weibull", shape=runif(1,0.1,10), scale=runif(1,0.1,10))


## -- Print statistics ------------------------------------------------------

cat("\nTests for continuous univariate distributions\n\n",
    "\tnumber of tests = ",length(pvals),"  (number of warnings = ",n.warns,")\n\n")
summary(pvals)


## -- End -------------------------------------------------------------------
