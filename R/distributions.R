#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007 Josef Leydold and Wolfgang Hoermann                          ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Wrapper for standard distributions                                    ##
##                                                                         ##
#############################################################################
##   Remark: Please sort function calls alphabetically!                    ##
#############################################################################


#############################################################################
## Continuous univariate Distributions                                      #
#############################################################################

## Beta distribution - (replacement for rbeta) ------------------------------
urbeta <- function (n,shape1,shape2,lb=0,ub=1) {
        unr<-new("unuran", paste("beta(",shape1,",",shape2,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)	
}

## Burr distribution ---------------------------------------------

urburr <- function (n,a,b,lb=0,ub=Inf) {
## works in theor for a >= 1 and b >= 2 
## numerical problems for a*b > 175 or so
        unr <- new("unuran", paste("distr=cont;pdf='",a*(b-1),"*x^(",a-1,")/(1+x^",a,")^",b,"'; domain=(",lb,",",ub,")"),"TDR")
        unuran.sample(unr,n)
}
                                                  
## Cauchy distribution - (replacement for rcauchy) --------------------------
urcauchy <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr<-new("unuran",paste("cauchy(",location,",",scale,"); domain=(",lb,",",ub,")"),"HINV")
        unuran.sample(unr,n)	
}

## Chi distribution ---------------------------------------------------------
urchi <- function (n,df,lb=0,ub=Inf) {
        unr <- new("unuran", paste("chi(",df,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Chi^2 distribution - (replacement for rchisq) ----------------------------
urchisq <- function (n,df,lb=0,ub=Inf) {
        unr <- new("unuran", paste("chisquare(",df,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Exponential distribution - (replacement for rexp) ------------------------
urexp <- function (n,rate=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("exponential(",1./rate,"); domain=(",lb,",",ub,")"), "CSTD")
        unuran.sample(unr,n)
}

## Extreme value type I (Gumbel-type) distribution --------------------------
urextremeI <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("extremeI(",location,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Extreme value type II (Frechet-type) distribution ------------------------
urextremeII <- function (n,shape,location=0,scale=1,lb=location,ub=Inf) {
        unr <- new("unuran", paste("extremeII(",shape,",",location,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## F distribution  - (replacement for rf) -----------------------------------
urf <- function (n,df1,df2,lb=0,ub=Inf) {
        unr <- new("unuran", paste("F(",df1,",",df2,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Gamma distribution  - (replacement for rgamma) ---------------------------
urgamma <- function (n,shape,scale=1,lb=0,ub=Inf) {
        unr<-new("unuran", paste("gamma(",shape,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Generalized inverse Gaussian ---------------------------------------------
urgig <- function (n,lambda,omega,lb=1.e-12,ub=Inf) { 
        ## works for lambda>=1 and omega>0 and for lambda>0 and omega>=0.5
        unr<-new("unuran",
                 paste("cont; pdf='x^(",lambda-1,")*exp(-(",omega/2,")*(x+1/x))'; domain=(",lb,",",ub,")"),"TDR")
        unuran.sample(unr,n)
}

## Hyperbolic distribution --------------------------------------------------
urhyperbolic <- function (n,shape,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran",
                   paste("cont; pdf='exp(-",shape,"*sqrt(1.+x*x/",(scale*scale),"))'; domain=(",lb,",",ub,")"),
                   "TDR")
	 unuran.sample(unr,n)
}

## Laplace (double exponential) distribution --------------------------------
urlaplace <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("laplace(",location,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Lognormal distribution  - (replacement for rlnorm) -----------------------
urlnorm <- function (n,meanlog=0,sdlog=1,lb=0,ub=Inf) {
        exp(urnorm(n,meanlog,sdlog,log(lb),log(ub)))
}

## Logistic distribution - (replacement for rlogistic) ----------------------
urlogis <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("logistic(",location,",",scale,"); domain=(",lb,",",ub,")"), "CSTD")
        unuran.sample(unr,n)
}

## Lomax distribution (Pareto distribution of second kind) ------------------
urlomax <- function (n,shape,scale=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("lomax(",shape,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Normal (Gaussian) distribution - (replacement for rnorm) -----------------
urnorm <- function (n,mean=0,sd=1,lb=-Inf,ub=Inf) {
        unr<-new("unuran",paste("normal(",mean,",",sd,"); domain=(",lb,",",ub,")"),"HINV")
        unuran.sample(unr,n)
}

## Pareto distribution ------------------------------------------------------
urpareto <- function (n,k,a,lb=k,ub=Inf) {
        unr <- new("unuran", paste("pareto(",k,",",a,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Planck distribution ------------------------------------------------------
urplanck <- function (n,a,lb=1.e-12,ub=Inf) { 
        ## works for a>=1 
        unr <- new("unuran", paste("cont; pdf='x^",a,"/(exp(x)-1)'; domain=(",lb,",",ub,")"), "TDR")
        unuran.sample(unr,n)
}

## Powerexponential (Subbotin) distribution ---------------------------------
urpowerexp <- function (n,shape,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("powerexponential(",shape,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Rayleigh distribution ----------------------------------------------------
urrayleigh <- function (n,scale=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("rayleigh(",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Student's t distribution - (replacement for rt) --------------------------
urt <- function (n,df,lb=-Inf,ub=Inf) { 
        unr <- new("unuran", paste("student(",df,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

## Triangular distribution with lower border a, mode m and upper border b ---
urtriang <- function (n,a,m,b,lb=a,ub=b) {
        if (a>=b || m<=a || m>=b)
                stop("Invalid arguments for a,m,b")
        l <- c(a*a, -2*a, 1) / ((a-b)*(a-m))
        r <- c(a*b-a*m+b*m, -2*b, 1) / ((a-b)*(b-m))
        cdfstring <- paste("'(x<=",m,")*(",l[1],"+(",l[2],")*x+(",l[3],")*x*x)+(x>",
                           m,")*(",r[1],"+(",r[2],")*x+(",r[3],")*x*x)';", sep="")
        domainstring <- paste("domain=(",max(lb,a),",",min(ub,b),")", sep="")
        unr <- new("unuran", paste("cont; cdf=",cdfstring,domainstring), "HINV")
        unuran.sample(unr,n)
}

## Weibull distribution - (replacement for rweibull) ------------------------
urweibull <- function (n,shape,scale=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("weibull(",shape,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}


#############################################################################
## Discrete univariate Distributions                                        #
#############################################################################


## Binomial distribution - (replacement for rbinom) -------------------------
urbinom <- function (n,size,prob,lb=0,ub=size) { 
        unr <- new("unuran", paste("binomial(",size,",",prob,"); domain=(",lb,",",ub,")"), "DGT")
        unuran.sample(unr,n)
}

## Geometric distribution - (replacement for rgeom) -------------------------
urgeom <- function (n,prob,lb=0,ub=Inf) {
        if (prob > 0.02) {
                ub  <- min(ub,2000);
                unr <- new("unuran", paste("geometric(",prob,"); domain=(",lb,",",ub,")"), "DGT");
	}
        else {
                unr <- new("unuran", paste("geometric(",prob,"); domain=(",lb,",",ub,")"), "DARI");
        }
        unuran.sample(unr,n)
}
 
## Hypergeometric distribution - (replacement for rhyper) -------------------
urhyper <- function (nn,m,n,k,lb=max(0,k-n),ub=min(k,m)) {
        unr <- new("unuran", paste("hypergeometric(",m+n,",",m,",",k,"); domain=(",lb,",",ub,")"), "DGT")
        unuran.sample(unr,nn)
}

## Logarithmic distribution -------------------------------------------------
urlogarithmic <- function (n,shape,lb=1,ub=Inf) {
        if(shape<0.98) {
                ub  <- min(ub,2000);
                unr <- new("unuran", paste("logarithmic(",shape,"); domain=(",lb,",",ub,")"), "DGT")
        }
        else {
                unr <- new("unuran", paste("logarithmic(",shape,"); domain=(",lb,",",ub,")"), "DARI")
        }
        unuran.sample(unr,n)
}

## Negative binomial distribution - (replacement for rnbinom) ---------------
urnbinom <- function (n,size,prob,lb=0,ub=Inf) {
        if (pnbinom(1000,size,prob,lower.tail=F) < 1.e-10){
                ub  <- min(ub,1000);
                unr <- new("unuran", paste("negativebinomial(",prob,",",size,"); domain=(",lb,",",ub,")"), "DGT")
        }
        else {
                unr <- new("unuran", paste("negativebinomial(",prob,",",size,"); domain=(",lb,",",ub,")"), "DARI")
        }
        unuran.sample(unr,n)
}

## Poisson distribution - (replacement for rpois) ---------------------------
urpois <- function (n,lambda,lb=0,ub=Inf) {
        if (ppois(1000,lambda,lower.tail=F) < 1.e-10) {
                ub <- min(ub,1000);
                unr <- new("unuran", paste("poisson(",lambda,"); domain=(",lb,",",ub,")"), "DGT")
        }
        else {
                unr <- new("unuran", paste("poisson(",lambda,"); domain=(",lb,",",ub,")"), "DARI")
        }
        unuran.sample(unr,n)
}


#############################################################################
## Continuous multivariate Distributions:                                   #
##    Yet net implemented.                                                  #
#############################################################################
