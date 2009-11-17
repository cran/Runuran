\name{Runuran.distributions}
\docType{package}

\alias{Runuran.distributions}

\title{UNU.RAN distribution objects}

\description{
  Create objects for particular distributions suitable for using with
  generation methods from the UNU.RAN library.
}

\details{
  \pkg{Runuran} provides an interface to the UNU.RAN library for
  universal non-uniform random number generators. This is a very
  flexible and powerful collection of sampling routines, where the user
  first has to specify the target distribution and then has to choose an
  appropriate sampling method.

  Creating an object for a particular distribution can be a bit
  tedious especially if the target distribution has a more complex
  density function. Thus we have compiled a set of functions that
  provides ready-to-use distribution objects.
  Moreover, using these object often results in faster setup time than
  objects created with pure \R code.

  These functions share a similar syntax and naming scheme 
  (only \code{ud} is prefixed) with analogous \R built-in functions that
  provide density, distribution function and quantile:

  \code{ud\ldots(distribution parameters, lb , ub)}
  
  Currently generators for the following distributions are implemented.
  
  %% -- begin: list of distributions --
  Continuous Univariate Distributions (24):

  \tabular{lcl}{ 
    \emph{Function} \tab \tab \emph{Distribution} \cr
    \code{\link{udbeta}} \tab \ldots \tab Beta \cr
    \code{\link{udcauchy}} \tab \ldots \tab Cauchy \cr
    \code{\link{udchi}} \tab \ldots \tab Chi \cr
    \code{\link{udchisq}} \tab \ldots \tab Chi-squared \cr
    \code{\link{udexp}} \tab \ldots \tab Exponential \cr
    \code{\link{udf}} \tab \ldots \tab F \cr
    \code{\link{udfrechet}} \tab \ldots \tab Frechet (Extreme value type II) \cr
    \code{\link{udgamma}} \tab \ldots \tab Gamma \cr
    \code{\link{udghyp}} \tab \ldots \tab Generalized Hyperbolic \cr
    \code{\link{udgig}} \tab \ldots \tab Generalized Inverse Gaussian \cr
    \code{\link{udgumbel}} \tab \ldots \tab Gumbel (Extreme value type I) \cr
    \code{\link{udhyperbolic}} \tab \ldots \tab Hyperbolic \cr
    \code{\link{udig}} \tab \ldots \tab Inverse Gaussian (Wald) \cr
    \code{\link{udlaplace}} \tab \ldots \tab Laplace (double exponential) \cr
    \code{\link{udlnorm}} \tab \ldots \tab Log Normal \cr
    \code{\link{udlogis}} \tab \ldots \tab Logistic \cr
    \code{\link{udlomax}} \tab \ldots \tab Lomax (Pareto of second kind) \cr
    \code{\link{udnorm}} \tab \ldots \tab Normal (Gaussian) \cr
    \code{\link{udpareto}} \tab \ldots \tab Pareto (of first kind) \cr
    \code{\link{udpowerexp}} \tab \ldots \tab Powerexponential (Subbotin) \cr
    \code{\link{udrayleigh}} \tab \ldots \tab Rayleigh \cr
    \code{\link{udslash}} \tab \ldots \tab Slash \cr
    \code{\link{udt}} \tab \ldots \tab t (Student) \cr
    \code{\link{udweibull}} \tab \ldots \tab Weibull (Extreme value type III) \cr
  }

  Discrete Distributions (6):

  \tabular{lcl}{ 
    \emph{Function} \tab \tab \emph{Distribution} \cr
    \code{\link{udbinom}} \tab \ldots \tab Binomial \cr
    \code{\link{udgeom}} \tab \ldots \tab Geometric \cr
    \code{\link{udhyper}} \tab \ldots \tab Hypergeometric \cr
    \code{\link{udlogarithmic}} \tab \ldots \tab Logarithmic \cr
    \code{\link{udnbinom}} \tab \ldots \tab Negative Binomial \cr
    \code{\link{udpois}} \tab \ldots \tab Poisson \cr
  }

  %% -- end: list of distributions --
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

%% \references{}

\seealso{
  \code{\link{Runuran-package}}.
}

\examples{
## Create an object for a gamma distribution with shape parameter 5.
distr <- udgamma(shape=5)
## Create the UNU.RAN generator object. use method PINV (inversion).
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)
## Compute some quantiles for Monte Carlo methods
x <- uq(gen, (1:9)/10)

## Analogous for half normal distribution
distr <- udnorm(lb=0, ub=Inf)
gen <- pinvd.new(distr)
x <- ur(gen,100)
x <- uq(gen, (1:9)/10)

## Analogous for a generalized hyperbolic distribution
distr <- udghyp(lambda=-1.0024, alpha=39.6, beta=4.14, delta=0.0118, mu=-0.000158)
gen <- pinvd.new(distr)
x <- ur(gen,100)
x <- uq(gen, (1:9)/10)

}

\keyword{distribution}
\keyword{datagen}
