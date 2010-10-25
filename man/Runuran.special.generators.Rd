\name{Runuran.special.generators}

\alias{Runuran.special.generators}
\alias{SpecialGenerator}

\title{
  Generators for distributions based on methods from the UNU.RAN library
}

\description{
  Generators for particular distributions. Their syntax is similar to
  the corresponding \R built-in functions.
}

\details{
  \pkg{Runuran} provides an interface to the UNU.RAN library for
  universal non-uniform random number generators. This is a very
  flexible and powerful collection of sampling routines, where the user
  first has to specify the target distribution and then has to choose an
  appropriate sampling method.
  However, we found that this approach is a little bit confusing for the
  beginner.

  Thus we have prepared easy-to-use sampling functions for
  standard distributions to facilitate the use of the package.
  All these functions share a similar syntax and naming scheme 
  (only \code{u} is prefixed) with their analogous \R
  built-in generating functions (if these exist) but have optional
  domain arguments \code{lb} and \code{ub},
  i.e., these calls also allow to draw samples from truncated
  distributions:

  \code{ur\ldots(n, distribution parameters, lb , ub)}
  
  These functions also show the interested user how we used the more
  powerful functions. We recommend to directly use these more flexible
  functions. Then one has faster marginal generation times and one may
  choose the best generation method for one's application.

  Currently generators for the following distributions are implemented.
  
  %% -- begin: list of distributions --
  Continuous Univariate Distributions (24):

  \tabular{lcl}{ 
    \emph{Function} \tab \tab \emph{Distribution} \cr
    \code{\link{urbeta}} \tab \ldots \tab Beta \cr
    \code{\link{urburr}} \tab \ldots \tab Burr \cr
    \code{\link{urcauchy}} \tab \ldots \tab Cauchy \cr
    \code{\link{urchi}} \tab \ldots \tab Chi \cr
    \code{\link{urchisq}} \tab \ldots \tab Chi-squared \cr
    \code{\link{urexp}} \tab \ldots \tab Exponential \cr
    \code{\link{urextremeI}} \tab \ldots \tab Gumbel (extreme value type I) \cr
    \code{\link{urextremeII}} \tab \ldots \tab Frechet (extreme value type II) \cr
    \code{\link{urf}} \tab \ldots \tab F \cr
    \code{\link{urgamma}} \tab \ldots \tab Gamma \cr
    \code{\link{urgig}} \tab \ldots \tab GIG (generalized inverse Gaussian) \cr
    \code{\link{urhyperbolic}} \tab \ldots \tab Hyperbolic \cr
    \code{\link{urlaplace}} \tab \ldots \tab Laplace \cr
    \code{\link{urlnorm}} \tab \ldots \tab Log-Normal \cr
    \code{\link{urlogis}} \tab \ldots \tab Logistic \cr
    \code{\link{urlomax}} \tab \ldots \tab Lomax \cr
    \code{\link{urnorm}} \tab \ldots \tab Normal (Gaussian) \cr
    \code{\link{urpareto}} \tab \ldots \tab Pareto (of first kind) \cr
    \code{\link{urplanck}} \tab \ldots \tab Planck \cr
    \code{\link{urpowerexp}} \tab \ldots \tab Powerexponential (Subbotin) \cr
    \code{\link{urrayleigh}} \tab \ldots \tab Rayleigh \cr
    \code{\link{urt}} \tab \ldots \tab t (Student) \cr
    \code{\link{urtriang}} \tab \ldots \tab Triangular \cr
    \code{\link{urweibull}} \tab \ldots \tab Weibull \cr
  }

  Discrete Distributions (6):

  \tabular{lcl}{ 
    \emph{Function} \tab \tab \emph{Distribution} \cr
    \code{\link{urbinom}} \tab \ldots \tab Binomial \cr
    \code{\link{urgeom}} \tab \ldots \tab Geometric \cr
    \code{\link{urhyper}} \tab \ldots \tab Hypergeometric \cr
    \code{\link{urlogarithmic}} \tab \ldots \tab Logarithmic \cr
    \code{\link{urnbinom}} \tab \ldots \tab Negative Binomial \cr
    \code{\link{urpois}} \tab \ldots \tab Poisson \cr
  }

  %% -- end: list of distributions --
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

%% \references{}

\seealso{
  \code{\link{Runuran-package}}, \code{\link{Runuran.distributions}}.
}

\examples{
## draw a sample of size 100 from a
## gamma distribution with shape parameter 5
x <- urgamma(n=100, shape=5)

## draw a sample of size 100 from a
## half normal distribution
x <- urnorm(n=100, lb=0, ub=Inf)

}

\keyword{distribution}
\keyword{datagen}
