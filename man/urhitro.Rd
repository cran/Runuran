\name{urhitro}
\docType{function}

\alias{urhitro}

\title{UNU.RAN transformend density rejection generator}

\description{
  UNU.RAN random variate generator for continuous multivariate
  distributions with given probability density function (PDF).
  It is based on the Hit-and-Run algorithm in combinaton with the
  Ratio-of-Uniforms method (HITRO).
}

\usage{
urhitro(n, dim=1, pdf, mode=NULL, center=NULL, ll=NULL, ur=NULL, thinning=1, burnin = 0, ...)
}
\arguments{
  \item{n}{size of required sample.}
  \item{dim}{number of dimensions of the distribution. (integer)}
  \item{pdf}{probability density function. (\R function)}
  \item{mode}{location of the mode -- optional. (numeric vector)}
  \item{center}{point in \dQuote{typical} region of distribution,
    e.g. the approximate location of the mode. If omitted the
    \code{mode} is implicitly used. If the \code{mode} is not given
    either, the origin is used. (numeric vector -- optional)}
  \item{ll,ur}{lower left and upper right vertex of a rectangular
    domain of the \code{pdf}. The domain is only set if both vertices
    are not \code{NULL}. Otherwise, the domain is unbounded by
    default. (numeric vectors)}
  \item{thinning}{thinning factor. (positive integer)}
  \item{burnin}{lenght of burin-in phase. (positive integer)}
  \item{...}{(optional) arguments for \code{pdf}}
}

\details{
  This is a Markov chain sampler that generates a sequence of continuous
  random vectors with given probability density function. This function
  must be provided by \code{pdf}, a function which must return
  non-negative numbers and which need not be normalized (i.e., it can be
  any multiple of a density function).

  The \code{center} is used as starting point of the Hit-and-Run
  algorithm. It is thus important, that the center is contained in the
  (closure of the) domain. 
  Alternatively, one could the \code{mode}. However, this requires its
  exact position whereas \code{center} allows any point in the
  \dQuote{typical} region of the distribution.
  
  The \code{mode} is used a starting point of the Hit-and-Run
  algorithm. It is thus important, that the mode is contained in the
  (closure of the) domain.
  
  Notice that MCMC samplers do not produce a sequence of independent
  variates. The drawn sample follows the target distribution only
  appropriately.
  \strong{Beware: MCMC sampling can be dangerous!}
      
  The algorithm works best with log-concave distributions. Other
  distributions work as well but convergence can be slower.
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  R. Karawatzki, J. Leydold, and K. P\"otzelberger (2005):
  Automatic Markov Chain Monte Carlo Procedures for Sampling from
  Multivariate Distributions.
  Research Report Series / Department of Statistics and Mathematics,
  Nr. 27, December 2005 Department of Statistics and Mathematics,
  Wien, Wirtschaftsuniv., 2005.
  \url{http://epub.wu-wien.ac.at/dyn/virlib/wp/showentry?ID=epub-wu-01_8cb}
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 100 for a 
## Gaussian distribution
mvpdf <- function (x) { exp(-sum(x^2)) }
x <- urhitro(20, dim=2, pdf=mvpdf)
x

## use:
##  mode at (0,0)
##  thinning factor 3
##    (only every 3rd vector in the sequence is returned)
##  burn-in of length 100
##    (the first 1000 vectors in the sequence are discarded)
x <- urhitro(20, dim=2, pdf=mvpdf, mode=c(0,0), thinning=3, burnin=100)
x

## Create a sample of size 100 for a 
## Gaussian distribution restricted to first quadrant
mvpdf <- function (x) { exp(-sum(x^2)) }
x <- urhitro(20, dim=2, pdf=mvpdf, ll=c(0,0), ur=c(Inf,Inf))
x

## Create a sample of size 100 for a 
## Gaussian distribution restricted to the rectangle [1,2]x[1,2]
## (don't forget to provide a starting point using 'center')

mvpdf <- function (x) { exp(-sum(x^2)) }
x <- urhitro(20, dim=2, pdf=mvpdf, center=c(1.1,1.1), ll=c(1,1), ur=c(2,2))
x

}

\keyword{distribution}
\keyword{datagen}

