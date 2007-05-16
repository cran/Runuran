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
urhitro(n, dim=1, pdf, mode=NULL, thinning=1, burnin = 0, ...)
}
\arguments{
  \item{n}{size of required sample.}
  \item{dim}{number of dimensions of the distribution. (integer)}
  \item{pdf}{probability density function. (\R function)}
  \item{mode}{location of the mode -- optional. (numeric vector)}
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
}

\keyword{distribution}
\keyword{datagen}

