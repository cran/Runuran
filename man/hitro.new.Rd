\name{hitro.new}
\alias{hitro.new}

\title{UNU.RAN generator based on Hit-and-Run sampler (HITRO)}

\description{
  UNU.RAN random variate generator for continuous multivariate
  distributions with given probability density function (PDF).
  It is based on the Hit-and-Run algorithm in combinaton with the
  Ratio-of-Uniforms method (\sQuote{HITRO}).

  [Universal] -- MCMC (Markov chain sampler).
}

\usage{
hitro.new(dim=1, pdf, ll=NULL, ur=NULL, mode=NULL, center=NULL,
          thinning=1, burnin=0, \dots)
}
\arguments{
  \item{dim}{number of dimensions of the distribution. (integer)}
  \item{pdf}{probability density function. (\R function)}
  \item{ll,ur}{lower left and upper right vertex of a rectangular
    domain of the \code{pdf}. The domain is only set if both vertices
    are not \code{NULL}. Otherwise, the domain is unbounded by
    default. (numeric vectors)}
  \item{mode}{location of the mode. (numeric vector)}
  \item{center}{point in \dQuote{typical} region of distribution,
    e.g. the approximate location of the mode. If omitted the
    \code{mode} is used. If the \code{mode} is not given
    either, the origin is used. (numeric vector)}
  \item{thinning}{thinning factor. (positive integer)}
  \item{burnin}{length of burnin-in phase. (positive integer)}
  \item{\dots}{(optional) arguments for \code{pdf}}
}

\details{
  \strong{Beware: MCMC sampling can be dangerous!}
  
  This function creates a \code{unuran} object based on the
  Hit-and-Run algorithm in combinaton with the Ratio-of-Uniforms method
  (\sQuote{HITRO}).
  It can be used to draw samples of a continuous random vector with
  given probability density function using \code{\link{ur}}.

  The algorithm works best with log-concave distributions. Other
  distributions work as well but convergence can be slower.
  
  The density must be provided by a function \code{pdf} which must
  return non-negative numbers and but need not be normalized (i.e., it
  can be any multiple of a density function).
  
  The \code{center} is used as starting point of the Hit-and-Run
  algorithm. It is thus important, that the \code{center} is contained
  in the (interior of the) domain. 
  Alternatively, one could provide the location of the
  \code{mode}. However, this requires its exact position whereas
  \code{center} allows any point in the \dQuote{typical} region of the
  distribution.

  If the \code{mode} is given, then it is used to obtain an upper bound
  on the \code{pdf} and thus its location should be given sufficiently
  accurate.
  
  The \sQuote{HITRO} algorithm is a MCMC samplers and thus it does not produce a
  sequence of independent variates. The drawn sample follows the target
  distribution only approximately.
  The dependence between consecutive vectors can be decreased when
  only a subsequence is returned (and the other elements are erased).
  This is called \dQuote{thinning} of the Markov chain and can be
  controlled by the \code{thinning} factor. A thinning factor \eqn{k}
  means that only every \eqn{k}-th element is returned.

  Markov chains also depend on the chosen starting point (i.e., the
  \code{center} in this implementation of the algorithm).
  This dependence can be decreased by erasing the first
  part of the chain. This is called the \dQuote{burn-in} of the Markov
  chain and its length is controlled by the argument \code{burnin}.
}

\seealso{
  \code{\link{ur}}, \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  R. Karawatzki, J. Leydold, and K. P\"otzelberger (2005):
  Automatic Markov Chain Monte Carlo Procedures for Sampling from
  Multivariate Distributions.
  Research Report Series / Department of Statistics and Mathematics,
  Nr. 27, December 2005 Department of Statistics and Mathematics,
  Wien, Wirtschaftsuniv., 2005.
  \url{https://epub.wu.ac.at/id/eprint/1400}
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a 
## Gaussian distribution
mvpdf <- function (x) { exp(-sum(x^2)) }
gen <- hitro.new(dim=2, pdf=mvpdf)
x <- ur(gen,100)

## Use mode of Gaussian distribution.
## Reduce auto-correlation by thinning and burn-in.
##  mode at (0,0)
##  thinning factor 3
##    (only every 3rd vector in the sequence is returned)
##  burn-in of length 1000
##    (the first 100 vectors in the sequence are discarded)
mvpdf <- function (x) { exp(-sum(x^2)) }
gen <- hitro.new(dim=2, pdf=mvpdf, mode=c(0,0), thinning=3, burnin=1000)
x <- ur(gen,100)

## Gaussian distribution restricted to the rectangle [1,2]x[1,2]
##  (don't forget to provide a starting point using 'center')
mvpdf <- function (x) { exp(-sum(x^2)) }
gen <- hitro.new(dim=2, pdf=mvpdf, center=c(1.1,1.1), ll=c(1,1), ur=c(2,2))
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
