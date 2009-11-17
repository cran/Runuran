\name{vnrou.new}
\docType{function}

\alias{vnrou.new}

\title{UNU.RAN generator based on Multivariate Naive Ratio-Of-Uniforms method (VNROU)}

\description{
  UNU.RAN random variate generator for continuous multivariate
  distributions with given probability density function (PDF).
  It is based on the Multivariate Naive Ratio-Of-Uniforms method
  (\sQuote{VNROU}).
  
  [Universal] -- Rejection Method.
}

\usage{
vnrou.new(dim=1, pdf, ll=NULL, ur=NULL, mode=NULL, center=NULL, \dots)
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
  \item{\dots}{(optional) arguments for \code{pdf}}
}

\details{
  This function creates a \code{unuran} object based on the naive
  ratio-of-uniforms method (\sQuote{VNROU}), i.e., a bounding rectangle
  for the acceptance region is estimated and use for sampling proposal
  points.
  It can be used to draw samples of a continuous random vector with
  given probability density function using \code{\link{ur}}.

  The algorithm works with unimodal distributions provided that the
  tails are not too \dQuote{high} in every direction.
  
  The density must be provided by a function \code{pdf} which must
  return non-negative numbers and which need not be normalized (i.e., it
  can be any multiple of a density function).
  
  The \code{center} is used as starting point for computing the bounding
  rectangle. Alternatively, one also could provide the location the
  \code{mode}. However, this requires its exact position whereas
  \code{center} allows any point in the \dQuote{typical} region of the
  distribution.
  
  The setup can be accelerated when the \code{mode} is given.
}

\seealso{
  \code{\link{ur}}, \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  Section 11.1.4.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gaussian distribution
mvpdf <- function (x) { exp(-sum(x^2)) }
gen <- vnrou.new(dim=2, pdf=mvpdf)
x <- ur(gen,100)

## Use mode of Gaussian distribution to accelerate set-up.
mvpdf <- function (x) { exp(-sum(x^2)) }
gen <- vnrou.new(dim=2, pdf=mvpdf, mode=c(0,0))
x <- ur(gen,100)

## Gaussian distribution restricted to the rectangle [1,2]x[1,2]
##  (don't forget to provide a point inside domain using 'center')
mvpdf <- function (x) { exp(-sum(x^2)) }
gen <- vnrou.new(dim=2, pdf=mvpdf, ll=c(1,1), ur=c(2,2), center=c(1.5,1.5))
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
