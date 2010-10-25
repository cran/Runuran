\name{unuran.new}
\alias{unuran.new}

\title{Create a UNU.RAN object}

\description{
  Create a new \code{unuran} object in package \pkg{Runuran} that can be used
  for sampling from the specified distribution.
  The function \code{ur} can then be used to draw a random
  sample.

  [Advanced] -- Create generator object.
}

\usage{
unuran.new(distr,method="auto")
}

\arguments{
  \item{distr}{a string or an S4 class describing the distribution.}
  \item{method}{a string describing the random variate generation
    method.}
}

\details{
  This function creates an instance of S4 class
  \code{\linkS4class{unuran}} which contains a generator for the target
  distribution. This distribution has to be provided as an instance of
  S4 class \code{\linkS4class{unuran.distr}}. Depending on the type of
  distribution such an instance can be created by 
  \describe{
    \item{\code{\link{unuran.cont.new}}}{
      for univariate continuous distributions,
    }
    \item{\code{\link{unuran.discr.new}}}{
      for discrete distributions, and
    }
    \item{\code{\link{unuran.cmv.new}}}{
      for multivariate continuous distributions.
    }
  }

  The generation can be chosen by passing \code{method} to the UNU.RAN
  String API. The default method, \code{"auto"} tries to find an
  appropriate method for the given distribution. However, this method is
  experimental and is yet not very powerfull.
  
  Once a \code{unuran} object has been created it can be used to draw random
  samples from the target distribution using \code{\link{ur}}.
}

\seealso{
  See \code{\linkS4class{unuran}} for the UNU.RAN class of generators.
  See \code{\link{unuran.details}} for printing details about the
  generator object, and \code{\link{ur}} and \code{\link{uq}} for
  sampling and quantile function, respectively.

  For distribution objects see \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}}, and \code{\linkS4class{unuran.cmv}}.
  
  \code{\link{runif}}, \code{\link{.Random.seed}} about random number
  generation in \R.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Use method 'TDR' (Transformed Density Rejection) to 
## draw a sample of size 10 from a hyperbolic distribution with PDF
##   f(x) = const * exp(-sqrt(1+x^2)) 
## restricted to domain [-1,2].

## We first have to define functions that return the log-density and
## its derivative, respectively. (We also could use the density itself.)
lf  <- function (x) { -sqrt(1+x^2) }
dlf <- function (x) { -x/sqrt(1+x^2) }

## Next create the continuous distribution object.
d <- unuran.cont.new(pdf=lf,dpdf=dlf,islog=TRUE,lb=-1,ub=2)

## Create 'unuran' object. We choose method 'TDR' with 
## immediate acceptance (IA) and parameter c=0.
gen <- unuran.new(distr=d, method="tdr; variant_ia; c=0")

## Now we can use this object to draw the sample.
## (Of course we can repeat this step as often as required.)
ur(gen,10)

## Here is some information about our generator object.
unuran.details(gen)
}

\keyword{distribution}
\keyword{datagen}

