\name{ud}
\docType{function}

\alias{ud}

\title{Density function for "unuran" distribution object}

\description{
  Evaluates the probability density function or
  probability mass function of a continuous and discrete "unuran" 
  distribution object.

  [ EXPERIMENTAL CODE ! ]
}

\usage{
ud(distr, x)
}
\arguments{
  \item{distr}{distribution object. (S4 object of
    class \code{"unuran.cont"} or \code{"unuran.discr"})}
  \item{x}{vector of x values. (numeric)}
}

\details{
  The routine evaluates the density of inside the UNU.RAN distribution
  object.
}

\seealso{
  \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}}.
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
## Create an UNU.RAN distribution object (for standard Gaussian)
distr <- udnorm()

## Evaluate density at some points
ud(distr, 1.5)
ud(distr, -3:3)

}
