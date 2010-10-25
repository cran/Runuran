\name{ud}
\alias{ud}

\title{Density function for "unuran" object}

\description{
  Evaluates the probability density function (PDF) or
  probability mass function (PMF) for a "unuran" 
  object for a continuous and discrete distribution,
  respectively.
}

\usage{
ud(obj, x)
}
\arguments{
  \item{obj}{one of
    \itemize{
      \item a distribution object of class \code{"unuran.cont"}
      that contains the PDF, or
      \item a distribution object of class \code{"unuran.discr"}
      that contains the PMF, or
      \item a generator object (class \code{"unuran"}) that
      contains the PDF and PMF, resp.
    }
  }
  \item{x}{vector of x values. (numeric)}
}

\details{
  The routine evaluates the probability density function of a
  distribution stored in a UNU.RAN distribution object or UNU.RAN
  generator object. If the PDF is not available in the object, then
  the computation is aborted.
}

\note{
  The generator object must not be packed (see
  \code{\link{unuran.packed}}).
}

\seealso{
  \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}},
  \code{\linkS4class{unuran}}.
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
## and evaluate density for some points
distr <- udnorm()
ud(distr, 1.5)
ud(distr, -3:3)

## Create an UNU.RAN generator object (for standard Gaussian)
## and evaluate density of underyling distribution
gen <- tdrd.new(udnorm())
ud(gen, 1.5)
ud(gen, -3:3)

}
