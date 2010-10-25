\name{up}
\alias{up}

\title{Distribution function for "unuran" object}

\description{
  Evaluates the (approximate) cumulative distribution function (CDF) 
  of a "unuran" object for a continuous or discrete distribution.
}

\usage{
up(obj, x)
}
\arguments{
  \item{obj}{one of
    \itemize{
      \item a distribution object of class \code{"unuran.cont"}
      that contains the CDF, or
      \item a distribution object of class \code{"unuran.discr"}
      that contains the CDF, or
      \item a generator object (class \code{"unuran"}) that
      contains a CDF or implements method \sQuote{PINV}.
    }
  }
  \item{x}{vector of x values. (numeric)}
}

\details{
  The routine evaluates the cumulative distribution function of a
  distribution stored in a UNU.RAN distribution object or UNU.RAN
  generator object. If the CDF is not available in the object, then
  \code{obj} must be a generator object that implements method
  \sQuote{PINV}. In this case an approximate value for the CDF is
  returned. The approximation error is about one tenth of the requested
  \code{uresolution} for method \sQuote{PINV}. 
  
  If a generator object is given that implements method PINV but also
  contain the CDF, the value of the (exact) CDF is returned. 
}

\note{
  The generator object must not be packed (see
  \code{\link{unuran.packed}}).
}

\seealso{
  \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}},
  \code{\linkS4class{unuran}},
  \code{\link{pinv.new}}.
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
up(distr, 1.5)
up(distr, -3:3)

## Create an UNU.RAN generator object (for standard Gaussian)
## and evaluate density of underyling distribution
unr <- tdrd.new(udnorm())
up(unr, 1.5)
up(unr, -3:3)

## Create an UNU.RAN generator object that does not contain
## the CDF but implements method PINV.
unr <- pinv.new(pdf=function(x){exp(-x)}, lb=0,ub=Inf)
up(unr, 0:5)

}
