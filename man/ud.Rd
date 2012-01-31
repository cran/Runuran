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
ud(obj, x, islog = FALSE)
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
  \item{islog}{if \code{TRUE}, the log-density is returnd. (boolean)}
}

\details{
  The routine evaluates the probability density function of a
  distribution stored in a UNU.RAN distribution object or UNU.RAN
  generator object. If \code{islog} is \code{TRUE}, then the logarithm
  of the density is returned.

  If the PDF (or its respective logarithm) is not available in the
  object, then \code{NA} is returned and a warning is thrown.

  Note: when the log-density is not given explicitly (by setting
  \code{islog=TRUE} in the corresponding routing like 
  \code{\link{unuran.cont.new}} or in an \pkg{Runuran} built-in
  distribution), then \code{NA} is returned even if the density is
  given.

  \strong{Important:}
  Routine \code{ud} just evaluates the density function that is stored
  in \code{obj}. It ignores the boundaries of the domain of the
  distribution, i.e., it does not return \code{0} outside the domain
  unless the implementation of the PDF handles this case correctly.
  This behavior is in particular important when \pkg{Runuran}
  built-in distributions are truncated by explicitly setting the domain
  boundaries.
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
