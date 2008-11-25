\name{unuran.distr-class}
\docType{class}

\alias{unuran.distr-class}
\alias{print,unuran.distr-method}
\alias{show,unuran.distr-method}
\alias{initialize,unuran.distr-method}

\title{Virtual class "unuran.distr"}

\description{%
  The virtual class \code{unuran.distr} provides an interface to UNU.RAN
  objects for distributions.

  The following classes extend this class:
  \describe{
    \item{class \code{unuran.cont}}{
      for univariate continuous distributions,
      see \code{\linkS4class{unuran.cont}}.
    }
    \item{class \code{unuran.discr}}{
      for discrete distributions,
      see \code{\linkS4class{unuran.discr}}.
    }
    \item{class \code{unuran.cmv}}{
      for multivariate continuous distributions,
      see \code{\linkS4class{unuran.cmv}}.
    }
  }
  
  [Advanced] -- Distribution Object.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\references{
  J. Leydold and W. H\"ormann (2000-2007):
  UNU.RAN User Manual, see
  \url{http://statmath.wu-wien.ac.at/unuran/}.
}

\seealso{
  \code{\linkS4class{unuran.cont}}, \code{\linkS4class{unuran.discr}},
  \code{\linkS4class{unuran.cmv}};
  \code{\linkS4class{unuran}}.
}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
