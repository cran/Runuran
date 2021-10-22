\name{unuran.cmv-class}
\docType{class}

\alias{unuran.cmv-class}
\alias{print,unuran.cmv-method}
\alias{show,unuran.cmv-method}
\alias{initialize,unuran.cmv-method}

\title{Class "unuran.cmv" for Continuous Multivariate Distribution}

\description{%
  Class \code{unuran.cmv} provides an interface to UNU.RAN objects
  for continuous multivariate distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Continuous Multivariate Distribution Object.
}

\details{
  Create a new instance of a \code{unuran.cmv} object using

  \code{new ("unuran.cmv", dim=1, pdf=NULL, ll=NULL, ur=NULL,
             mode=NULL, center=NULL, name=NA)}.

  \describe{
    \item{dim}{number of dimensions of the distribution. (integer)}
    \item{pdf}{probability density function. (\R function)}
    \item{ll,ur}{lower left and upper right vertex of a rectangular
      domain of the \code{pdf}. The domain is only set if both vertices
      are not \code{NULL}. Otherwise, the domain is unbounded by default.
      (numeric vectors)}
    \item{mode}{location of the mode. (numeric vector -- optional)}
    \item{center}{point in \dQuote{typical} region of distribution,
      e.g. the approximate location of the mode. It is used by several
      methods to locate the main region of the distribution.
      If omitted the \code{mode} is implicitly used. If the \code{mode}
      is not given either, the origin is used. (numeric vector --
      optional)}
    \item{name}{name of distribution. (string)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used. It is important, that the mode is contained in the (closure of
  the) domain.
} 

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\references{
  J. Leydold and W. H\"ormann (2000-2007):
  UNU.RAN User Manual, see
  \url{https://statmath.wu.ac.at/unuran/}.
}

\seealso{
  \code{\link{unuran.cmv.new}}, \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\examples{
## Create distribution with given PDF
mvpdf <- function (x) { exp(-sum(x^2)) }
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf)

## Restrict domain to rectangle [0,1]x[0,1] and set
## mode to (0,0)
mvpdf <- function (x) { exp(-sum(x^2)) }
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf, ll=c(0,0), ur=c(1,1), mode=c(0,0))

}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
