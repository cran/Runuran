\name{unuran.cmv.new}
\docType{function}

\alias{unuran.cmv.new}

\title{Create a UNU.RAN continuous multivariate distribution object}

\description{
  Create a new UNU.RAN object for a continuous multivariate distribution.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Continuous Multivariate Distribution.
}

\usage{
unuran.cmv.new(dim=1, pdf=NULL, ll=NULL, ur=NULL,
               mode=NULL, center=NULL, name=NA)
}

\arguments{
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

\details{
  Creates an instance of class \code{\linkS4class{unuran.cmv}}.

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
}

\note{
  \code{unuran.cmv.new(\dots)} is an alias for 
  \code{new("unuran.cmv", \dots)}.
}

\seealso{
  \code{\linkS4class{unuran.cmv}}, \code{\link{unuran.new}},
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
## Get a distribution object with given pdf and mode
mvpdf <- function (x) { exp(-sum(x^2)) }
mvd <- unuran.cmv.new(dim=2, pdf=mvpdf, mode=c(0,0))

## Restrict domain to rectangle [0,1]x[0,1] and set
## mode to (0,0)
mvpdf <- function (x) { exp(-sum(x^2)) }
mvd <- unuran.cmv.new(dim=2, pdf=mvpdf, ll=c(0,0), ur=c(1,1), mode=c(0,0))

}

\keyword{distribution}
\keyword{datagen}
