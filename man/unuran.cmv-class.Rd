\name{unuran.cmv-class}
\docType{class}

\alias{unuran.cmv-class}
\alias{print,unuran.cmv-method}
\alias{show,unuran.cmv-method}
\alias{initialize,unuran.cmv-method}

\title{Experimental class "unuran.cmv"}

\description{%
  The class "unuran.cmv" provides am interface to UNU.RAN objects for
  continuous multivariate distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
}

\details{
  Create a new instance of an "unuran.cmv" object using

  \code{new ("unuran.cmv", dim=1, pdf=NULL, mode=NULL)}.

  \describe{
    \item{dim}{number of dimensions of the distribution. (integer)}
    \item{pdf}{probability density function. (\R function)}
    \item{mode}{location of the mode -- optional. (numeric vector)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
} 

\examples{
## Create discrete distribution with given probability density function
mvpdf <- function (x) { exp(-sum(x^2)) }
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf)
## Make generator (using method HITRO)
mvunr <- unuran.new(mvdist, "hitro")
## Draw sample of size 100
x <- unuran.sample(mvunr, 100)
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
  \code{\linkS4class{unuran}}, \code{\link{unuran.new}}.
}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
