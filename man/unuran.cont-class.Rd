\name{unuran.cont-class}
\docType{class}

\alias{unuran.cont-class}
\alias{print,unuran.cont-method}
\alias{show,unuran.cont-method}
\alias{initialize,unuran.cont-method}

\title{Experimental class "unuran.cont"}

\description{%
  The class "unuran.cont" provides am interface to UNU.RAN objects for
  continuous distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
}

\details{
  Create a new instance of an "unuran.cont" object using

  \code{new ("unuran.cont", cdf=NULL, pdf=NULL, dpdf=NULL, islog=TRUE, lb=-Inf, ub=Inf)}.

  \describe{
    \item{cdf}{cumulative distribution function. (\R function)}
    \item{pdf}{probability density function. (\R function)}
    \item{dpdf}{derivative of the \code{pdf}. (\R function)}
    \item{islog}{whether the given \code{cdf} and \code{pdf} are given
      by their logarithms (the \code{dpdf} is then the derivative of the
      logarithm). (boolean)}
    \item{lb, ub}{lower and upper bound of domain. (numeric)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
} 

\examples{
## Create discrete distribution with given probability density function
pdf <- function (x) { -0.5*x^2 }
dpdf <- function (x) { -x }
distpdf <- new("unuran.cont", pdf=pdf, dpdf=dpdf)
## Make generator (using method TDR)
unr <- unuran.new(distpdf, "tdr")
## Draw sample of size 100
x <- unuran.sample(unr, 100)

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
