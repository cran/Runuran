\name{unuran.cont-class}
\docType{class}

\alias{unuran.cont-class}
\alias{print,unuran.cont-method}
\alias{show,unuran.cont-method}
\alias{initialize,unuran.cont-method}

\title{Class "unuran.cont" for Continuous Distribution}

\description{%
  Class \code{unuran.cont} provides an interface to UNU.RAN objects for
  continuous distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Continuous Distribution Object.
}

\details{
  Create a new instance of a \code{unuran.cont} object using

  \code{new ("unuran.cont", cdf=NULL, pdf=NULL, dpdf=NULL, islog=FALSE,
             lb=NA, ub=NA, mode=NA, center=NA, area=NA, name=NA)}.

  \describe{
    \item{cdf}{cumulative distribution function. (\R function)}
    \item{pdf}{probability density function. (\R function)}
    \item{dpdf}{derivative of the \code{pdf}. (\R function)}
    \item{islog}{whether the given \code{cdf} and \code{pdf} are given
      by their logarithms (the \code{dpdf} is then the derivative of the
      logarithm). (boolean)}
    \item{lb}{lower bound of domain;
      use \code{-Inf} if unbounded from left. (numeric)}
    \item{ub}{upper bound of domain;
      use \code{Inf} if unbounded from right. (numeric)}
    \item{mode}{mode of distribution. (numeric)}
    \item{center}{typical point (\dQuote{center}) of distribution. If
      not given the \code{mode} is used. (numeric)}
    \item{area}{area below \code{pdf}; used for computing normalization
      constants if required. (numeric)}
    \item{name}{name of distribution. (string)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
} 

\examples{
## Create continuous distribution with given logPDF and its derivative
pdf <- function (x) { -0.5*x^2 }
dpdf <- function (x) { -x }
distr <- new("unuran.cont", pdf=pdf, dpdf=dpdf, islog=TRUE, lb=-Inf, ub=Inf)

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
  \code{\link{unuran.cont.new}}, \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
