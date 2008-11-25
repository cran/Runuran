\name{unuran.cont.new}
\docType{function}

\alias{unuran.cont.new}

\title{Create a UNU.RAN continuous univariate distribution object}

\description{
  Create a new UNU.RAN object for a continuous univariate distribution.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Continuous Distribution.
}

\usage{
unuran.cont.new( cdf=NULL, pdf=NULL, dpdf=NULL, islog=FALSE, mode=NA,
                 center=NA, lb=-Inf, ub=Inf, area=NA, name=NA)
}
\arguments{
  \item{cdf}{cumulative distribution function. (\R function)}
  \item{pdf}{probability density function. (\R function)}
  \item{dpdf}{derivative of the \code{pdf}. (\R function)}
  \item{islog}{whether the given \code{cdf} and \code{pdf} are given
    by their logarithms (the \code{dpdf} is then the derivative of the
    logarithm). (boolean)}
  \item{mode}{mode of distribution. (numeric)}
  \item{center}{typical point (\dQuote{center}) of distribution. If
    not given the \code{mode} is used. (numeric)}
  \item{lb, ub}{lower and upper bound of domain. (numeric)}
  \item{area}{area below \code{pdf}; used for computing normalization
    constants if required. (numeric)}
  \item{name}{name of distribution. (string)}
}

\details{
  Creates an instance of class \code{\linkS4class{unuran.cont}}.

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
}

\note{
  \code{unuran.cont.new(...)} is an alias for 
  \code{new("unuran.cont", ...)}.
}

\seealso{
  \code{\linkS4class{unuran.cont}}, \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Get a distribution object with given pdf, domain and mode
mypdf <- function (x) { exp(-x) }
distr <- unuran.cont.new(pdf=mypdf, islog=FALSE, mode=0, lb=0, ub=Inf)

}

\keyword{distribution}
\keyword{datagen}

