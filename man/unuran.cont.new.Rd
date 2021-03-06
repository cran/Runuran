\name{unuran.cont.new}
\alias{unuran.cont.new}

\title{Create a UNU.RAN continuous univariate distribution object}

\description{
  Create a new UNU.RAN object for a continuous univariate distribution.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Continuous Distribution.
}

\usage{
unuran.cont.new( cdf=NULL, pdf=NULL, dpdf=NULL, islog=FALSE,
                 lb=NA, ub=NA, mode=NA, center=NA, area=NA, name=NA)
}

\arguments{
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

\details{
  Creates an instance of class \code{\linkS4class{unuran.cont}}.

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
}

\note{
  \code{unuran.cont.new(\dots)} is an alias for 
  \code{new("unuran.cont", \dots)}.
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
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Get a distribution object with given pdf, domain and mode
mypdf <- function (x) { exp(-x) }
distr <- unuran.cont.new(pdf=mypdf, islog=FALSE, lb=0, ub=Inf, mode=0)

## This object can now be used to create an generator object.
## 1. select a method using a Runuran function:
gen <- pinvd.new(distr, uresolution=1e-12)

## 2. directly use the UNU.RAN string API
gen <- unuran.new(distr, method="pinv; u_resolution=1e-12")

}

\keyword{distribution}
\keyword{datagen}
