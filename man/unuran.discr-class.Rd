\name{unuran.discr-class}
\docType{class}

\alias{unuran.discr-class}
\alias{print,unuran.discr-method}
\alias{show,unuran.discr-method}
\alias{initialize,unuran.discr-method}

\title{Class "unuran.discr" for Discrete Distribution}

\description{%
  Class \code{unuran.discr} provides an interface to UNU.RAN objects for
  discrete distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Discrete Distribution Object.
}

\details{
  Create a new instance of a \code{unuran.discr} object using

  \code{new ("unuran.discr", pv=NULL, pmf=NULL, lb=NA, ub=NA,
             mode=NA, sum=NA, name=NA)}.

  \describe{
    \item{pv}{probability vector. (numeric vector)}
    \item{pmf}{probability mass function. (\R function)}
    \item{lb}{lower bound of domain;
      use \code{-Inf} if unbounded from left. (numeric, integer)}
    \item{ub}{upper bound of domain;
      use \code{Inf} if unbounded from right;
      when \code{pmf} is not given, the default \code{ub=Inf} is
      used. (numeric, integer)}
    \item{mode}{mode of distribution. (numeric, integer)}
    \item{sum}{sum over \code{pv} / \code{pmf}; used for computing
      normalization constants if required. (numeric)}
    \item{name}{name of distribution. (string)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
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
  \code{\link{unuran.discr.new}}, \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\examples{
## Create discrete distribution with given probability vector
## (the PV need not be normalized)
pv <- c(1.,2.,1.5,0.,3.,1.2)     
dpv <- new("unuran.discr", pv=pv, lb=1)

## Create discrete distribution with given PMF
pmf <- function(x) dbinom(x,100,0.3)
dpmf <- new("unuran.discr",pmf=pmf,lb=0,ub=100)

}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
