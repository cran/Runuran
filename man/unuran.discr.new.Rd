\name{unuran.discr.new}
\alias{unuran.discr.new}

\title{Create a UNU.RAN discrete univariate distribution object}

\description{
  Create a new UNU.RAN object for a discrete univariate distribution.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
  
  [Advanced] -- Discrete Distribution.
}

\usage{
unuran.discr.new(pv=NULL, pmf=NULL, lb=NA, ub=NA,
                 mode=NA, sum=NA, name=NA)
}

\arguments{
    \item{pv}{probability vector. (numeric vector)}
    \item{pmf}{probability mass function. (\R function)}
    \item{mode}{mode of distribution. (numeric, integer)}
    \item{lb}{lower bound of domain;
      use \code{-Inf} if unbounded from left. (numeric, integer)}
    \item{ub}{upper bound of domain;
      use \code{Inf} if unbounded from right;
      when \code{pmf} is not given, the default \code{ub=Inf} is
      used. (numeric, integer)}
    \item{sum}{sum over \code{pv} / \code{pmf}; used for computing
      normalization constants if required. (numeric)}
    \item{name}{name of distribution. (string)}
}

\details{
  Creates an instance of class \code{\linkS4class{unuran.discr}}.
  For more details see also \code{\link{unuran.new}}.

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
}

\note{
  \code{unuran.discr.new(\dots)} is an alias for 
  \code{new("unuran.discr", \dots)}.
}

\seealso{
  \code{\linkS4class{unuran.discr}}, \code{\link{unuran.new}},
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
## Create a distribution object with given PV and mode
mypv <- dbinom(0:100,100,0.3)
distr <- new("unuran.discr", pv=mypv, lb=0, mode=30)

## Create discrete distribution with given probability vector
## (the PV need not be normalized)
pv <- c(1.,2.,1.5,0.,3.,1.2)     
dpv <- new("unuran.discr", pv=pv, lb=1)

## Create discrete distribution with given PMF
pmf <- function(x) dbinom(x,100,0.3)
dpmf <- new("unuran.discr", pmf=pmf, lb=0, ub=100)

}

\keyword{distribution}
\keyword{datagen}

