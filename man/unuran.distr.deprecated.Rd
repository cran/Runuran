\name{unuran.cont}
\docType{function}

\alias{unuran.cmv}
\alias{unuran.cont}
\alias{unuran.discr}

\title{Create a UNU.RAN distribution objects}

\description{
  Create a new UNU.RAN object for distribution.
  
  [ \strong{Deprecated}! ]
}

\usage{
unuran.cmv(dim=1, pdf=NULL, mode=NULL, center=NULL, ll=NULL, ur=NULL)

unuran.cont(cdf=NULL, pdf=NULL, dpdf=NULL, islog=TRUE, mode=NA,
            center=NA, lb=-Inf, ub=Inf, area=NA)

unuran.discr(pv=NULL, pmf=NULL, mode=NA, lb=0, ub=Inf, sum=NA)
}

\arguments{
  \item{dim}{number of dimensions of the distribution. (integer)}
  \item{cdf}{cumulative distribution function. (\R function)}
  \item{pdf}{probability density function. (\R function)}
  \item{dpdf}{derivative of the \code{pdf}. (\R function)}
  \item{pv}{probability vector. (numeric vector)}
  \item{pmf}{probability mass function. (\R function)}
  \item{ll,ur}{lower left and upper right vertex of a rectangular
    domain of the \code{pdf}. The domain is only set if both vertices
    are not \code{NULL}. Otherwise, the domain is unbounded by default.
    (numeric vectors)}
  \item{lb, ub}{lower and upper bound of domain. (numeric)}
  \item{islog}{whether the given \code{cdf} and \code{pdf} are given
    by their logarithms (the \code{dpdf} is then the derivative of the
    logarithm). (boolean)}
  \item{mode}{mode of distribution. (numeric)}
  \item{center}{typical point (\dQuote{center}) of distribution. If
    not given the \code{mode} is used. (numeric)}
  \item{area}{area below \code{pdf}; used for computing normalization
    constants if required. (numeric)}
    \item{sum}{sum over \code{pv} / \code{pmf}; used for computing
      normalization constants if required. (numeric)}
}

\details{
  Creates an instance of class \code{\linkS4class{unuran.distr}}.
  
  These functions are deprecated. Use respective functions
  \code{\link{unuran.cmv.new}}, 
  \code{\link{unuran.cont.new}}, and
  \code{\link{unuran.discr.new}}. 
}

\seealso{
  \code{\linkS4class{unuran.cmv}},
  \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}}, 
  \code{\linkS4class{unuran.distr}}.
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

\keyword{distribution}
\keyword{datagen}

