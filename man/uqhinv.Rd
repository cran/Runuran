\name{uqhinv}
\docType{function}

\alias{uqhinv}

\title{Quantile function using UNU.RAN method HINV}

\description{
  Evaluates quantile of distribution approximately using
  UNU.RAN method \sQuote{HINV}.
  
  [ \strong{Deprecated}.
  Use \code{\link{uq}} instead. ]
}

\usage{
uqhinv(unr, U)
}
\arguments{
  \item{unr}{a \code{unuran} object using method \sQuote{HINV}.}
  \item{U}{vector of probabilities.}
}

\details{
  This routine evaluates the quantiles (inverser CDF) for a given
  (vector of) probabilities approximately.
  It requires a \code{unuran} object that uses method \sQuote{HINV} and thus
  contains the information about the target distribution.

  This method is \strong{deprecated}. Please use 
  \code{\link{uq}} instead.
}

\seealso{
  \code{\link{uq}}.
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

%%\examples{}

\keyword{distribution}
\keyword{datagen}

