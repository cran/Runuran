\name{unuran.is.inversion}
\alias{unuran.is.inversion}

\title{Test whether a "unuran" generator object implements an inversion method}

\description{
  Test whether a given \code{unuran} generator object implements an
  inversion method.
  
  [Advanced] -- Test type of method.
}

\usage{
unuran.is.inversion(unr)
}

\arguments{
  \item{unr}{a \code{unuran} object.}
}

\details{
  A \code{unuran} object may implement quite a couple of generation
  methods. This function tests whether the method used in \code{unr}
  is an (approximate) inversion method.

  It returns \code{TRUE} when \code{unr} implements an inversion method
  and \code{FALSE} otherwise.
}

%% \value{}

\seealso{%
  \code{\linkS4class{unuran}}.
}


\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## PINV is an inversion method
unr <- pinvd.new(udnorm())
unuran.is.inversion(unr)

## TDR is a rejection method
unr <- tdrd.new(udnorm())
unuran.is.inversion(unr)

}

\keyword{datagen}

