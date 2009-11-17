\name{ur}
\docType{function}

\alias{ur}
\alias{unuran.sample}

\title{Sample from a distribution specified by a "unuran" object}

\description{
  Get random sample from a \code{unuran} object in package
  \pkg{Runuran}.
  
  [Universal] -- Sampling Function.
}

\usage{
ur(unr, n=1)
unuran.sample(unr, n=1)
}

\arguments{
  \item{unr}{a \code{unuran} object.}
  \item{n}{sample size.}
}

\seealso{%
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\note{
  \code{unuran.sample} is just an (older) longer name for \code{ur}.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Draw random sample of size 10 from normal distribution using
## method 'TDR'
unr <- unuran.new("normal","tdr")
x <- ur(unr,n=10)
}

\keyword{distribution}
\keyword{datagen}

