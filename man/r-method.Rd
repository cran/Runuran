\name{r-methods}
\docType{methods}

\alias{r}
\alias{r,unuran-method}
\alias{r,unuran,numeric-method}

\title{Sample from a distribution specified by a "unuran" object}
\description{
  Get random sample from a \code{unuran} object in package
  \pkg{Runuran}.
  
  [Universal] -- Sampling Function.
}

\usage{
\S4method{r}{unuran,numeric}(unur,n=1)
}
\arguments{
  \item{unur}{a \code{unuran} object.}
  \item{n}{sample size.}
}

\seealso{%
  \code{\linkS4class{unuran}}.
}

\note{
  This S4 method is equivalent to the function
  \code{\link{ur}}. However, it is (slightly) slower due to some
  overhead for the S4 mechanism.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

%%\examples{ }

\keyword{methods}
\keyword{distribution}
\keyword{datagen}

