\name{unuran.sample-methods}
\docType{methods}

\alias{unuran.sample-methods}
\alias{unuran.sample}
\alias{unuran.sample,unuran-method}
\alias{unuran.sample,unuran,numeric-method}
\alias{r}
\alias{r,unuran-method}
\alias{r,unuran,numeric-method}

\title{Sample from a distribution specified by an UNU.RAN object}
\description{
  Get random sample from an UNU.RAN object in package
  \pkg{Runuran}.
}

\usage{
\S4method{unuran.sample}{unuran,numeric}(unur,n=1)
\S4method{r}{unuran,numeric}(unur,n=1)
}
\arguments{
  \item{unur}{an "unuran" object.}
  \item{n}{sample size.}
}

\seealso{%
  \code{\link{unuran-class}}.
}

\note{
  For performance reasons \code{unuran.sample} is implemented as a
  simple \code{function}.

  \code{r} is an S4 method and equivalent to \code{unuran.sample}.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

%%\examples{ }

\keyword{methods}
\keyword{distribution}
\keyword{datagen}

