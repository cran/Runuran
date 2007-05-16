\name{urexp}
\docType{function}

\alias{urexp}

\title{UNU.RAN Exponential random variate generator}

\description{
  UNU.RAN random variate generator for the Exponential distribution with
  rate \code{rate} (i.e., mean \code{1/rate}).
  It also allows sampling from the truncated distribution.
}

\usage{
urexp(n, rate=1, lb=0, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{rate}{(strictly positive) rate parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{rate} is not specified, it assumes the default value of
  \code{1}.

  The Exponential distribution with rate \eqn{\lambda} has density
  \deqn{
    f(x) = \lambda {e}^{- \lambda x}}{
    f(x) = lambda e^(- lambda x)}
  for \eqn{x \ge 0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Exponential distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\link{unuran-class}} for the UNU.RAN class, and
  \code{\link{rexp}} for the \R built-in generator.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
  Compared to \code{rexp}, \code{urexp} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rexp} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urexp(n=1000)
}

\keyword{distribution}
\keyword{datagen}

