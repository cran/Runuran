\name{urweibull}
\docType{function}

\alias{urweibull}

\title{UNU.RAN Weibull random variate generator}

\description{
  UNU.RAN random variate generator for the Weibull distribution with
  with parameters \code{shape} and \code{scale}.
  It also allows sampling from the truncated distribution.
}

\usage{
urweibull(n, shape, scale=1, lb=0, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{shape}{(strictly positive) shape parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Weibull distribution with \code{shape} parameter \eqn{a} and
  \code{scale} parameter \eqn{\sigma}{b} has density given by
  \deqn{f(x) = (a/\sigma) {(x/\sigma)}^{a-1} \exp (-{(x/\sigma)}^{a})}{%
        f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a)}
  for \eqn{x \ge 0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Weibull distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\link{unuran-class}} for the UNU.RAN class, and
  \code{\link{rweibull}} for the \R built-in generator.
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
  Compared to \code{rweibull}, \code{urweibull} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rweibull} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urweibull(n=1000,shape=3)
}

\keyword{distribution}
\keyword{datagen}

