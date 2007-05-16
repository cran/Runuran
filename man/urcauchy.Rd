\name{urcauchy}
\docType{function}

\alias{urcauchy}

\title{UNU.RAN Cauchy random variate generator}

\description{
  UNU.RAN random variate generator for the Cauchy distribution with
  location parameter \code{location} and scale parameter \code{scale}.
  It also allows sampling from the truncated distribution.
}

\usage{
urcauchy(n, location=0, scale=1, lb = -Inf, ub = Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{location} or \code{scale} are not specified, they assume
  the default values of \code{0} and \code{1} respectively.

  The Cauchy distribution with location \eqn{l} and scale \eqn{s} has
  density
  \deqn{f(x) = \frac{1}{\pi s}
    \left( 1 + \left(\frac{x - l}{s}\right)^2 \right)^{-1}%
  }{f(x) = 1 / (pi s (1 + ((x-l)/s)^2))}
  for all \eqn{x}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Cauchy distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\link{unuran-class}} for the UNU.RAN class, and
  \code{\link{rcauchy}} for the \R built-in generator.
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
  Compared to \code{rcauchy}, \code{urcauchy} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rcauchy} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urcauchy(n=1000)
}

\keyword{distribution}
\keyword{datagen}

