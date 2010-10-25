\name{urpois}
\alias{urpois}

\title{UNU.RAN Poisson random variate generator}

\description{
  UNU.RAN random variate generator for the Poisson distribution
  with parameter \code{lambda}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Poisson.  %% Discrete
}

\usage{
urpois(n, lambda, lb = 0, ub = Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{lambda}{(non-negative) mean.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Poisson distribution has density
  \deqn{
    p(x) = \frac{\lambda^x e^{-\lambda}}{x!}
  }{p(x) = lambda^x exp(-lambda)/x!}
  for \eqn{x = 0, 1, 2, \ldots}.

  The generation algorithm uses guide table based inversion when the
  tails are not too heavy and method \sQuote{DARI} otherwise.
  The parameters \code{lb} and \code{ub} can be used to generate
  variates from the Poisson distribution truncated to the interval
  (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rpois}} for the \R built-in generator.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
  Compared to \code{rpois}, \code{urpois} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rpois} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000 from Poisson distribution with lamda=2.3
x <- urpois(n=1000,lambda=2.3)
}

\keyword{distribution}
