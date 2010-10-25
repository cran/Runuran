\name{urlogarithmic}
\alias{urlogarithmic}

\title{UNU.RAN Logarithmic random variate generator}

\description{
  UNU.RAN random variate generator for the Logarithmic distribution with
  shape parameter \code{shape}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Logarithmic.  %% Discrete
}

\usage{
urlogarithmic(n, shape, lb = 1, ub = Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{shape}{shape parameter. Must be between 0 and 1.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Logarithmic distribution with parameters \code{shape} \eqn{=
    \theta}{= t} has density
  \deqn{
    f(x) = -\log(1-\theta) \theta^x / x
  }{f(x) = -log(1-t) * t^x / x}
  for \eqn{x = 1, 2, \ldots} and \eqn{0 < \theta < 1}{0 < t < 1}.

  The generation algorithm uses guide table based inversion when the
  tails are not too heavy and method \sQuote{DARI} otherwise.
  The parameters \code{lb} and \code{ub} can be used to generate
  variates from the Logarithmic distribution truncated to the interval
  (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.

  N.L. Johnson, S. Kotz, and A.W. Kemp (1992):
  Univariate Discrete Distributions, 2nd edition.
  John Wiley & Sons, Inc., New York.
  Chap.7, p.285.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is a wrapper for the UNU.RAN class in \R.
}

\examples{
## Create a sample of size 1000
x <- urlogarithmic(n=1000,shape=0.3)
}

\keyword{distribution}

