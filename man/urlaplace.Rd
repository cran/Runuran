\name{urlaplace}
\alias{urlaplace}

\title{UNU.RAN Laplace random variate generator}

\description{
  UNU.RAN random variate generator for the Laplace (double exponential)
  distribution with location parameter \code{location} and scale
  parameter \code{scale}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Laplace.  %% Continuous
}

\usage{
urlaplace(n, location=0, scale=1, lb = -Inf, ub = Inf)
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

  The Laplace distribution with location \eqn{l} and scale \eqn{s} has
  density
  \deqn{
    f(x) = \exp( -\frac{|x-l|}{s} )
  }{f(x) = exp( -|x-l|/s )}
  for all \eqn{x}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Laplace distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap.24, p.164. 
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
x <- urlaplace(n=1000)
}

\keyword{distribution}
