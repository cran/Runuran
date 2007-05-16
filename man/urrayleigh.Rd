\name{urrayleigh}
\docType{function}

\alias{urrayleigh}

\title{UNU.RAN Rayleigh random variate generator}

\description{
  UNU.RAN random variate generator for the Rayleigh distribution with
  scale parameter \code{scale}.
  It also allows sampling from the truncated distribution.
}

\usage{
urrayleigh(n, scale=1, lb = 0, ub = Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{scale} is omitted, it assumes the default value of \code{1}.

  The Rayleigh distribution with scale parameter \code{scale}
  \eqn{=\sigma}{= s} has density 
  \deqn{
    f(x) = \frac{1}{sigma^2} x \exp( -\frac{1}{2}(\frac{x}{sigma})^2 )}{%
    f(x) = 1/s^2 * x * exp(-1/2*(x/s)^2)}
  for \eqn{x \ge 0} and \eqn{\sigma > 0}{s > 0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Rayleigh distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions,
  Volume 1. 2nd edition, 
  John Wiley & Sons, Inc., New York.
  Chap.18, p.456.  
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
}

\examples{
## Create a sample of size 1000
x <- urrayleigh(n=1000,scale=1)
}

\keyword{distribution}
\keyword{datagen}

