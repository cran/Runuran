\name{urextremeII}
\docType{function}

\alias{urextremeII}

\title{UNU.RAN Extreme value type II (Frechet-type) random variate generator}

\description{
  UNU.RAN random variate generator for the Extreme value type II
  (Frechet-type) distribution with shape parameter \code{shape},
  location parameter \code{location} and scale parameter \code{scale}.
  It also allows sampling from the truncated distribution.
}

\usage{
urextremeII(n, shape, location=0, scale=1, lb=location, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{shape}{(strictly positive) shape parameter.}
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{location} or \code{scale} are not specified, they assume
  the default values of \code{0} and \code{1} respectively.

  The Frechet distribution with shape \eqn{k}, location \eqn{l} and
  scale \eqn{s} has density
  \deqn{
    f(x) = \exp(-(\frac{x-l}{s})^{-k})\,(\frac{x-l}{s})^{-k-1}}{%
    f(x) = exp(-((x-l)/s)^(-k))*((x-l)/s)^(-k-1)}
  for \eqn{x \ge l}.
  
  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Frechet distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\link{unuran-class}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.

  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap.22, p.2. 
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
x <- urextremeII(n=1000,shape=2)
}

\keyword{distribution}
\keyword{datagen}

