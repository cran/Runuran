\name{urpowerexp}
\docType{function}

\alias{urpowerexp}

\title{UNU.RAN Powerexponential random variate generator}

\description{
  UNU.RAN random variate generator for the Powerexponential (Subbotin)
  distribution with shape parameter \code{shape}.
  It also allows sampling from the truncated distribution.
  
  [Distribution] -- Sampling Function: Powerexponential (Subbotin).  %% Continuous
}

\usage{
urpowerexp(n, shape, lb = -Inf, ub = Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{shape}{(strictly positive) shape parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Powerexponential distribution with parameter \code{shape}
  \eqn{=\tau}{= t} has density
  \deqn{
    f(x) = \frac{1}{2\Gamma(1+1/\tau)} \exp(-|x|^\tau)}{%
    f(x) = 1/(2 * Gamma(1+1/t)) * exp(-|x|^t )}
  for all \eqn{x} and \eqn{\tau > 0}{t > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented by \R's
  \code{\link{gamma}()} and defined in its help.)

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Powerexponential distribution truncated to the interval (\code{lb},\code{ub}).
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
  Chap.24, p.195. 
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
x <- urpowerexp(n=1000, shape=4)
}

\keyword{distribution}
\keyword{datagen}

