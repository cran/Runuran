\name{udpowerexp}
\alias{udpowerexp}

\title{UNU.RAN object for Powerexponential distribution}

\description{
  Create UNU.RAN object for a Powerexponential (Subbotin)
  distribution with shape parameter \code{shape}.
  
  [Distribution] -- Powerexponential (Subbotin).  %% Continuous
}

\usage{
udpowerexp(shape, lb=-Inf, ub=Inf)
}

\arguments{
  \item{shape}{(strictly positive) shape parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Powerexponential distribution with parameter \code{shape}
  \eqn{=\tau}{= t} has density
  \deqn{
    f(x) = \frac{1}{2\Gamma(1+1/\tau)} \exp(-|x|^\tau)
  }{f(x) = 1/(2 * Gamma(1+1/t)) * exp(-|x|^t )}
  for all \eqn{x} and \eqn{\tau > 0}{t > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented by \R's
  \code{\link{gamma}()} and defined in its help.)


  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.cont"}.
}

\seealso{
  \code{\linkS4class{unuran.cont}}.
}

\references{
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 24, p. 195. 
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
}

\examples{
## Create distribution object for powerexponential distribution
distr <- udpowerexp(shape=4)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
