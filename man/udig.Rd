\name{udig}
\alias{udig}

\title{UNU.RAN object for Inverse Gaussian distribution}

\description{
  Create UNU.RAN object for a Inverse Gaussian (Wald) distribution
  with mean \code{mu} and shape parameter \code{lambda}.
  
  [Distribution] -- Inverse Gaussian (Wald).  %% Continuous
}

\usage{
udig(mu, lambda, lb=0, ub=Inf)
}

\arguments{
  \item{mu}{mean (strictly positive).}
  \item{lambda}{shape parameter (strictly positive).}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The inverse Gaussian distribution with mean \eqn{\mu}{mu} and shape
  parameter \eqn{\lambda}{lambda}
  has density
  \deqn{
    f(x) =
    \sqrt{\frac{\lambda}{2 \pi x^3} }
    \exp( -\frac{\lambda (x-\mu)^2}{2\mu^2 x} )
  }{f(x) = sqrt(lambda / (2*pi*x^3)) * exp( -(lambda*(x-mu)^2) / (2*mu^2*x) )}
  where \eqn{\mu>0}{mu>0} and \eqn{\lambda>0}{lambda>0}.

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
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions, Volume 1.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 15, p. 259.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for inverse Gaussian distribution
distr <- udig(mu=3, lambda=2)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
