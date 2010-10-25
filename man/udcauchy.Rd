\name{udcauchy}
\alias{udcauchy}

\title{UNU.RAN object for Cauchy distribution}

\description{
  Create UNU.RAN object for a Cauchy distribution with
  location parameter \code{location} and scale parameter \code{scale}.

  [Distribution] -- Cauchy.  %% Continuous
}

\usage{
udcauchy(location=0, scale=1, lb=-Inf, ub=Inf)
}

\arguments{
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Cauchy distribution with location \eqn{l} and scale \eqn{s} has
  density
  \deqn{f(x) = \frac{1}{\pi s}
    \left( 1 + \left(\frac{x - l}{s}\right)^2 \right)^{-1}
  }{f(x) = 1 / (pi s (1 + ((x-l)/s)^2))}
  for all \eqn{x}.

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
  Chap. 16, p. 299.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Cauchy distribution
distr <- udcauchy()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
