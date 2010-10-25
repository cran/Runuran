\name{udgumbel}
\alias{udgumbel}

\title{UNU.RAN object for Gumbel distribution}

\description{
  Create UNU.RAN object for a Gumbel (Extreme value type I) distribution
  location parameter \code{location} and scale parameter \code{scale}.
  
  [Distribution] -- Gumbel (Extreme value type I).  %% Continuous
}

\usage{
udgumbel(location=0, scale=1, lb=-Inf, ub=Inf)
}

\arguments{
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Gumbel distribution function with location \eqn{l} and
  scale \eqn{s} is
  \deqn{
    F(x) = \exp(-\exp(-\frac{x-l}{s}))
  }{F(x) = exp(-exp(-(x-l)/s))}
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
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 22, p. 2. 
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Gumbel distribution
distr <- udgumbel()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
