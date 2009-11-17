\name{udrayleigh}
\docType{function}

\alias{udrayleigh}

\title{UNU.RAN object for Rayleigh distribution}

\description{
  Create UNU.RAN object for a Rayleigh distribution with
  scale parameter \code{scale}.
  
  [Distribution] -- Rayleigh.  %% Continuous
}

\usage{
udrayleigh(scale=1, lb=0, ub=Inf)
}

\arguments{
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Rayleigh distribution with scale parameter \code{scale}
  \eqn{=\sigma}{= s} has density 
  \deqn{
    f(x) = \frac{1}{sigma^2} x \exp( -\frac{1}{2}(\frac{x}{sigma})^2 )
  }{f(x) = 1/s^2 * x * exp(-1/2*(x/s)^2)}
  for \eqn{x \ge 0} and \eqn{\sigma > 0}{s > 0}.

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
  Chap. 18, p. 456.  
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for standard Rayleigh distribution
distr <- udrayleigh()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
