\name{udlaplace}
\alias{udlaplace}

\title{UNU.RAN object for Laplace distribution}

\description{
  Create UNU.RAN object for a Laplace (double exponential) distribution
  with location parameter \code{location} and scale
  parameter \code{scale}.
  
  [Distribution] -- Laplace (double exponential).  %% Continuous
}

\usage{
udlaplace(location=0, scale=1, lb = -Inf, ub = Inf)
}

\arguments{
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Laplace distribution with location \eqn{l} and scale \eqn{s} has
  density
  \deqn{
    f(x) = \exp( -\frac{|x-l|}{s} )
  }{f(x) = exp( -|x-l|/s )}
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
  Chap. 24, p. 164. 
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for standard Laplace distribution
distr <- udlaplace()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
