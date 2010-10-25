\name{udexp}
\alias{udexp}

\title{UNU.RAN object for Exponential distribution}

\description{
  Create UNU.RAN object for an Exponential distribution with 
  rate \code{rate} (i.e., mean \code{1/rate}).
  
  [Distribution] -- Exponential.  %% Continuous
}

\usage{
udexp(rate=1, lb=0, ub=Inf)
}

\arguments{
  \item{rate}{(strictly positive) rate parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Exponential distribution with rate \eqn{\lambda} has density
  \deqn{
    f(x) = \lambda {e}^{- \lambda x}
  }{f(x) = lambda exp(-lambda x)}
  for \eqn{x \ge 0}.

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
  Chap. 19, p. 494.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for standard exponential distribution
distr <- udexp()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
