\name{udpois}
\alias{udpois}

\title{UNU.RAN object for Poisson distribution}

\description{
  Create UNU.RAN object for a Poisson distribution
  with parameter \code{lambda}.
  
  [Distribution] -- Poisson.  %% Discrete
}

\usage{
udpois(lambda, lb = 0, ub = Inf)
}

\arguments{
  \item{lambda}{(non-negative) mean.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Poisson distribution has density
  \deqn{
    p(x) = \frac{\lambda^x e^{-\lambda}}{x!}
  }{p(x) = lambda^x exp(-lambda)/x!}
  for \eqn{x = 0, 1, 2, \ldots}.

  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.discr"}.
}

\seealso{
  \code{\linkS4class{unuran.discr}}.
}

\references{
  N.L. Johnson, S. Kotz, and A.W. Kemp (1992):
  Univariate Discrete Distributions.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 4, p. 151.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Poisson distribution
dist <- udpois(lambda=2.3)
## Generate generator object; use method DARI
gen <- darid.new(dist)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
