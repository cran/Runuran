\name{udnbinom}
\alias{udnbinom}

\title{UNU.RAN object for Negative Binomial distribution}

\description{
  Create UNU.RAN object for a Negative Binomial distribution
  with parameters \code{size} and \code{prob}.
  
  [Distribution] -- Negative Binomial.  %% Discrete
}

\usage{
udnbinom(size, prob, lb = 0, ub = Inf)
}

\arguments{
  \item{size}{target for number of successful trials, or dispersion
    parameter (the shape parameter of the gamma mixing distribution).
    Must be strictly positive.}
  \item{prob}{probability of success in each trial. \code{0 < prob <= 1}.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Negative Binomial distribution with \code{size} \eqn{= n} and
  \code{prob} \eqn{= p} has density
  \deqn{
    p(x) = \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x
  }{p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x}
  for \eqn{x = 0, 1, 2, \ldots}, \eqn{n > 0} and \eqn{0 < p \le 1}.
  This represents the number of failures which occur in a sequence of
  Bernoulli trials before a target number of successes is reached.

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
  Sect. 5.1, p. 200.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Negative Binomial distribution
dist <- udnbinom(size=100, prob=0.33)
## Generate generator object; use method DARI
gen <- darid.new(dist)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
