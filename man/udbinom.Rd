\name{udbinom}
\docType{function}

\alias{udbinom}

\title{UNU.RAN object for Binomial distribution}

\description{
  Create UNU.RAN object for a Binomial distribution
  with parameters \code{size} and \code{prob}.
  
  [Distribution] -- Binomial.  %% Discrete
}

\usage{
udbinom(size, prob, lb=0, ub=size)
}

\arguments{
  \item{size}{number of trials (one or more).}
  \item{prob}{probability of success on each trial.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Binomial distribution with \code{size} \eqn{= n} and
  \code{prob} \eqn{= p} has probability mass function
  \deqn{p(x) = {n \choose x} {p}^{x} {(1-p)}^{n-x}
  }{p(x) = choose(n,x) p^x (1-p)^(n-x)}
  for \eqn{x = 0, \ldots, n}.

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
  Chap. 3, p. 105.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Binomial distribution
dist <- udbinom(size=100, prob=0.33)
## Generate generator object; use method DGT (inversion)
gen <- dgtd.new(dist)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
