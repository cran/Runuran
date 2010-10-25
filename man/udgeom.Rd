\name{udgeom}
\alias{udgeom}

\title{UNU.RAN object for Geometric distribution}

\description{
  Create UNU.RAN object for a Geometric distribution
  with parameter \code{prob}.
  
  [Distribution] -- Geometric.  %% Discrete
}

\usage{
udgeom(prob, lb = 0, ub = Inf)
}

\arguments{
  \item{prob}{probability of success in each trial. \code{0 < prob <= 1}.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Geometric distribution with \code{prob} \eqn{= p} has density
  \deqn{p(x) = p (1-p)^x}
  for \eqn{x = 0, 1, 2, \ldots}{x = 0, 1, 2, \dots}, \eqn{0 < p \le 1}.

  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.discr"}.
}

\references{
  N.L. Johnson, S. Kotz, and A.W. Kemp (1992):
  Univariate Discrete Distributions.
  2nd edition, John Wiley & Sons, Inc., New York.
  Sect. 5.2, p. 201
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Geometric distribution
dist <- udgeom(prob=0.33)
## Generate generator object; use method DARI
gen <- darid.new(dist)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
