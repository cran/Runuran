\name{udlogarithmic}
\docType{function}

\alias{udlogarithmic}

\title{UNU.RAN object for Logarithmic distribution}

\description{
  Create UNU.RAN object for a Logarithmic distribution with
  shape parameter \code{shape}.
  
  [Distribution] -- Logarithmic.  %% Discrete
}

\usage{
udlogarithmic(shape, lb = 1, ub = Inf)
}

\arguments{
  \item{shape}{shape parameter. Must be between 0 and 1.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Logarithmic distribution with parameters \code{shape}
  \eqn{= \theta}{= t} has density
  \deqn{
    f(x) = -\log(1-\theta) \theta^x / x
  }{f(x) = -log(1-t) * t^x / x}
  for \eqn{x = 1, 2, \ldots} and \eqn{0 < \theta < 1}{0 < t < 1}.

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
  Chap. 7, p. 285.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Logarithmic distribution
dist <- udlogarithmic(shape=0.3)
## Generate generator object; use method DARI
gen <- darid.new(dist)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
