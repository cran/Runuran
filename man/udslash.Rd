\name{udslash}
\alias{udslash}

\title{UNU.RAN object for Slash distribution}

\description{
  Create UNU.RAN object for a Slash distribution.
  
  [Distribution] -- Slash.  %% Continuous
}

\usage{
udslash(lb=-Inf, ub=Inf)
}

\arguments{
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The slash distribution has density
  \deqn{
    f(x) = \frac{1}{\sqrt{2 \pi}} (1 - exp(-x^2/2)) / x^2
  }{f(x) = 1/sqrt(2*pi) (1 - exp(-x^2/2)) / x^2}
  for \eqn{x\not= 0}{x != 0} and
  \eqn{\frac{1}{\sqrt{2 \pi}}}{1/sqrt(2*pi)} otherwise.
  It is the distribution of the ratio of a unit normal variable to an
  independent standard uniform (0,1) variable.

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
  Chap. 12, p. 63.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for a slash distribution
distr <- udslash()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
