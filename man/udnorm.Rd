\name{udnorm}
\docType{function}

\alias{udnorm}

\title{UNU.RAN object for Normal distribution}

\description{
  Create UNU.RAN object for a Normal (Gaussian) distribution with mean
  equal to \code{mean} and standard deviation to \code{sd}.
  
  [Distribution] -- Normal (Gaussian).  %% Continuous
}

\usage{
udnorm(mean=0, sd=1, lb=-Inf, ub=Inf)
}

\arguments{
  \item{mean}{mean of distribution.}
  \item{sd}{standard deviation.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The normal distribution with mean \eqn{\mu}{mu} and standard deviation
  \eqn{\sigma}{sigma} has density
  \deqn{
    f(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-(x-\mu)^2/2\sigma^2}
  }{f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))}
  where \eqn{\mu}{mu} is the mean of the distribution and
  \eqn{\sigma}{sigma} the standard deviation.

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
  Chap. 13, p. 80.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for standard normal distribution
distr <- udnorm()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

## Create distribution object for positive normal distribution
distr <- udnorm(lb=0, ub=Inf)
## ... and draw a sample
gen <- pinvd.new(distr)
x <- ur(gen,100)

}

\keyword{distribution}
