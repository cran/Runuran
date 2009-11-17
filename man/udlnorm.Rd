\name{udlnorm}
\docType{function}

\alias{udlnorm}

\title{UNU.RAN object for Log Normal distribution}

\description{
  Create UNU.RAN object for a Log Normal distribution
  whose logarithm has mean equal to \code{meanlog} and standard
  deviation equal to \code{sdlog}.
  
  [Distribution] -- Log Normal.  %% Continuous
}

\usage{
udlnorm(meanlog=0, sdlog=1, lb=0, ub=Inf)
}

\arguments{
  \item{meanlog}{mean of the distribution on the log scale.}
  \item{sdlog}{standard deviation of the distribution on the log scale.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The log normal distribution has density
  \deqn{
    f(x) =
    \frac{1}{\sqrt{2 \pi} \sigma x} \exp{- (\log(x)-\mu)^2 / (2 sigma^2)}
  }{f(x) = 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))}
  where \eqn{\mu}{mu} is the mean and \eqn{\sigma}{sigma} the standard
  deviation of the logarithm.
  
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
  Chap. 14, p. 207.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for log normal distribution
distr <- udlnorm()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
