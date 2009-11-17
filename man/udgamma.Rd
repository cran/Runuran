\name{udgamma}
\docType{function}

\alias{udgamma}

\title{UNU.RAN object for Gamma distribution}

\description{
  Create UNU.RAN object for a Gamma distribution
  with parameters \code{shape} and \code{scale}.
  
  [Distribution] -- Gamma.  %% Continuous
}

\usage{
udgamma(shape, scale=1, lb=0, ub=Inf)
}

\arguments{
  \item{shape}{(strictly positive) shape parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution}
  \item{ub}{upper bound of (truncated) distribution}
}

\details{
  The Gamma distribution with parameters \code{shape} \eqn{=\alpha}{= a}
  and \code{scale} \eqn{=\sigma}{= s} has density
  \deqn{
    f(x) = \frac{1}{{\sigma}^{\alpha}\Gamma(\alpha)} {x}^{\alpha-1} e^{-x/\sigma}
  }{f(x) = 1/(s^a Gamma(a)) x^(a-1) exp(-x/s)}
  for \eqn{x \ge 0}, \eqn{\alpha > 0}{a > 0} and \eqn{\sigma > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented by \R's
  \code{\link{gamma}()} and defined in its help.)

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
  Chap. 17, p. 337.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for gamma distribution
distr <- udgamma(shape=4)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
