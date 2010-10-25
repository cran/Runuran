\name{udweibull}
\alias{udweibull}

\title{UNU.RAN object for Weibull distribution}

\description{
  Create UNU.RAN object for a Weibull (Extreme value type III)
  distribution with with parameters \code{shape} and \code{scale}.
  
  [Distribution] -- Weibull (Extreme value type III).  %% Continuous
}

\usage{
udweibull(shape, scale=1, lb=0, ub=Inf)
}

\arguments{
  \item{shape}{(strictly positive) shape parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Weibull distribution with \code{shape} parameter \eqn{a} and
  \code{scale} parameter \eqn{\sigma}{b} has density given by
  \deqn{
    f(x) = (a/\sigma) {(x/\sigma)}^{a-1} \exp (-{(x/\sigma)}^{a})
  }{f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a)}
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
  Chap. 21, p. 628.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Weibull distribution
distr <- udweibull(shape=3)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
