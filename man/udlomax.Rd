\name{udlomax}
\docType{function}

\alias{udlomax}

\title{UNU.RAN object for Lomax distribution}

\description{
  Create UNU.RAN object for a Lomax distribution (Pareto
  distribution of second kind) with shape parameter \code{shape} and
  scale parameter \code{scale}.
  
  [Distribution] -- Lomax (Pareto of second kind).  %% Continuous
}

\usage{
udlomax(shape, scale=1, lb=0, ub=Inf)
}

\arguments{
  \item{shape}{(strictly positive) shape parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Lomax distribution with parameters \code{shape} \eqn{=\alpha}{=a}
  and \code{scale} \eqn{=\sigma}{= s} has density
  \deqn{
    f(x) = \alpha \sigma^\alpha (x+\sigma)^{-(\alpha+1)}
  }{f(x) = a s^a (x+s)^(-a-1)}
  for \eqn{x \ge 0}, \eqn{\alpha > 0}{a > 0} and \eqn{\sigma > 0}{s > 0}.

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
  Chap. 20, p. 575.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Lomax distribution
distr <- udlomax(shape=2)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
