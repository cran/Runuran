\name{udlogis}
\alias{udlogis}

\title{UNU.RAN object for Logistic distribution}

\description{
  Create UNU.RAN object for a Logistic distribution with
  parameters \code{location} and \code{scale}.
  
  [Distribution] -- Logistic.  %% Continuous
}

\usage{
udlogis(location=0, scale=1, lb=-Inf, ub=Inf)
}

\arguments{
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Logistic distribution with \code{location} \eqn{= \mu}{= m} and
  \code{scale} \eqn{= \sigma}{= s} has distribution function
  \deqn{
    F(x) = \frac{1}{1 + e^{-(x-\mu)/\sigma}}
  }{F(x) = 1 / (1 + exp(-(x-m)/s))}
  and density
  \deqn{
    f(x) = \frac{1}{\sigma}\frac{e^{(x-\mu)/\sigma}}{(1 + e^{(x-\mu)/\sigma})^2}
  }{f(x) = 1/s exp((x-m)/s) (1 + exp((x-m)/s))^-2.}

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
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 23, p. 115.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for standard logistic distribution
distr <- udlogis()
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
