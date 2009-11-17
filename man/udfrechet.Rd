\name{udfrechet}
\docType{function}

\alias{udfrechet}

\title{UNU.RAN object for Frechet distribution}

\description{
  Create UNU.RAN object for a Frechet (Extreme value type II) distribution
  with shape parameter \code{shape}, location parameter \code{location}
  and scale parameter \code{scale}.

  [Distribution] -- Frechet (Extreme value type II).  %% Continuous
}

\usage{
udfrechet(shape, location=0, scale=1, lb=location, ub=Inf)
}

\arguments{
  \item{shape}{(strictly positive) shape parameter.}
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Frechet distribution function with shape \eqn{k}, location \eqn{l}
  and scale \eqn{s} is
  \deqn{
    F(x) = \exp(-(\frac{x-l}{s})^{-k})
  }{F(x) = exp(-((x-l)/s)^(-k))}
  for \eqn{x \ge l}.
  
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
  Chap. 22, p. 2. 
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is a wrapper for the UNU.RAN class in \R.
}

\examples{
## Create distribution object for Frechet distribution
distr <- udfrechet(shape=2)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
