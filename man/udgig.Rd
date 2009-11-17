\name{udgig}
\docType{function}

\alias{udgig}

\title{UNU.RAN object for Generalized Inverse Gaussian distribution}

\description{
  Create UNU.RAN object for a Generalized Inverse Gaussian distribution
  with parameters \code{theta}, \code{psi}, and \code{chi}.
  
  [Distribution] -- Generalized Inverse Gaussian.  %% Continuous
}

\usage{
udgig(theta, psi, chi, lb=0, ub=Inf)
}

\arguments{
  \item{theta}{shape parameter.}
  \item{psi}{shape parameter (must be strictly positive).}
  \item{chi}{shape parameter (must be strictly positive).}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The generalized inverse Gaussian distribution with parameters
  \eqn{\theta}{theta}, \eqn{\psi}{psi}, and \eqn{\chi}{chi}
  has density proportional to
  \deqn{
    f(x) = x^{\theta-1} \exp( -\frac{1}{2} (\frac{\chi}{x} + \psi x))
  }{f(x) = x^(theta-1) * exp( -1/2 * (chi/x + psi*x) )}
  where \eqn{\psi>0}{psi>0} and \eqn{\chi>0}{chi>0}.

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
  Chap. 15, p. 284.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for GIG distribution
distr <- udgig(theta=3, psi=1, chi=1)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
