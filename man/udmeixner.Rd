\name{udmeixner}
\alias{udmeixner}

\title{UNU.RAN object for Meixner distribution}

\description{
  Create UNU.RAN object for a Meixner distribution
  with scale parameter \code{alpha},
  asymmetry (shape) parameter \code{beta},
  shape parameter \code{delta} and
  location parameter \code{mu}.
  
  [Distribution] -- Meixner.  %% Continuous
}

\usage{
udmeixner(alpha, beta, delta, mu, lb=-Inf, ub=Inf)
}

\arguments{
  \item{alpha}{scale parameter (must be strictly positive).}
  \item{beta}{asymmetry (shape) parameter
    (must be larger than \eqn{-\pi}{-pi} and smaller than \eqn{\pi}{pi}).}
  \item{delta}{shape parameter (must be strictly positive).}
  \item{mu}{location parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Mexiner distribution with parameters
 \eqn{\alpha}{alpha}, \eqn{\beta}{beta}, \eqn{\delta}{delta}, and \eqn{\mu}{mu}
  has density 
  \deqn{
    f(x) = \kappa \,\exp(\beta(x-\mu)/\alpha)
    \, |\Gamma\left(\delta+ i(x-\mu)/\alpha\right)|^2
  }{f(x) = kappa * exp(beta*(x-mu)/alpha) * |Gamma(delta + i * (x-mu)/alpha)|^2}
  where the normalization constant is given by
  \deqn{
    \kappa =
    \frac{\left(2\cos\left(\beta/2\right)\right)^{2\delta}}{
      2 \alpha \pi \, \Gamma\left(2 \delta\right)}
  }{kappa = (2*cos(beta/2))^(2*delta) / (2 * alpha * pi * Gamma(2*delta))}

  The symbol \eqn{i} denotes the imaginary unit, that is, we have to
  evaluate the gamma function \eqn{\Gamma(z)}{Gamma(z)} for complex
  arguments \eqn{z= x + i\,y}{z = x + i*y}.

  Notice that \eqn{\alpha>0}{alpha>0}, \eqn{|\beta| < \pi}{|beta| < pi}
  and \eqn{\delta>0}{delta>0}.

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
  Grigelionis, B., 1999. Processes of Meixner type. Lithuanian 
  Mathematical Journal, Vol. 39, p. 33--41.
	  
  Schoutens, W., 2001. The Meixner Processes in Finance. Eurandom 
  Report 2001-002, Eurandom, Eindhoven.
}

\author{
  Josef Leydold and Kemal Dingec
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for meixner distribution
distr <- udmeixner(alpha=0.0298, beta=0.1271, delta=0.5729, mu=-0.0011)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
