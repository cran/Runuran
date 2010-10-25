\name{udghyp}
\alias{udghyp}

\title{UNU.RAN object for Generalized Hyperbolic distribution}

\description{
  Create UNU.RAN object for a Generalized Hyperbolic distribution
  with shape parameter \code{lambda}, shape parameter
  \code{alpha}, asymmetry (shape) parameter \code{beta}, scale
  parameter \code{delta}, and location parameter \code{mu}.
  
  [Distribution] -- Generalized Hyperbolic.  %% Continuous
}

\usage{
udghyp(lambda, alpha, beta, delta, mu, lb=-Inf, ub=Inf)
}

\arguments{
  \item{lambda}{shape parameter.}
  \item{alpha}{shape parameter (must be strictly larger than
    absolute value of \code{beta}).}
  \item{beta}{asymmetry (shape) parameter.}
  \item{delta}{scale parameter (must be strictly positive).}
  \item{mu}{location parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The generalized hyperbolic distribution with parameters
  \eqn{\lambda}{lambda}, \eqn{\alpha}{alpha}, \eqn{\beta}{beta},
  \eqn{\delta}{delta}, and \eqn{\mu}{mu}
  has density
  \deqn{
    f(x) = \kappa \;
    (\delta^2+(x-\mu)^2)^{1/2 (\lambda-1/2)}
    \cdot \exp(\beta(x-\mu))
    \cdot K_{\lambda-1/2}\left(\alpha\sqrt{\delta^2+(x-\mu)^2}\right)
  }{f(x) = kappa * (delta^2+(x-mu)^2)^(1/2*(lambda-1/2)) * exp(beta*(x-mu)) * K_{lambda-1/2}(alpha * sqrt(delta^2+(x-mu)^2))}
  where the normalization constant is given by
  \deqn{
    \kappa =
    \frac{\left(\sqrt{\alpha^2 - \beta^2}/\delta\right)^{\lambda}}{
      \sqrt{2\pi} \, \alpha^{\lambda-1/2} \,
      K_{\lambda}\left(\delta \sqrt{\alpha^2-\beta^2}\right)}
  }{kappa = (sqrt(alpha^2 - beta^2)/delta)^lambda / (sqrt(2*pi) * alpha^(lambda-1/2) * K_{lambda}(delta*sqrt(alpha^2-beta^2))}
        
  \eqn{K_{\lambda}(t)}{K_(lambda)(t)} is the modified Bessel function of
  the third kind with index \eqn{\lambda}{lambda}.

  Notice that \eqn{\alpha>|\beta|}{alpha > |beta|} and
  \eqn{\delta>0}{delta>0}.

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
  Barndorff-Nielsen, O., Blaesild, P., 1983. Hyperbolic 
  distributions. In: Johnson, N. L., Kotz, S., Read, C. B. 
  (Eds.), Encyclopedia of Statistical Sciences. Vol. 3. 
  Wiley, New York, p. 700--707.

  Prause, K., 1997. Modelling financial data using generalized 
  hyperbolic distributions. FDM preprint 48, University of Freiburg.

  Prause, K., 1999. The generalized hyperbolic model: Estimation, 
  financial derivatives, and risk measures. Ph.D. thesis, 
  University of Freiburg.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for generalized hyperbolic distribution
distr <- udghyp(lambda=-1.0024, alpha=39.6, beta=4.14, delta=0.0118, mu=-0.000158)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
