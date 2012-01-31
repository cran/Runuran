\name{udvg}
\alias{udvg}

\title{UNU.RAN object for Variance Gamma distribution}

\description{
  Create UNU.RAN object for a Variance Gamma distribution
  with shape parameter \code{lambda}, shape parameter
  \code{alpha}, asymmetry (shape) parameter \code{beta}, and
  location parameter \code{mu}.
  
  [Distribution] -- Variance Gamma.  %% Continuous
}

\usage{
udvg(lambda, alpha, beta, mu, lb=-Inf, ub=Inf)
}

\arguments{
  \item{lambda}{shape parameter (must be strictly positive).}
  \item{alpha}{shape parameter (must be strictly larger than
    absolute value of \code{beta}).}
  \item{beta}{asymmetry (shape) parameter.}
  \item{mu}{location parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The variance gamma distribution with parameters
  \eqn{\lambda}{lambda}, \eqn{\alpha}{alpha}, \eqn{\beta}{beta}, and \eqn{\mu}{mu}
  has density 
  \deqn{
    f(x) = \kappa \;
    |x-\mu|^{\lambda-1/2}
    \cdot \exp(\beta(x-\mu))
    \cdot K_{\lambda-1/2}\left(\alpha|x-\mu|\right)
  }{f(x) = kappa * |x-mu|^(lambda-1/2) * exp(beta*(x-mu)) * K_{lambda-1/2}(alpha * |x-mu|)}
  where the normalization constant is given by
  \deqn{
    \kappa =
    \frac{\left(\alpha^2 - \beta^2\right)^{\lambda}}{
      \sqrt{\pi} \, (2 \alpha)^{\lambda-1/2} \,
      \Gamma\left(\lambda\right)}
  }{kappa = (alpha^2 - beta^2)^lambda / (sqrt(pi) * (2*alpha)^(lambda-1/2) * Gamma(lambda))}
        
  \eqn{K_{\lambda}(t)}{K_(lambda)(t)} is the modified Bessel function of
  the third kind with index \eqn{\lambda}{lambda}.
  \eqn{\Gamma(t)}{Gamma(t)} is the Gamma function.
    
  Notice that \eqn{\alpha>|\beta|}{alpha > |beta|} and
  \eqn{\lambda>0}{lambda > 0}.

  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.cont"}.
}

\note{
  For \eqn{\lambda \le 0.5}{lambda <= 0.5}, the density has a pole at
  \eqn{\mu}{mu}.
}

\seealso{
  \code{\linkS4class{unuran.cont}}.
}

\references{
  Eberlein, E., von Hammerstein, E. A., 2004. Generalized hyperbolic 
  and inverse Gaussian distributions: limiting cases and approximation 
  of processes. In Seminar on Stochastic Analysis, Random Fields and 
  Applications IV, Progress in Probability 58, R. C. Dalang, M. Dozzi, 
  F. Russo (Eds.), Birkhauser Verlag, p. 221--264.
  
  Madan, D. B., Seneta, E., 1990. The variance gamma (V.G.) model for 
  share market returns. Journal of Business, Vol. 63, p. 511--524.
  
  Raible, S., 2000. L\'evy Processes in Finance: Theory, Numerics, and 
  Empirical Facts. Ph.D. thesis, University of Freiburg.
}

\author{
  Josef Leydold and Kemal Dingec
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for variance gamma distribution
distr <- udvg(lambda=2.25, alpha=210.5, beta=-5.14, mu=0.00094)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)
}

\keyword{distribution}

