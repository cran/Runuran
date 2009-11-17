\name{udhyperbolic}
\docType{function}

\alias{udhyperbolic}

\title{UNU.RAN object for Hyperbolic distribution}

\description{
  Create UNU.RAN object for a Hyperbolic distribution
  with location parameter \code{mu}, tail (shape) parameter
  \code{alpha}, asymmetry (shape) parameter \code{beta}, and scale
  parameter \code{delta}.
  
  [Distribution] -- Hyperbolic.  %% Continuous
}

\usage{
udhyperbolic(alpha, beta, delta, mu, lb=-Inf, ub=Inf)
}

\arguments{
  \item{alpha}{tail (shape) parameter (must be strictly larger than
  absolute value of \code{beta}).}
  \item{beta}{asymmetry (shape) parameter.}
  \item{delta}{scale parameter (must be strictly positive).}
  \item{mu}{location parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The hyperbolic distribution with parameters
  \eqn{\mu}{mu},\eqn{\alpha}{alpha},\eqn{\beta}{beta}, and
  \eqn{\delta}{delta}
  has density proportional to
  \deqn{
    f(x) = \exp( -\alpha \sqrt(\delta^2 + (x - \mu)^2) + \beta*(x-\mu) )
  }{f(x) = exp( -alpha * sqrt(delta^2 + (x - mu)^2) + beta*(x-mu) )}
  where \eqn{\alpha>|\beta|}{alpha > |beta|} and
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

%% \references{}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for hyperbolic distribution
distr <- udhyperbolic(alpha=3,beta=2,delta=1,mu=0)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
