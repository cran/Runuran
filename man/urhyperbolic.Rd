\name{urhyperbolic}
\alias{urhyperbolic}

\title{UNU.RAN Hyperbolic random variate generator}

\description{
  UNU.RAN random variate generator for the Hyperbolic distribution
  with parameters \code{shape} and \code{scale}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Hyperbolic.  %% Continuous
}

\usage{
urhyperbolic(n, shape, scale=1, lb = -Inf, ub = Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{shape}{(strictly positive) shape parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{scale} is omitted, it assumes the default value of \code{1}.

  The Hyperbolic distribution with parameters \code{shape} \eqn{=\alpha}{= a}
  and \code{scale} \eqn{=\sigma}{= s} has density proportional to
  \deqn{
    f(x) \sim \exp(-\alpha \sqrt{1+(\frac{x}{s})^2})
  }{f(x) ~ exp(-a * sqrt(1+(x/s)^2))}
  for all \eqn{x}, \eqn{\alpha > 0}{a > 0} and \eqn{\sigma > 0}{s > 0}.

  The generation algorithm uses transformed density rejection \sQuote{TDR}. The
  parameters \code{lb} and \code{ub} can be used to generate variates from 
  the Hyperbolic distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.

  Do not confuse with \code{\link{rhyper}}
  %% or \code{\link{urhyper}}
  that samples from the (discrete) \emph{hypergeometric} distribution.
}

\examples{
## Create a sample of size 1000 from Hyperbolic distribution with shape=3
x <- urhyperbolic(n=1000,shape=3)
}

\keyword{distribution}
