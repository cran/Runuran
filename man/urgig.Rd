\name{urgig}
\docType{function}

\alias{urgig}

\title{UNU.RAN Generalized Inverse Gaussian Distribution variate generator}

\description{
  UNU.RAN random variate generator for the Generalized Inverse Gaussian
  Distribution with parameters \code{lambda} and \code{omega}.
  It also allows sampling from the truncated distribution.
  
  [Distribution] -- Sampling Function: GIG (generalized inverse Gaussian).  %% Continuous
}

\usage{
urgig(n, lambda, omega, lb=1.e-12, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{lambda}{(strictly positive) shape parameter.}
  \item{omega}{(strictly positive) shape parameter.}
  \item{lb}{lower bound of (truncated) distribution}
  \item{ub}{upper bound of (truncated) distribution}
}

\details{
  The Generalized Inverse Gaussian distribution with parameters 
  \code{lambda} \eqn{=\lambda}{=a} and \code{omega} \eqn{=\omega}{=b}
  has a density proportional to 
  \deqn{
    f(x)\sim x^{\lambda-1}\exp(-(\omega/2)(x+1/x))}{%
    f(x) ~ x^(a-1) exp(-(b/2)(x+1/x))}
  for \eqn{x \ge 0}, \eqn{\lambda > 0}{a > 0} and \eqn{\omega > 0}{b > 0}.

  The generation algorithm uses transformed density rejection "TDR". The
  parameters \code{lb} and \code{ub} can be used to generate variates from 
  the distribution truncated to the interval (\code{lb},\code{ub}).

  The generation algorithm works for
  \eqn{\lambda \ge 1} and \eqn{\omega>0} and 
  for \eqn{\lambda>0} and \eqn{\omega \ge 0.5}. 
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.

  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions,
  Volume 1. 2nd edition, 
  John Wiley & Sons, Inc., New York.
  Chap.15, p.284.  
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
}

\examples{
## Create a sample of size 1000
x <- urgig(n=1000,lambda=2,omega=3)
}

\keyword{distribution}
\keyword{datagen}

