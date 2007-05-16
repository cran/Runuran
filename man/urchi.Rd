\name{urchi}
\docType{function}

\alias{urchi}

\title{UNU.RAN Chi random variate generator}

\description{
  UNU.RAN random variate generator for the Chi
  (\eqn{\chi}{chi}) distribution with \code{df} degrees of freedom.
  It also allows sampling from the truncated distribution.
  (Do not confuse with the Chi-Squared (\eqn{\chi^2}{chi^2})
  distribution!)
}

\usage{
urchi(n, df, lb=0, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{df}{degrees of freedom (strictly positive, but can be non-integer).}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Chi distribution with \code{df}\eqn{= n > 0} degrees of
  freedom has density
  \deqn{
    f(x) = x^{n-1} e^{-x^2/2}}{%
    f(x) = x^(n-1) * exp(-x^2/2)}
  for \eqn{x > 0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Chi distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions,
  Volume 1. 2nd edition, 
  John Wiley & Sons, Inc., New York.
  Chap.18, p.417.  
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
x <- urchi(n=1000,df=3)
}

\keyword{distribution}
\keyword{datagen}

