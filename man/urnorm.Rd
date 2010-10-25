\name{urnorm}
\alias{urnorm}

\title{UNU.RAN Normal random variate generator}

\description{
  UNU.RAN random variate generator for the Normal distribution with mean
  equal to \code{mean} and standard deviation to \code{sd}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Normal (Gaussian).  %% Continuous
}

\usage{
urnorm(n, mean = 0, sd = 1, lb = -Inf, ub = Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{mean}{mean of distribution.}
  \item{sd}{standard deviation.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{mean} or \code{sd} are not specified they assume the default
  values of \code{0} and \code{1}, respectively.

  The normal distribution has density
  \deqn{
    f(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-(x-\mu)^2/2\sigma^2}
  }{f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))}
  where \eqn{\mu}{mu} is the mean of the distribution and
  \eqn{\sigma}{sigma} the standard deviation.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Normal distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rnorm}} for the \R built-in normal random variate
  generator.
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
  This function is a wrapper for the UNU.RAN class in \R.
  Compared to \code{rnorm}, \code{urnorm} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rnorm} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urnorm(n=1000)
}

\keyword{distribution}
