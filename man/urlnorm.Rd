\name{urlnorm}
\docType{function}

\alias{urlnorm}

\title{UNU.RAN Log-Normal random variate generator}

\description{
  UNU.RAN random variate generator for the Log-Normal distribution
  whose logarithm has mean equal to \code{meanlog} and standard
  deviation equal to \code{sdlog}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Log-Normal.  %% Continuous
}

\usage{
urlnorm(n, meanlog=0, sdlog=1, lb=0, ub=Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{meanlog, sdlog}{mean and standard deviation of the distribution
    on the log scale. If not not specified they assume the default
    values of \code{0} and \code{1}, respectively.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Log-Normal distribution has density
  \deqn{
    f(x) = \frac{1}{\sqrt{2\pi}\sigma x} e^{-(\log(x) - \mu)^2/2 \sigma^2}%
  }{f(x) = 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))}
  where \eqn{\mu} and \eqn{\sigma} are the mean and standard
  deviation of the logarithm.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Log-Normal distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rlnorm}} for the \R built-in generator.
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
  Compared to \code{rlnorm}, \code{urlnorm} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rlnorm} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urlnorm(n=1000)
}

\keyword{distribution}
