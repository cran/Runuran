\name{urlogis}
\alias{urlogis}

\title{UNU.RAN Logistic random variate generator}

\description{
  UNU.RAN random variate generator for the Logistic distribution with
  parameters \code{location} and \code{scale}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Logistic.  %% Continuous
}

\usage{
urlogis(n, location=0, scale=1, lb=-Inf, ub=Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{location}{location parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  If \code{location} or \code{scale} are omitted, they assume the
  default values of \code{0} and \code{1} respectively.

  The Logistic distribution with \code{location} \eqn{= \mu}{= m} and
  \code{scale} \eqn{= \sigma}{= s} has distribution function
  \deqn{
    F(x) = \frac{1}{1 + e^{-(x-\mu)/\sigma}}%
  }{F(x) = 1 / (1 + exp(-(x-m)/s))}
  and density
  \deqn{
    f(x)= \frac{1}{\sigma}\frac{e^{(x-\mu)/\sigma}}{(1 + e^{(x-\mu)/\sigma})^2}%
  }{f(x) = 1/s exp((x-m)/s) (1 + exp((x-m)/s))^-2}

  The generation algorithm uses inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Logistic distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rlogis}} for the \R built-in generator.
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
  Compared to \code{rlogis}, \code{urlogis} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rlogis} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urlogis(n=1000)
}

\keyword{distribution}
